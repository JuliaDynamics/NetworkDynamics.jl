"""
    struct InitConstraint{F}
    InitConstraint(f, sym, dim)

A representation of an additional constraint that is applied during the initialization phase of a component.
It contains a function `f` that defines the constraint, a vector of symbols `sym` that are involved in the constraint,
and the dimension `dim` of the constraint.

    InitConstraint([:x, :y], 2) do res, u
        res[1] = u[:x]^2 + u[:y]^2 - 1
    end

See also [`@initconstraint`](@ref) for a macro to create such constraints.
"""
struct InitConstraint{F}
    f::F
    sym::Vector{Symbol}
    dim::Int
    prettyprint::Union{Nothing,String}
end
InitConstraint(f, sym, dim) = InitConstraint(f, sym, dim, nothing)

"""
    InitConstraint(subconstraints::InitConstraint...)

Combine multiple `InitConstraint` objects into a single constraint.

The resulting constraint will have:
- Combined symbols from all subconstraints (deduplicated)
- Total dimension equal to sum of individual dimensions
- Function that evaluates all subconstraints sequentially

Limitation: All subconstraints need to use `Symbol` indices internally, i.e. u[:x] rather than u[2]!
"""
function InitConstraint(subconstraints::InitConstraint...)
    isempty(subconstraints) && throw(ArgumentError("At least one subconstraint must be provided"))
    if length(subconstraints) == 1
        return only(subconstraints)
    end

    all_syms = unique!(reduce(vcat, c.sym for c in subconstraints))
    total_dim = sum(dim(c) for c in subconstraints)
    f_range = []
    offset = 0
    for c in subconstraints
        tup = (f=c.f, range=(offset+1):(offset+dim(c)))
        offset += dim(c)
        push!(f_range, tup)
    end

    f_range_tup = Tuple(f_range)
    function combined_f(res, u)
        offset = 0
        unrolled_foreach(f_range_tup) do (cf, cr)
            res_view = view(res, cr)
            cf(res_view, u)
        end
        nothing
    end

    try
        su = SymbolicView(zeros(length(all_syms)), all_syms, false)
        res = zeros(total_dim)
        combined_f(res, su)
    catch e
        if e isa IllegalIntIndexingError
            throw(ArgumentError(
                "Cannot combine multiple init constraints, because at least \
                 on of them uses `u[::Int]` indexing internally. Use `u[::Symbol]` \
                 within the `InitConstraint` to fix this error."
            ))
        else
            rethrow(e)
        end
    end

    if any(isnothing, c.prettyprint for c in subconstraints)
        prettyprint = nothing
    else
        # the combined constraint is itself expressible as one `@initconstraint` block:
        # concatenate the subconstraints' bodies, dropping each one's macro header/footer
        header = "@initconstraint begin"
        footer = "end"
        bodylines = mapreduce(vcat, subconstraints) do c
            split(c.prettyprint, '\n')[2:end-1]
        end
        prettyprint = join([header, join(bodylines, "\n"), footer], "\n")
    end

    InitConstraint(combined_f, all_syms, total_dim, prettyprint)
end

dim(c::InitConstraint) = c.dim

(c::InitConstraint)(res, u) = c.f(res, SymbolicView(u, c.sym))

function Base.show(io::IO, ::MIME"text/plain", @nospecialize(c::InitConstraint))
    _show_recipe(io, c)
end

# `prettyprint` is the complete authored recipe — a weak `InitFormula`'s header already carries
# `weak=true` (baked in at construction, see `_formula_macro`/`_build_formula`), so show just
# prints it. Without one (raw constructor) fall back to `repr`.
function _show_recipe(io::IO, @nospecialize(c))
    if isnothing(c.prettyprint)
        # strip trailing default-valued fields so the result stays a valid constructor call:
        # `nothing` (prettyprint/derived_from) and `weak=false`. A non-default `weak=true` stays,
        # keeping the preceding `nothing`s so the positional call still lines up.
        print(io, replace(repr(c), r"(, (nothing|false))+\)$" => ")"))
    else
        print(io, c.prettyprint)
    end
end

"""
   @initconstraint

Generate an [`InitConstraint`](@ref) from an expression using symbols.

    @initconstraint begin
        :x + :y
        :z^2
    end

is equal to

    InitConstraint([:x, :y, :z], 2) do out, u
        out[1] = u[:x] + u[:y]
        out[2] = u[:z]^2
    end
"""
macro initconstraint(ex)
    if ex isa QuoteNode || ex.head != :block
        ex = Base.remove_linenums!(Expr(:block, ex))
    end

    # capture the macro-form source before the symbols get wrapped into `u[:sym]`
    s = _macro_source_string("@initconstraint", ex)

    sym = Symbol[]
    out = :__out__
    u = :__u__
    body = Expr[]

    dim = 0
    for constraint in ex.args
        constraint isa Union{Expr, QuoteNode} || continue # skip line number nodes
        dim += 1
        wrapped = _wrap_symbols!(constraint, sym, u)
        push!(body, :($out[$dim] = $wrapped))
    end
    unique!(sym)

    body_esc = _escape_all.(body)

    quote
        InitConstraint($sym, $dim, $s) do $(esc(out)), $(esc(u))
            $(body_esc...)
            nothing
        end
    end
end
function _wrap_symbols!(ex, sym, u)
    postwalk(ex) do x
        if x isa QuoteNode && x.value isa Symbol
            push!(sym, x.value)
            :($u[$x])
        else
            x
        end
    end
end
function _escape_all(ex::Expr)
    postwalk(ex) do x
        if x isa Symbol
            esc(x)
        else
            x
        end
    end
end

# Reproduce the `@name begin … end` source the user wrote, so a constraint/formula prints
# as its macro call rather than as the expanded inline `do`-block — it reads nicer. Only the
# macros carry this; direct-constructor objects have no source and fall back to `repr` (see
# `_show_recipe`). Must be called before the symbols are wrapped into `u[:sym]`.
function _macro_source_string(macroname, ex)
    lines = [_strip_locations(string("    ", a)) for a in ex.args if a isa Union{Expr, QuoteNode}]
    string(macroname, " begin\n", join(lines, "\n"), "\nend")
end

# Nested macro calls (e.g. `@pf(:x)`) print with a `#= file:line =#` location comment that
# `remove_linenums!` cannot reach; drop it so the source reads cleanly.
_strip_locations(s) = replace(s, r"#=.*?=# " => "")

# Best-effort one-line `label`, only 1 line formulas
function _auto_formula_label(ex)
    exprs = filter(a -> a isa Union{Expr,QuoteNode}, ex.args)
    length(exprs) == 1 || return nothing
    a = only(exprs)
    (a isa Expr && a.head === :(=) && a.args[1] isa QuoteNode) || return nothing
    # `\p{L}` rather than `[A-Za-z]`: symbol names are commonly Greek (`:θ`, `:ω`)
    replace(_strip_locations(string(a)), r":(?=[\p{L}_])" => "")
end

# for metadata check, just passes down the
function assert_initconstraint_compat(cf::ComponentModel, c::InitConstraint)
    allowed_symbols = Set(vcat(
        sym(cf),
        psym(cf),
        insym_flat(cf),
        outsym_flat(cf),
        obssym(cf)
    ))
    missmatch = setdiff(c.sym, allowed_symbols)
    if !isempty(missmatch)
        throw(ArgumentError("InitConstraint requires symbols that are not part of the component model: $missmatch"))
    end
    c
end

"""
prep_initconstrasint allocates all the buffers and so on in order to be called
within the NonlienarProblem during initialization.
"""
function prep_initiconstraint(cm::ComponentModel, c::InitConstraint, chunksize)
    obscache = if !isempty(c.sym ∩ obssym(cm))
        DiffCache(zeros(length(obssym(cm))), chunksize)
    else
        nothing
    end
    symcache = DiffCache(zeros(length(c.sym)), chunksize)
    symtup = Tuple(c.sym)
    obsf! = obsf(cm)

    mapping! = generate_init_input_mapping(cm, c)
    initf! = c.f

    (res, outbufs, ubuf, inbufs, pbuf, t) -> begin
        if !isnothing(obscache)
            obsbuf = PreallocationTools.get_tmp(obscache, res)
            obsf!(obsbuf, ubuf, inbufs..., pbuf, t)
            obsf!
        else
            obsbuf = nothing
        end

        symbuf = PreallocationTools.get_tmp(symcache, res)
        mapping!(symbuf, outbufs, ubuf, inbufs, pbuf, obsbuf)

        symv = SymbolicView(symbuf, symtup)
        initf!(res, symv)

        nothing
    end
end
prep_initiconstraint(cm::ComponentModel, c::Nothing, _) = (args...) -> nothing

"""
Returns a function, which collects all the symbols required for InitConstraint.
For each symbol, it checks if it should be copied from flat outputs, inputs, parameters, or state
buffers.
"""
function generate_init_input_mapping(cm::ComponentModel, c::InitConstraint)
    outmapping = NTuple{3,Int}[]
    umapping   = NTuple{2,Int}[]
    inmapping  = NTuple{3,Int}[]
    pmapping   = NTuple{2,Int}[]
    obsmapping = NTuple{2,Int}[]
    for (iidx, s) in enumerate(c.sym)
        if s in outsym_flat(cm)
            candidates = findfirst.(Ref(isequal(s)), outsym_normalized(cm))
            candidx = findfirst(!isnothing, candidates)
            push!(outmapping, (iidx, candidx, candidates[candidx]))
        elseif s in sym(cm)
            unidx = findfirst(isequal(s), sym(cm))
            push!(umapping, (iidx, unidx))
        elseif s in insym_flat(cm)
            candidates = findfirst.(Ref(isequal(s)), insym_normalized(cm))
            candidx = findfirst(!isnothing, candidates)
            push!(inmapping, (iidx, candidx, candidates[candidx]))
        elseif s in psym(cm)
            pnidx = findfirst(isequal(s), psym(cm))
            push!(pmapping, (iidx, pnidx))
        elseif s in obssym(cm)
            obsidx = findfirst(isequal(s), obssym(cm))
            push!(obsmapping, (iidx, obsidx))
        else
            error("Could not locate symbol $s needed for constraint.")
        end
    end

    (syms, outs, u, ins, p, obs) -> begin
        for (iidx, outnr, outidx) in outmapping
            syms[iidx] = outs[outnr][outidx]
        end
        for (iidx, unidx) in umapping
            syms[iidx] = u[unidx]
        end
        for (iidx, inr, inidx) in inmapping
            syms[iidx] = ins[inr][inidx]
        end
        for (iidx, pnidx) in pmapping
            syms[iidx] = p[pnidx]
        end
        for (iidx, obsidx) in obsmapping
            syms[iidx] = obs[obsidx]
        end
    end
end

####
#### InitFormula
####
"""
    InitFormula(f, outsym, sym)

A representation of initialization formulas that are applied during the initialization phase of a component.
InitFormulas act earlier in the initialization pipeline than InitConstraints - they essentially set additional
defaults rather than adding equations to the nonlinear system.

It contains a function `f` that defines the formulas, a vector of output symbols `outsym` that will be set by the formulas,
a vector of input symbols `sym` that are used in the formulas, and an optional pretty-print string.

    InitFormula([:Vset], [:u_r, :u_i]) do out, u
        out[:Vset] = sqrt(u[:u_r]^2 + u[:u_i]^2)
    end

See also [`@initformula`](@ref) for a macro to create such formulas.
"""
struct InitFormula{F}
    f::F
    outsym::Vector{Symbol}   # output symbols (from LHS of assignments)
    sym::Vector{Symbol}      # input symbols (from RHS of assignments)
    prettyprint::Union{Nothing,String}
    # set by `normalize` to the untouched user-written formula this one was derived from,
    # `nothing` for user-written formulas themselves. See [`normalize`](@ref).
    derived_from::Union{Nothing,InitFormula}
    weak::Bool # don't overwrite set defaults
    label::Union{Nothing,String} # short line identifier (just verbose application)

    function InitFormula(f::F, outsym::Vector{Symbol}, sym::Vector{Symbol}, prettyprint::Union{Nothing,String}, derived_from::Union{Nothing,InitFormula}, weak::Bool=false, label::Union{Nothing,String}=nothing) where F
        # Check for self-dependencies (formula depending on its own output)
        self_deps = intersect(sym, outsym)
        if !isempty(self_deps)
            throw(ArgumentError("InitFormula cannot depend on its own output symbols: $self_deps"))
        end
        # weak defaulting is single-target: a weak formula is dropped whole when its target is
        # already backed, so a multi-output weak could strand a uniquely-pinned sibling output.
        if weak && length(outsym) != 1
            throw(ArgumentError("A weak InitFormula must have exactly one output symbol (got $outsym)."))
        end
        new{F}(f, outsym, sym, prettyprint, derived_from, weak, label)
    end
end
InitFormula(f, outsym, sym; weak::Bool=false, label=nothing) = InitFormula(f, outsym, sym, nothing, nothing, weak, label)
InitFormula(f, outsym, sym, prettyprint; weak::Bool=false, label=nothing) = InitFormula(f, outsym, sym, prettyprint, nothing, weak, label)

dim(c::InitFormula) = length(c.outsym)

(c::InitFormula)(out, u) = c.f(out, SymbolicView(u, c.sym))

function Base.show(io::IO, ::MIME"text/plain", @nospecialize(c::InitFormula))
    _show_formula(io, c)
end

"""
    GuessFormula(f, outsym, sym)

A representation of guess formulas that improve initial guesses during component initialization.
GuessFormulas are applied after InitFormulas in the pipeline, reading from both defaults and
guesses to compute improved guess values for free variables.

Unlike InitFormulas which set defaults (reducing free variables), GuessFormulas only update
the `guesses` dictionary to improve solver convergence without changing the problem dimension.

It contains a function `f` that defines the formulas, a vector of output symbols `outsym` that
will be set by the formulas, a vector of input symbols `sym` that are used in the formulas,
and an optional pretty-print string.

# Input Lookup Behavior
When evaluating a GuessFormula, input symbols are looked up with this priority:
1. First check `defaults` dict (fixed/known values take precedence)
2. Then check `guesses` dict (free variable current guesses)
3. Error if symbol not found in either

# Examples

    GuessFormula([:V, :theta], [:u_r, :u_i]) do out, u
        out[:V] = sqrt(u[:u_r]^2 + u[:u_i]^2)
        out[:theta] = atan(u[:u_i], u[:u_r])
    end

See also [`@guessformula`](@ref) for a macro to create such formulas.
"""
struct GuessFormula{F}
    f::F
    outsym::Vector{Symbol}   # output symbols (from LHS of assignments)
    sym::Vector{Symbol}      # input symbols (from RHS of assignments)
    prettyprint::Union{Nothing,String}
    # set by `normalize` to the untouched user-written formula this one was derived from,
    # `nothing` for user-written formulas themselves. See [`normalize`](@ref).
    derived_from::Union{Nothing,GuessFormula}
    label::Union{Nothing,String}  # one-line `(via …)` provenance for the verbose log, see `InitFormula.label`

    function GuessFormula(f::F, outsym::Vector{Symbol}, sym::Vector{Symbol}, prettyprint::Union{Nothing,String}, derived_from::Union{Nothing,GuessFormula}, label::Union{Nothing,String}=nothing) where F
        # Check for self-dependencies (formula depending on its own output)
        self_deps = intersect(sym, outsym)
        if !isempty(self_deps)
            throw(ArgumentError("GuessFormula cannot depend on its own output symbols: $self_deps"))
        end
        new{F}(f, outsym, sym, prettyprint, derived_from, label)
    end
end
GuessFormula(f, outsym, sym; label=nothing) = GuessFormula(f, outsym, sym, nothing, nothing, label)
GuessFormula(f, outsym, sym, prettyprint; label=nothing) = GuessFormula(f, outsym, sym, prettyprint, nothing, label)

dim(c::GuessFormula) = length(c.outsym)

(c::GuessFormula)(out, u) = c.f(out, SymbolicView(u, c.sym))

function Base.show(io::IO, ::MIME"text/plain", @nospecialize(c::GuessFormula))
    _show_formula(io, c)
end

# A normalized formula carries symbol lists the user never wrote, so printing its own
# `prettyprint` (which spells the original symbols) alone would misrepresent it. Show the
# original recipe and append the symbols actually in play. `derived_from` is never nested,
# so this recurses once.
function _show_formula(io::IO, @nospecialize(c))
    if isnothing(c.derived_from)
        _show_recipe(io, c)
    else
        _show_formula(io, c.derived_from)
        orig = c.derived_from
        print(io, "\n(normalized: $(c.sym) → $(c.outsym), \
                   derived from $(orig.sym) → $(orig.outsym))")
    end
end

"""
   @initformula

Generate an [`InitFormula`](@ref) from an expression using symbols.

    @initformula begin
        :Vset = sqrt(:u_r^2 + :u_i^2)
        :Pset = :u_r * :i_r + :u_i * :i_i
    end

is equal to

    InitFormula([:Vset, :Pset], [:u_r, :u_i, :i_r, :i_i]) do out, u
        out[:Vset] = sqrt(u[:u_r]^2 + u[:u_i]^2)
        out[:Pset] = u[:u_r] * u[:i_r] + u[:u_i] * u[:i_i]
    end
"""
macro initformula(args...)
    isempty(args) && throw(ArgumentError("@initformula expects a formula block"))
    ex = last(args)
    weak = false
    for opt in args[1:end-1]
        if opt isa Expr && opt.head == :(=) && opt.args[1] === :weak
            weak = opt.args[2]
        else
            throw(ArgumentError("@initformula: unexpected option `$opt`, only `weak=true/false` is supported"))
        end
    end
    _formula_macro(InitFormula, ex; weak)
end


"""
   @guessformula

Generate a [`GuessFormula`](@ref) from an expression using symbols.

This macro provides convenient syntax for creating guess formulas using assignment
expressions with quoted symbols. Each assignment computes a guess value for a free variable.

# Syntax

```julia
@guessformula begin
    :output1 = expression_with(:input_symbols)
    :output2 = other_expression(:more_inputs)
end
```

# Input Symbol Lookup
Input symbols (on RHS) are looked up with this priority:
1. First from `defaults` dict (fixed values)
2. Then from `guesses` dict (current guesses)
3. Error if not found in either

# Output Symbol Target
Output symbols (on LHS) must be:
- Valid component symbols (states, parameters, inputs, outputs)
- NOT observables (observables are computed, not guessed)

See also: [`GuessFormula`](@ref), [`@initformula`](@ref), [`initialize_component`](@ref)
"""
macro guessformula(ex)
    _formula_macro(GuessFormula, ex)
end

# combined macro code for InitFormula and GuessFormula
#=
In a nutshell: wrap every QuoteNote symbol either in u[:sym] or out[:sym]
- if it appears alone in the LHS of an assignment, it is an output symbol
- otherwise it is an input symbol
some thinks will break!
for example: set!(:out, :in) -> set!(u[:out], u[:in])
=#
function _formula_macro(type, ex; weak=false)
    if ex isa QuoteNode || ex.head != :block
        ex = Base.remove_linenums!(Expr(:block, ex))
    end

    macroname = type === InitFormula ? "@initformula" : "@guessformula"
    header = (type === InitFormula && weak === true) ? "$macroname weak=true" : macroname
    s = _macro_source_string(header, ex)
    lbl = _auto_formula_label(ex)

    input_syms = Symbol[]    # RHS symbols
    output_syms = Symbol[]   # LHS symbols
    out_var = :__out__
    u = :__u__
    body = Expr[]

    for formula in ex.args
        formula isa Union{Expr, QuoteNode} || continue # skip line number nodes
        if formula isa Expr && formula.head == :(=)
            rhs = formula.args[2]
            # Wrap symbols in the RHS
            wrapped_rhs = _wrap_symbols!(rhs, input_syms, u)

            lhs = formula.args[1]
            # Extract the target symbol from the LHS (should be a QuoteNode)
            if lhs isa QuoteNode && lhs.value isa Symbol
                # :output = ... assigmend
                target_sym = lhs.value
                push!(output_syms, target_sym)
                push!(body, :($out_var[$(QuoteNode(target_sym))] = $wrapped_rhs))
            else
                # "normal" assigmend
                push!(body, :($lhs = $wrapped_rhs))
            end
        else
            wrapped_expr = _wrap_symbols!(formula, input_syms, u)
            push!(body, :($wrapped_expr))  # standalone expression
        end
    end
    unique!(input_syms)

    body_esc = _escape_all.(body)

    closure = quote
        function ($(esc(out_var)), $(esc(u)))
            $(body_esc...)
            nothing
        end
    end
    # only InitFormula carries `weak`; GuessFormula's constructor has no such kwarg
    if type === InitFormula
        :($(type)($closure, $output_syms, $input_syms, $s; weak = $(esc(weak)), label = $lbl))
    else
        :($(type)($closure, $output_syms, $input_syms, $s; label = $lbl))
    end
end

# for metadata check, validates both input and output symbols
assert_initformula_compat(cf::ComponentModel, c::InitFormula) = _assert_formula_compat(cf, c)
assert_guessformula_compat(cf::ComponentModel, c::GuessFormula) = _assert_formula_compat(cf, c)

# Inputs may be observables: they are expanded to their settable roots at init time, see
# `normalize`. Outputs may be observables too — as pure aliases they canonicalize onto a
# settable symbol, and as genuinely algebraic observables the write *pins* them as an
# init-time dataflow node, see `pinned_obssyms`.
function _assert_formula_compat(cf::ComponentModel, c::Union{InitFormula,GuessFormula})
    label = c isa InitFormula ? "InitFormula" : "GuessFormula"
    settable = settable_symbols(cf)

    input_mismatch = setdiff(c.sym, settable ∪ obssym(cf))
    if !isempty(input_mismatch)
        throw(ArgumentError("$label uses input symbols not part of component model: $input_mismatch"))
    end

    output_mismatch = setdiff(c.outsym, settable ∪ keys(get_aliasmap(cf)) ∪ obssym(cf))
    if !isempty(output_mismatch)
        throw(ArgumentError("$label output symbols must be settable component symbols, \
                             aliases thereof or observables (pins), but these are not: \
                             $output_mismatch"))
    end

    c
end

"""
    topological_sort_formulas(formulas::Vector{Init/GuessFormula}) -> Vector{Init/GuessFormula}

Sort formulas in topological order based on their dependencies. This ensures that
formulas are applied in the correct order, where each formula's input symbols are
available before it executes.

A formula B depends on formula A if any of B's input symbols are produced by A's output symbols.
The function returns the formulas reordered such that dependencies are satisfied.
"""
function topological_sort_formulas(formulas)
    n = length(formulas)
    n <= 1 && return copy(formulas)

    type = if first(formulas) isa InitFormula
        "InitFormula"
    elseif first(formulas) isa GuessFormula
        "GuessFormula"
    else
        throw(ArgumentError("Expected a vector of InitFormula or GuessFormula, got $(typeof(first(formulas)))"))
    end

    # Check for output symbol conflicts
    all_outputs = Symbol[]
    for formula in formulas
        append!(all_outputs, formula.outsym)
    end

    if !allunique(all_outputs)
        conflicts = [s for s in unique(all_outputs) if count(==(s), all_outputs) > 1]
        # A weak writer yields to a `default` or a *strong* co-writer (both drop it earlier), so
        # reaching here means the colliding writers are all weak — a genuine ambiguity.
        hint = if any(f -> f isa InitFormula && f.weak && !isdisjoint(f.outsym, conflicts), formulas)
            "\nThe colliding formulas are `weak`: a weak formula yields to a default or a strong \
             writer, but here only weak writers target the symbol. Give the target a default, or \
             drop one of the writers."
        else
            ""
        end
        throw(ArgumentError("Multiple $(type)s set the same symbol(s): $conflicts$hint"))
    end

    # Build dependency graph using Graphs.jl
    g = SimpleDiGraph(n)

    # Add edges: i → j means formula j depends on formula i
    # (formula j must run after formula i)
    for i in 1:n, j in 1:n
        if i != j
            # Check if formula j depends on formula i
            # i.e., j's input symbols intersect with i's output symbols
            if !isdisjoint(formulas[j].sym, formulas[i].outsym)
                add_edge!(g, i, j)
            end
        end
    end

    # Perform topological sort
    try
        sorted_indices = Graphs.topological_sort(g)
        if length(sorted_indices) != n
            # This shouldn't happen with a proper topological sort implementation,
            # but let's be safe
            throw(ArgumentError(_circular_dependency_msg(formulas, g, type)))
        end
        return formulas[sorted_indices]
    catch e
        if e isa ErrorException && occursin("graph contains at least one loop", string(e))
            throw(ArgumentError(_circular_dependency_msg(formulas, g, type)))
        else
            rethrow(e)
        end
    end
end

function _circular_dependency_msg(formulas, g, type)
    # cycle count is bounded by the (small) number of formulas on one component; take the
    # shortest witness, one is enough to explain the failure
    cycles = Graphs.simplecycles_limited_length(g, Graphs.nv(g))
    isempty(cycles) && return "Circular dependency detected in $(type)s"
    cycle = argmin(length, cycles)

    rows = map(eachindex(cycle)) do i
        from = formulas[cycle[i]]
        to = formulas[cycle[mod1(i + 1, length(cycle))]]
        shared = intersect(to.sym, from.outsym)
        "  $(_cycle_ref(from)) writes $shared, read by $(_cycle_ref(to))"
    end
    normalized = filter(i -> !isnothing(formulas[i].derived_from), cycle)
    note = if isempty(normalized)
        ""
    else
        "\nSome of them were normalized, i.e. their inputs are the settable roots of the \
         observables they were written against: \
         $(join(("$(_cycle_ref(formulas[i])) was written to read $(formulas[i].derived_from.sym)" for i in normalized), ", ")). \
         An observable which is written by a formula of this initialization becomes a pin and \
         is read directly instead of being expanded to its roots — so an edge above may \
         disappear once the observable it runs through is written by a formula."
    end
    "Circular dependency detected between $(length(cycle)) $(type)s:\n" * join(rows, "\n") * note
end
# `_formula_ref`'s prose form ("for [:x]") does not survive mid-sentence next to a `writes`,
# so name a formula by its quoted recipe where it has one
_cycle_ref(f) = isnothing(f.label) ? "the formula for $(f.outsym)" : "`$(f.label)`"

"""
    extend_knowns_by_formulas!(knowns, cf, formulas; am, t, targets, error_unresolvable, verbose, io)

Extend `knowns` in place by everything the resolution graph derives from it: the component's
`:obsrules` together with the user's init `formulas`, resolved by [`resolve_rules`](@ref).

`targets` narrows what is worth computing — an `NWState` holds states and parameters only.
Which symbols have a *slot* is the component's full settable set either way.

`t=NaN` means "no time given" and is translated here: `:t` is then not seeded at all, so a rule
reading it stays unfired instead of quietly producing `NaN`.

See [`_resolution_rules`](@ref) for when the pass runs at all.
"""
function extend_knowns_by_formulas!(knowns, cf, formulas;
                                    am=get_aliasmap(cf), t=NaN, targets=nothing,
                                    error_unresolvable=true, verbose=false, io=stdout)
    settable = settable_symbols(cf)
    if isnothing(targets)
        targets = settable
    end
    # collect rules (user formulas + obs rules)
    rules = _resolution_rules(cf, formulas, am, knowns, settable)
    isempty(rules) && return knowns

    seed = Dict{Symbol,Float64}(knowns)
    isnan(t) || (seed[:t] = t) # `t` so that formulas and observables can read it
    # final targets: prune away any symbol already known and known to be not overwritten
    targets = _prune_resolution_targets(targets, rules, seed)
    # actually resolve rules (result stored in res)
    res = resolve_rules(seed, rules; targets)
    _check_resolution(res, rules; error_unresolvable, verbose, io)

    verbose && _print_resolution(res, rules, knowns, settable; io)
    for (s, v) in res.vals
        p = res.provenance[s]
        p === :provided && continue # an untouched seed, already in `knowns`
        # a `:derived` observable is a scratch value of the walk; one written by a formula is a
        # pin, which does belong in `knowns`
        (s ∈ settable || p !== :derived) || continue
        knowns[s] = v
    end
    return knowns
end

"""
    extend_guesses_by_formulas!(guesses, defaults, cf, formulas; am, t, targets, error_unresolvable, verbose, io)

Guess-formula sibling of [`extend_knowns_by_formulas!`](@ref): extend `guesses` in place by
everything the resolution graph derives from the component's `:obsrules` together with the
user's guess `formulas`.

Same executor, different roots. The **full output of the init pass** (`defaults`) is seeded as
`:provided` — everything it determined outranks any guess and is never moved — while a value
already sitting in `guesses` is seeded at `:guess`, the very rank a `GuessFormula` writes at.
With `overwrite_equal=true` a guess formula therefore *refines* an existing guess instead of
colliding with it: a guess is a seed, not an assertion. Which of two guess formulas writing one
symbol wins is deliberately unspecified.

An observed equation stays `:derived` and hence *below* the guesses: it may fill a symbol
nothing else reaches, but it never replaces a guess that was handed to us.
"""
function extend_guesses_by_formulas!(guesses, defaults, cf, formulas;
                                     am=get_aliasmap(cf), t=NaN, targets=nothing,
                                     error_unresolvable=true, verbose=false, io=stdout)
    settable = settable_symbols(cf)
    if isnothing(targets)
        targets = settable
    end
    # layered seed, `defaults` over `guesses` — the same lookup the guess pass always had
    seed = Dict{Symbol,Float64}(guesses)
    merge!(seed, defaults)
    # collect rules (user formulas + obs rules)
    rules = _resolution_rules(cf, formulas, am, seed, settable)
    isempty(rules) && return guesses

    # the init output is bound, a pre-existing guess is not; `:guess` is what a `GuessFormula`
    # rule writes at, so with `overwrite_equal` a formula may refine it
    bound = Set{Symbol}(keys(defaults))
    if !isnan(t)
        seed[:t] = t # `t` so that formulas and observables can read it
        push!(bound, :t)
    end
    seedprov = s -> s ∈ bound ? :provided : :guess
    # final targets: prune away any symbol already known and known to be not overwritten
    targets = _prune_resolution_targets(targets, rules, seed; seedprov, overwrite_equal=true)
    # actually resolve rules (result stored in res)
    res = resolve_rules(seed, rules; targets, seedprov, overwrite_equal=true)
    _check_resolution(res, rules; error_unresolvable, verbose, io, type="GuessFormula")

    verbose && _print_resolution(res, rules, guesses, settable;
                                 io, type="GuessFormula", op="≈", fixed=defaults)
    for (s, v) in res.vals
        p = res.provenance[s]
        p === :provided && continue # bound by the init pass (or the seeded `t`), not a guess
        (s ∈ settable || p !== :derived) || continue
        guesses[s] = v
    end
    return guesses
end

# Must run *before* the writeback, so `knowns` still holds the value each row reports as "(was
# …)". Two groups: what the user's formulas asserted, and what the model's own equations
# determined on the way — the latter is new with the graph.
function _print_resolution(res, rules, knowns, settable; io, type="InitFormula", op="=", fixed=nothing)
    formula_rows = String[]
    for (i, r) in enumerate(rules)
        _is_user_rule(r) || continue
        for s in r.outsym
            # `res.writer`, not the provenance: two same-tier formulas share a provenance, and a
            # formula that never fired keeps none of its outputs — neither may claim the write
            get(res.writer, s, 0) == i || continue
            push!(formula_rows, _formula_row(s, res.vals[s], knowns;
                                             op, fixed, pin=s ∉ settable, label=r.label))
        end
    end
    print_aligned_group(io, "$(type)s set:", formula_rows)

    derived_rows = [_formula_row(s, v, knowns; op, fixed)
                    for (s, v) in sort!(collect(res.vals); by=first)
                    if res.provenance[s] === :derived && s ∈ settable]
    print_aligned_group(io, "Derived from the component's own equations:", derived_rows)
    nothing
end

"""
    _resolution_rules(cf, formulas, am, knowns, settable) -> Vector{ResolutionRule}

The rule bucket to apply. Init stays cheap when there is nothing to resolve, so the bucket is
non-empty only if

- the user provided a formula, or
- some value sits on a non-settable symbol — a `default` on the observable `y` of `y ~ -x`
  reaches the state `x` only through the inverse rule.

Obs rules go first so a user formula wins the arbitration between two rules ready in the same
pass.
"""
function _resolution_rules(cf, formulas, am, knowns, settable)
    hasformulas = !isnothing(formulas) && !isempty(formulas)
    if !hasformulas && all(s -> s ∈ settable, keys(knowns))
        return ResolutionRule[]
    end

    rules = collect(ResolutionRule, get_obsrules(cf))
    if hasformulas
        for f in formulas
            push!(rules, ResolutionRule(f, am))
        end
    end
    rules
end

"""
    _prune_resolution_targets(targets, rules, vals; seedprov, overwrite_equal) -> Set{Symbol}

Targets are what the rules are pruned for. Narrows `targets` down to:
- the ones written by some rule
- MINUS the symbols which are already in `vals` and outrank the rules writing them

(for example, a user provided default can only be changed by a strong formula, so
we can remove it as target UNLESS there is such a strong formula in the ruleset)

`seedprov` is what the seed's entries are worth, mirroring `resolve_rules`' `seed_provenance`:
the init pass provides everything, the guess pass only the init output (the rest is a guess,
which a guess rule of equal rank may still refine). The comparison is therefore the same one
`_write_value!` makes — equal rank changes nothing unless `overwrite_equal`.
"""
function _prune_resolution_targets(targets, rules, vals;
                                   seedprov=Returns(:provided), overwrite_equal=false)
    best = Dict{Symbol,Int}()
    for r in rules, s in r.outsym
        best[s] = max(get(best, s, 0), _precedence(r.provenance))
    end

    kept = Set{Symbol}()
    for (s, rank) in best
        s ∈ targets || continue
        # `haskey`, matching `resolve_rules`: a seeded NaN is a value like any other
        if haskey(vals, s)
            held = _precedence(seedprov(s))
            (held > rank || (held == rank && !overwrite_equal)) && continue
        end
        push!(kept, s)
    end
    kept
end

# Caller-side policy, the half `resolve_rules` deliberately does not have: a rule that left one
# of its outputs unknown, two rules of equal precedence disagreeing about one, a required rule
# that turned out not to be needed, and values that came out non-finite.
function _check_resolution(res, rules; error_unresolvable, verbose, io, type="InitFormula")
    _check_pruned(res, rules; error_unresolvable, verbose, io, type)
    _check_nonfinite(res, rules; error_unresolvable, verbose, io)

    if !isempty(res.conflicts)
        rows = ["$(_rule_ref(rules[c.rule])) computed :$(c.sym) = $(c.offered), \
                 but it already held $(c.held)" for c in res.conflicts]
        msg = "Inconsistent initialization values:\n" * join("  - " .* rows, "\n")
        if error_unresolvable
            throw(ArgumentError(msg))
        else
            verbose && printstyled(io, " - $msg\n"; color=:yellow)
        end
    end

    for u in res.unfired
        r = rules[u.rule]
        # not firing is a failure only if something is still unknown — with duplicate writers
        # allowed, the other one may simply have got there first
        all(s -> haskey(res.vals, s), r.outsym) && continue
        msg = "$type $(_rule_ref(r)) could not be resolved: its inputs are unknown: \
               $(_resolution_hint(res, rules, u.unknown))."
        if error_unresolvable
            throw(ArgumentError(msg))
        else
            verbose && printstyled(io, " - skipping: $msg\n")
        end
    end
    nothing
end
# warn for on nonfininte values
function _check_nonfinite(res, rules; error_unresolvable, verbose, io)
    isempty(res.nonfinite) && return nothing
    (error_unresolvable || verbose) || return nothing
    rows = map(res.nonfinite) do nf
        src = iszero(nf.rule) ? "provided" : "computed by $(_rule_ref(rules[nf.rule]))"
        ":$(nf.sym) &= $(res.vals[nf.sym]) &($src)"
    end
    printstyled(io, " - WARNING: ", color=:yellow)
    printstyled(io, "non-finite values, which poison everything derived from them. A common \
                     cause is an observable with an explicit time dependence and no `t` \
                     passed to the initialization.\n")
    print_aligned_rows(io, rows)
    nothing
end
# warn if non-optional rule pruned away, weak formulas are optional
function _check_pruned(res, rules; error_unresolvable, verbose, io, type)
    for (i, r) in enumerate(rules)
        (res.pruned[i] && _is_user_rule(r)) || continue
        # both kinds that legitimately yield: a weak init formula, and any guess rule whose
        # target the init pass already determined
        if r.provenance === :weak_formula || r.provenance === :guess
            verbose && printstyled(io, " - $type: $(_rule_ref(r)) yields, \
                                        its target $(r.outsym) is already determined\n")
        elseif error_unresolvable # the init path; `NWState` reconstruction narrows the slot
            # set to states and parameters, where a formula writing an output prunes routinely
            printstyled(io, " - WARNING: ", color=:yellow)
            printstyled(io, "$(_rule_ref(r)) had no effect: nothing in this component reads \
                             $(r.outsym), so it cannot change any state, parameter or \
                             output.\n")
        else
            verbose && printstyled(io, " - skipping: $(_rule_ref(r)), nothing reads \
                                        $(r.outsym)\n")
        end
    end
    nothing
end
# One level backwards from each missing input: either nothing computes it, or something does
# but is itself starved — the "to compute :δ, provide one of {…}" material.
function _resolution_hint(res, rules, unknown)
    parts = map(unknown) do s
        src = [r for r in rules if s ∈ r.outsym]
        isempty(src) && return ":$s (nothing computes it, provide it directly)"
        need = unique(Symbol[x for r in src for x in r.sym if !haskey(res.vals, x)])
        isempty(need) && return ":$s"
        # `t` is the one missing input nobody can provide a default for — it is a call argument
        :t ∈ need && return ":$s (equations that depend explicitly on time need a concrete \
                               `t`, which this initialization was not given)"
        ":$s (computable from $need)"
    end
    join(parts, ", ")
end

# A weak formula yields — and is dropped — when its (canonical) target already carries a
# `default` or is written by a *strong* formula (`strong_outputs`): an InitFormula always fires,
# so a strong writer pins the target and the weak default is redundant. The default check is on
# `default` only, never `init` — a weak formula persists its own output as an `init`, and testing
# `init` here would self-block it on reinit. Weak formulas are single-output by construction, so
# a drop never strands a uniquely-pinned sibling output.
function drop_weak_formulas(formulas, defaults, strong_outputs=(); verbose=false, io=stdout)
    any(f -> f isa InitFormula && f.weak, formulas) || return formulas
    kept = empty(formulas)
    for f in formulas
        if !f.weak
            push!(kept, f)
            continue
        end
        @assert length(f.outsym) == 1 "weak InitFormula must be single-output (enforced at construction)"
        s = only(f.outsym)
        if haskey(defaults, s) || s in strong_outputs
            ref = isnothing(f.label) ? "weak formula for :$s" : f.label
            verbose && printstyled(io, " - InitFormula: $ref yields to \
                $(haskey(defaults, s) ? "existing default :$s = $(defaults[s])" : "a strong formula writing :$s")\n")
        else
            push!(kept, f)
        end
    end
    kept
end

function apply_init_formulas!(defaults, formulas_unsorted; verbose=false, io=stdout,
                              error_unresolvable=true, pinned=Set{Symbol}())
    # Convert tuple to vector if necessary
    formulas_vec = formulas_unsorted isa Tuple ? collect(formulas_unsorted) : formulas_unsorted
    formulas = topological_sort_formulas(formulas_vec)

    rows = String[]
    for f in formulas
        out = SymbolicView(zeros(length(f.outsym)), f.outsym)
        # ensure all input symbols are in defaults
        invals = map(f.sym) do s
            haskey(defaults, s) ? defaults[s] : NaN
        end
        if any(v -> ismissing(v) || isnothing(v) || isnan(v), invals)
            if error_unresolvable
                throw(ArgumentError("InitFormula requires all input symbols to be initialized, but found NaN/missing/nothing in inputs: $(f.sym .=> invals)" * _unresolved_note(f)))
            else
                verbose && printstyled(io, " - InitFormula: skipping formula $(_formula_ref(f)) with unresolvable inputs: $(f.sym .=> invals)$(_unresolved_note(f))\n")
                continue
            end
        end
        in = SymbolicView(invals, f.sym)
        if !_run_formula!(out, f, in, "InitFormula"; verbose, io, error_unresolvable)
            continue
        end
        for s in f.outsym
            val = out[s]
            verbose && push!(rows, _formula_row(s, val, defaults; op="=", pin=s ∈ pinned, label=f.label))
            defaults[s] = val
        end
    end
    verbose && print_aligned_group(io, "InitFormulas set:", rows)
    return defaults
end
function apply_guess_formulas!(guesses, defaults, formulas_unsorted; verbose=false, io=stdout,
                               error_unresolvable=true, pinned=Set{Symbol}())
    # Convert tuple to vector if necessary
    formulas_vec = formulas_unsorted isa Tuple ? collect(formulas_unsorted) : formulas_unsorted
    formulas = topological_sort_formulas(formulas_vec)

    rows = String[]
    for f in formulas
        out = SymbolicView(zeros(length(f.outsym)), f.outsym)
        # Layered lookup: defaults (fixed) take precedence over guesses
        invals = map(f.sym) do s
            if haskey(defaults, s)
                defaults[s]  # Use fixed default value if available
            elseif haskey(guesses, s)
                guesses[s]   # Otherwise use guess value
            else
                NaN
            end
        end
        # Validate inputs are not NaN/missing/nothing
        if any(v -> ismissing(v) || isnothing(v) || isnan(v), invals)
            if error_unresolvable
                throw(ArgumentError("GuessFormula requires all input symbols to be initialized, but found NaN/missing/nothing in inputs: $(f.sym .=> invals)" * _unresolved_note(f)))
            else
                verbose && printstyled(io, " - GuessFormula: skipping formula $(_formula_ref(f)) with unresolvable inputs: $(f.sym .=> invals)$(_unresolved_note(f))\n")
                continue
            end
        end
        in = SymbolicView(invals, f.sym)
        if !_run_formula!(out, f, in, "GuessFormula"; verbose, io, error_unresolvable)
            continue
        end
        # Update guesses dictionary (NOT defaults!)
        for s in f.outsym
            val = out[s]
            verbose && push!(rows, _formula_row(s, val, guesses; op="≈", fixed=defaults, pin=s ∈ pinned, label=f.label))
            guesses[s] = val
        end
    end
    verbose && print_aligned_group(io, "GuessFormulas set:", rows)
    return guesses
end

# One aligned row for a formula that wrote `val` onto `:s`, annotated with what it did
# relative to what was there: nothing for a fresh write, the previous value if it changed
# one, or a no-effect note when a fixed default (guesses only) shadows the write. `label` (the
# formula's short name) is appended as its own trailing `(via …)` column so provenance lines up
# across rows; rows whose formula has no label leave the column empty (rstripped away on print).
function _formula_row(s, val, prev; op, fixed=nothing, pin=false, label=nothing)
    v = str_significant(val; sigdigits=5, phantom_minus=true)
    note = if !isnothing(fixed) && haskey(fixed, s)
        "(no effect, fixed at $(str_significant(fixed[s]; sigdigits=5)))"
    elseif haskey(prev, s)
        prev[s] ≈ val ? "(unchanged)" : "(was $(str_significant(prev[s]; sigdigits=5)))"
    else
        ""
    end
    pin && (note = strip(note * " (pinned observable)"))
    via = isnothing(label) ? "" : "(via $label)"
    ":$s &$op $v &$note &$via"
end

# Runs `f`, returning whether its outputs may be used. A normalized formula only discovers
# during the gather that an observable input expands to NaN, which is an unresolvable input
# like any other — so it lands here rather than escaping, and is treated exactly as the
# up-front dict lookup treats a missing value.
function _run_formula!(out, f, u, label; verbose, io, error_unresolvable)
    try
        f(out, u)
        true
    catch e
        e isa UnresolvableExpansionError || rethrow()
        if error_unresolvable
            throw(ArgumentError("$label for $(f.outsym): $(e.msg)" * _unresolved_note(f)))
        end
        verbose && printstyled(io, " - $label: skipping formula $(_formula_ref(f)), $(e.msg)$(_unresolved_note(f))\n")
        false
    end
end

# Short human reference to a formula for the verbose log's prose lines (skips, weak-drops): its
# `label` when it has one, else a `for [:out…]` fallback naming the output symbols.
_formula_ref(f) = isnothing(f.label) ? "for $(f.outsym)" : f.label

# A normalized formula reads the settable *roots* of what the user asked for, so the symbols
# named in an unresolvable-input message are not the ones they wrote. Say where they came
# from, and name the fix: a default on the observable itself would not help, since formulas
# read observables from the model rather than from the dict.
function _unresolved_note(f)
    orig = f.derived_from
    isnothing(orig) && return ""
    "\nNote: this formula was normalized, its inputs $(f.sym) are the settable roots of the \
     originally requested $(orig.sym). Defaults on observables are not consumed as formula \
     inputs — provide defaults for the roots instead."
end

_vcattable(t::Tuple) = collect(t)
_vcattable(x) = x # Formula types and AbstractVector

# Generic formula collection function (works for InitFormula, GuessFormula, etc.)
collect_formulas(::Nothing, ::Nothing) = nothing
collect_formulas(::Nothing, f) = vcat(_vcattable(f))
collect_formulas(f, ::Nothing) = vcat(_vcattable(f))
collect_formulas(a, b) = vcat(_vcattable(a), _vcattable(b))

merge_initconstraints(::Nothing, ::Nothing) = nothing
merge_initconstraints(::Nothing, c::InitConstraint) = c
merge_initconstraints(::Nothing, cs) = InitConstraint(cs...)

merge_initconstraints(c ::InitConstraint, ::Nothing) = c
merge_initconstraints(cA::InitConstraint, cB::InitConstraint) = InitConstraint(cA, cB)
merge_initconstraints(c ::InitConstraint, cs) = InitConstraint(c, cs...)

merge_initconstraints(cs, ::Nothing) = InitConstraint(cs...)
merge_initconstraints(cs, c::InitConstraint) = InitConstraint(cs..., c)
merge_initconstraints(csA, csB) = InitConstraint(csA..., csB...)
