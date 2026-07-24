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

    function InitFormula(f::F, outsym::Vector{Symbol}, sym::Vector{Symbol}, prettyprint::Union{Nothing,String}, derived_from::Union{Nothing,InitFormula}, weak::Bool=false) where F
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
        new{F}(f, outsym, sym, prettyprint, derived_from, weak)
    end
end
InitFormula(f, outsym, sym; weak::Bool=false) = InitFormula(f, outsym, sym, nothing, nothing, weak)
InitFormula(f, outsym, sym, prettyprint; weak::Bool=false) = InitFormula(f, outsym, sym, prettyprint, nothing, weak)

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

    function GuessFormula(f::F, outsym::Vector{Symbol}, sym::Vector{Symbol}, prettyprint::Union{Nothing,String}, derived_from::Union{Nothing,GuessFormula}) where F
        # Check for self-dependencies (formula depending on its own output)
        self_deps = intersect(sym, outsym)
        if !isempty(self_deps)
            throw(ArgumentError("GuessFormula cannot depend on its own output symbols: $self_deps"))
        end
        new{F}(f, outsym, sym, prettyprint, derived_from)
    end
end
GuessFormula(f, outsym, sym) = GuessFormula(f, outsym, sym, nothing, nothing)
GuessFormula(f, outsym, sym, prettyprint) = GuessFormula(f, outsym, sym, prettyprint, nothing)

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
        :($(type)($closure, $output_syms, $input_syms, $s; weak = $(esc(weak))))
    else
        :($(type)($closure, $output_syms, $input_syms, $s))
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
            throw(ArgumentError("Circular dependency detected in $(type)s"))
        end
        return formulas[sorted_indices]
    catch e
        if e isa ErrorException && occursin("graph contains at least one loop", string(e))
            throw(ArgumentError("Circular dependency detected in $(type)s"))
        else
            rethrow(e)
        end
    end
end

"""
    extend_knowns_by_formulas!(knowns, cf, formulas; am, t, pinned, error_unresolvable, verbose, io)

Extend a dict of known component values `knowns` in place by applying init `formulas`, mutating and
returning `knowns`. This is the shared formula pass behind both `initialize_component` and the
`NWState` reconstruction (`_get_appropriate_dict`): compute the pin-set from the *complete* formula
set, [`normalize`](@ref) each formula onto canonical settable symbols, [`drop_weak_formulas`](@ref)
whose target is already backed (a value in `knowns`, or a strong co-writer), then
[`apply_init_formulas!`](@ref) the survivors. A no-op when `formulas` is `nothing`.

The one seam between the callers is *what they pass as `knowns`*: `initialize_component` passes the
default-only dict (so weak formulas re-fire on reinit and the solve refines the rest); `NWState`
passes defaults-and-inits (so a reconstruction reproduces the post-init state). Keeping the
normalize / weak-drop / pin-set machinery here — not copied into each caller — is what stops the two
paths from silently drifting apart (as they once did on the weak-drop, see git history).
"""
function extend_knowns_by_formulas!(knowns, cf, formulas;
                                    am=get_aliasmap(cf), t=NaN,
                                    pinned=pinned_obssyms(formulas, cf),
                                    error_unresolvable=true, verbose=false, io=stdout)
    isnothing(formulas) && return knowns
    normed = [normalize(f, am, cf; t, pinned) for f in formulas]
    # a weak formula also yields to a *strong* co-writer on the same target (see
    # `drop_weak_formulas`); computed post-normalize so the outputs are canonical
    strong_out = Set(s for f in normed if !f.weak for s in f.outsym)
    kept = drop_weak_formulas(normed, knowns, strong_out; verbose, io)
    isempty(kept) || apply_init_formulas!(knowns, kept; error_unresolvable, verbose, io, pinned)
    return knowns
end

"""
    extend_guesses_by_formulas!(guesses, defaults, cf, formulas; am, t, init_pinned, error_unresolvable, verbose, io)

Guess-formula sibling of [`extend_knowns_by_formulas!`](@ref), shared by `initialize_component` and
the `NWState` reconstruction (`_get_appropriate_dict`). Extends the `guesses` dict in place by
`normalize`-ing and applying the guess `formulas`, mutating and returning `guesses`. A no-op when
`formulas` is `nothing`.

Differs from the init pass in two structural ways, which is why it is a separate function rather
than the same one: guess application is a *layered* two-dict write — [`apply_guess_formulas!`](@ref)
writes `guesses` but reads `defaults`-before-`guesses` as inputs (fixed values win) — and there is
no weak-drop (`GuessFormula` has no `weak`).

The guess frontier is `(init_pinned ∩ keys(defaults)) ∪ guess_pins`. The caller passes the full
static *init* pin-set in (`initialize_component` reuses the one it already computed — which covers
its `additional_initformula`s too; `NWState` passes `pinned_obssyms(cm)`), and it is **intersected
with what actually landed in `defaults`**: an init formula that did not run (skipped on unresolvable
inputs under `error_unresolvable=false`, i.e. the NWState path) leaves its pinned observable
*unwritten*, so a guess reading it must expand to roots rather than stop at a value that is not
there. This filtering is only valid for the init pins — the init pass has already run and
materialized its outputs, so `keys(defaults)` is ground truth for "did this pin happen." The guess
pins cannot be filtered the same way: guess formulas are all normalized before any is applied, so
their pin-set must stay the full static set (order-independent *within* the guess pass, exactly as
init pins are within the init pass). Kept here so the pin-set/normalize rule is single-sourced and
the two callers cannot drift.
"""
function extend_guesses_by_formulas!(guesses, defaults, cf, formulas;
                                     am=get_aliasmap(cf), t=NaN, init_pinned=Set{Symbol}(),
                                     error_unresolvable=true, verbose=false, io=stdout)
    isnothing(formulas) && return guesses
    guess_pinned = (init_pinned ∩ keys(defaults)) ∪ pinned_obssyms(formulas, cf)
    normed = [normalize(f, am, cf; t, pinned=guess_pinned) for f in formulas]
    apply_guess_formulas!(guesses, defaults, normed; error_unresolvable, verbose, io, pinned=guess_pinned)
    return guesses
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
            verbose && printstyled(io, " - InitFormula: weak formula for :$s yields to \
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
                verbose && printstyled(io, " - InitFormula: skipping formula for $(f.outsym) with unresolvable inputs: $(f.sym .=> invals)$(_unresolved_note(f))\n")
                continue
            end
        end
        in = SymbolicView(invals, f.sym)
        if !_run_formula!(out, f, in, "InitFormula"; verbose, io, error_unresolvable)
            continue
        end
        for s in f.outsym
            val = out[s]
            verbose && push!(rows, _formula_row(s, val, defaults; op="=", pin=s ∈ pinned))
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
                verbose && printstyled(io, " - GuessFormula: skipping formula for $(f.outsym) with unresolvable inputs: $(f.sym .=> invals)$(_unresolved_note(f))\n")
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
            verbose && push!(rows, _formula_row(s, val, guesses; op="≈", fixed=defaults, pin=s ∈ pinned))
            guesses[s] = val
        end
    end
    verbose && print_aligned_group(io, "GuessFormulas set:", rows)
    return guesses
end

# One aligned row for a formula that wrote `val` onto `:s`, annotated with what it did
# relative to what was there: nothing for a fresh write, the previous value if it changed
# one, or a no-effect note when a fixed default (guesses only) shadows the write.
function _formula_row(s, val, prev; op, fixed=nothing, pin=false)
    v = str_significant(val; sigdigits=5, phantom_minus=true)
    note = if !isnothing(fixed) && haskey(fixed, s)
        "(no effect, fixed at $(str_significant(fixed[s]; sigdigits=5)))"
    elseif haskey(prev, s)
        prev[s] ≈ val ? "(unchanged)" : "(was $(str_significant(prev[s]; sigdigits=5)))"
    else
        ""
    end
    pin && (note = strip(note * " (pinned observable)"))
    ":$s &$op $v &$note"
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
        verbose && printstyled(io, " - $label: skipping formula for $(f.outsym), $(e.msg)$(_unresolved_note(f))\n")
        false
    end
end

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
