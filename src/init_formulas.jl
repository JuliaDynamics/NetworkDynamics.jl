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

# `prettyprint` is the complete authored recipe — a weak or optional `InitFormula`'s header
# already carries the flag (baked in at construction, see `_formula_macro`/`_build_formula`), so
# show just prints it. Without one (raw constructor) fall back to `repr`.
function _show_recipe(io::IO, @nospecialize(c))
    if isnothing(c.prettyprint)
        # strip trailing default-valued fields so the result stays a valid constructor call:
        # `nothing` (prettyprint) and the `false` flags. A flag that is set stays, keeping the
        # preceding `nothing` so the positional call still lines up.
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

Two independent flags say what may go wrong with it. `weak` yields on the target side: the
value is already there, so step aside. `optional` yields on the source side: the inputs never
became known, so skip the formula instead of failing the initialization.

See also [`@initformula`](@ref) for a macro to create such formulas.
"""
struct InitFormula{F}
    f::F
    outsym::Vector{Symbol}   # output symbols (from LHS of assignments)
    sym::Vector{Symbol}      # input symbols (from RHS of assignments)
    prettyprint::Union{Nothing,String}
    weak::Bool # don't overwrite set defaults
    optional::Bool # may stay unresolved
    label::Union{Nothing,String} # short line identifier (just verbose application)

    function InitFormula(f::F, outsym::Vector{Symbol}, sym::Vector{Symbol}, prettyprint::Union{Nothing,String}, weak::Bool=false, optional::Bool=false, label::Union{Nothing,String}=nothing) where F
        # Check for self-dependencies (formula depending on its own output)
        self_deps = intersect(sym, outsym)
        if !isempty(self_deps)
            throw(ArgumentError("InitFormula cannot depend on its own output symbols: $self_deps"))
        end
        new{F}(f, outsym, sym, prettyprint, weak, optional, label)
    end
end
InitFormula(f, outsym, sym; weak::Bool=false, optional::Bool=false, label=nothing) = InitFormula(f, outsym, sym, nothing, weak, optional, label)
InitFormula(f, outsym, sym, prettyprint; weak::Bool=false, optional::Bool=false, label=nothing) = InitFormula(f, outsym, sym, prettyprint, weak, optional, label)

dim(c::InitFormula) = length(c.outsym)

(c::InitFormula)(out, u) = c.f(out, SymbolicView(u, c.sym))

function Base.show(io::IO, ::MIME"text/plain", @nospecialize(c::InitFormula))
    _show_recipe(io, c)
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
    label::Union{Nothing,String}  # one-line `(via …)` provenance for the verbose log, see `InitFormula.label`

    function GuessFormula(f::F, outsym::Vector{Symbol}, sym::Vector{Symbol}, prettyprint::Union{Nothing,String}, label::Union{Nothing,String}=nothing) where F
        # Check for self-dependencies (formula depending on its own output)
        self_deps = intersect(sym, outsym)
        if !isempty(self_deps)
            throw(ArgumentError("GuessFormula cannot depend on its own output symbols: $self_deps"))
        end
        new{F}(f, outsym, sym, prettyprint, label)
    end
end
GuessFormula(f, outsym, sym; label=nothing) = GuessFormula(f, outsym, sym, nothing, label)
GuessFormula(f, outsym, sym, prettyprint; label=nothing) = GuessFormula(f, outsym, sym, prettyprint, label)

dim(c::GuessFormula) = length(c.outsym)

(c::GuessFormula)(out, u) = c.f(out, SymbolicView(u, c.sym))

function Base.show(io::IO, ::MIME"text/plain", @nospecialize(c::GuessFormula))
    _show_recipe(io, c)
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

Takes the [`InitFormula`](@ref) flags as leading options:

    @initformula weak=true :Sn = :S_b          # yield if :Sn already has a value
    @initformula optional=true :y = :K * :u    # skip if :u never becomes known
"""
macro initformula(args...)
    isempty(args) && throw(ArgumentError("@initformula expects a formula block"))
    ex = last(args)
    weak = false
    optional = false
    for opt in args[1:end-1]
        if opt isa Expr && opt.head == :(=) && opt.args[1] === :weak
            weak = opt.args[2]
        elseif opt isa Expr && opt.head == :(=) && opt.args[1] === :optional
            optional = opt.args[2]
        else
            throw(ArgumentError("@initformula: unexpected option `$opt`, only `weak=true/false` \
                                 and `optional=true/false` are supported"))
        end
    end
    _formula_macro(InitFormula, ex; weak, optional)
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
function _formula_macro(type, ex; weak=false, optional=false)
    if ex isa QuoteNode || ex.head != :block
        ex = Base.remove_linenums!(Expr(:block, ex))
    end

    macroname = type === InitFormula ? "@initformula" : "@guessformula"
    # the header is the recipe `_show_recipe` prints back, so it repeats the options as given
    opts = [o for (o, on) in (("weak=true", weak), ("optional=true", optional))
            if type === InitFormula && on === true]
    header = isempty(opts) ? macroname : "$macroname $(join(opts, " "))"
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
    # only InitFormula carries the flags; a GuessFormula is weak and optional by nature
    if type === InitFormula
        :($(type)($closure, $output_syms, $input_syms, $s;
                  weak = $(esc(weak)), optional = $(esc(optional)), label = $lbl))
    else
        :($(type)($closure, $output_syms, $input_syms, $s; label = $lbl))
    end
end

# for metadata check, validates both input and output symbols
assert_initformula_compat(cf::ComponentModel, c::InitFormula) = _assert_formula_compat(cf, c)
assert_guessformula_compat(cf::ComponentModel, c::GuessFormula) = _assert_formula_compat(cf, c)

# Observables are allowed on both ends. As an input the graph computes them, as an output the
# formula either lands on an aliased symbol or asserts something that is checked after the solve.
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
    extend_knowns_by_formulas!(knowns, cf, formulas; am, t, targets, error_unresolvable, verbose, io)

Extend `knowns` in place by whatever the component's `:obsrules` and the user's init `formulas`
can derive from it, using [`resolve_rules`](@ref).

`targets` narrows down what is worth computing, e.g. an `NWState` only cares about states and
parameters. With `t=NaN` no time is given, and a rule reading `:t` simply does not fire.

See [`_resolution_rules`](@ref) for when the pass runs at all.
"""
function extend_knowns_by_formulas!(knowns, cf, formulas;
                                    am=get_aliasmap(cf), t=NaN, targets=nothing,
                                    error_unresolvable=true, verbose=false, io=stdout)
    # collect rules (user formulas + obs rules)
    rules = _resolution_rules(cf, formulas, am, s -> haskey(knowns, s))
    isempty(rules) && return knowns
    settable = settable_symbols(cf)   # after the guard: nothing above it needs the set
    if isnothing(targets)
        # include outputs of user formulas in targets so no user rule is dropped
        targets = settable ∪ Set(s for r in rules if _is_user_rule(r) for s in r.outsym)
    end

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

Guess-formula sibling of [`extend_knowns_by_formulas!`](@ref), same executor but seeded with both
`defaults` and `guesses`.

Anything the init pass determined is off limits here. A guess formula may overwrite an existing
guess, an observed equation may only fill a symbol nothing else reached. Conflicts are ignored,
guesses are inconsistent by construction.
"""
function extend_guesses_by_formulas!(guesses, defaults, cf, formulas;
                                     am=get_aliasmap(cf), t=NaN, targets=nothing,
                                     error_unresolvable=true, verbose=false, io=stdout)
    # collect rules (user formulas + obs rules) before anything else, so a component with
    # nothing to resolve pays for neither the seed nor the settable set
    rules = _resolution_rules(cf, formulas, am, s -> haskey(guesses, s) || haskey(defaults, s))
    isempty(rules) && return guesses
    settable = settable_symbols(cf)
    isnothing(targets) && (targets = settable)

    # layered seed, `defaults` over `guesses` — the same lookup the guess pass always had
    seed = Dict{Symbol,Float64}(guesses)
    merge!(seed, defaults)

    # the init output is bound, a pre-existing guess is not
    bound = Set{Symbol}(keys(defaults))
    if !isnan(t)
        seed[:t] = t # `t` so that formulas and observables can read it
        push!(bound, :t)
    end
    seedprov = s -> s ∈ bound ? :provided : :guess
    # final targets: prune away any symbol already known and known to be not overwritten
    targets = _prune_resolution_targets(targets, rules, seed; seedprov)
    # actually resolve rules (result stored in res)
    res = resolve_rules(seed, rules; targets, seedprov)
    _check_resolution(res, rules; error_unresolvable, verbose, io,
                      type="GuessFormula", check_conflicts=false)

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
    _resolution_rules(cf, formulas, am, hasvalue) -> Vector{ResolutionRule}

The rule bucket to apply. Init stays cheap when there is nothing to resolve, so the bucket is
non-empty only if

- the user provided a formula, or
- some value sits on an observable — a `default` on the observable `y` of `y ~ -x` reaches the
  state `x` only through the inverse rule.

`hasvalue(s)` answers "does `s` carry a value", asked about the observables rather than the other
way round: this runs per component on every `NWState`, so it must not build anything.

Obs rules go first so a user formula wins the arbitration between two rules ready in the same
pass.
"""
function _resolution_rules(cf, formulas, am, hasvalue)
    hasformulas = !isnothing(formulas) && !isempty(formulas)
    if !hasformulas && !any(hasvalue, obssym(cf))
        return ResolutionRule[]
    end

    obsrules = get_obsrules(cf)
    rules = ResolutionRule[]
    sizehint!(rules, length(obsrules) + (hasformulas ? length(formulas) : 0))
    append!(rules, obsrules)
    if hasformulas
        for f in formulas
            push!(rules, ResolutionRule(f, am))
        end
    end
    rules
end

"""
    _prune_resolution_targets(targets, rules, vals; seedprov) -> Set{Symbol}

Targets are what the rules are pruned for. Narrows `targets` down to:
- the ones written by some rule
- MINUS the symbols which are already in `vals` and outrank the rules writing them

(for example, a user provided default can only be changed by a strong formula, so
we can remove it as target UNLESS there is such a strong formula in the ruleset)

`seedprov` says what the entries of `vals` are worth, same as in `resolve_rules`.
"""
function _prune_resolution_targets(targets, rules, vals; seedprov=Returns(:provided))
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
            held >= rank && continue
        end
        push!(kept, s)
    end
    kept
end

# Decide what the recorded problems mean — `resolve_rules` itself never errors or warns.
# `check_conflicts=false` for the guess pass, guesses are inconsistent by construction.
function _check_resolution(res, rules; error_unresolvable, verbose, io,
                           type="InitFormula", check_conflicts=true)
    _check_pruned(res, rules; error_unresolvable, verbose, io, type)
    _check_skipped(res, rules; verbose, io, type)
    _check_yielded(res, rules; verbose, io, type)
    _check_nonfinite(res, rules; error_unresolvable, verbose, io)

    if check_conflicts && !isempty(res.conflicts)
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
# An optional user rule that never ran. `res.unfired` leaves these out — not running is what
# optional means — but the log should still show that the formula was there and did nothing.
function _check_skipped(res, rules; verbose, io, type)
    verbose || return nothing
    for (i, r) in enumerate(rules)
        (r.optional && _is_user_rule(r) && !res.fired[i] && !res.pruned[i]) || continue
        unknown = [s for s in r.sym if !haskey(res.vals, s)]
        printstyled(io, " - $type: $(_rule_ref(r)) skipped, it is optional and \
                         $unknown never became known\n")
    end
    nothing
end

# A user rule whose write was outranked. Pruning catches the same thing earlier where the rule
# writes *only* outranked targets, but a rule that also feeds something live survives pruning and
# yields at the write instead — reporting both keeps a weak formula's yield visible either way.
function _check_yielded(res, rules; verbose, io, type)
    verbose || return nothing
    seen = Set{Tuple{Int,Symbol}}()
    for y in res.yields
        r = rules[y.rule]
        (_is_user_rule(r) && (y.rule, y.sym) ∉ seen) || continue
        push!(seen, (y.rule, y.sym))
        # the *final* value, which is the one that superseded this rule whichever way round the
        # two fired — `offered` is what this rule wanted, and saying both is the point
        printstyled(io, " - $type: $(_rule_ref(r)) yields on :$(y.sym) = \
                         $(str_significant(res.vals[y.sym]; sigdigits=5)) \
                         (would have set $(str_significant(y.offered; sigdigits=5)))\n")
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
# Complain about user formulas that got pruned away. Obs rules are pruned all the time, that is
# the whole point of `targets`, so they never show up here.
function _check_pruned(res, rules; error_unresolvable, verbose, io, type)
    for (i, r) in enumerate(rules)
        (res.pruned[i] && _is_user_rule(r)) || continue
        # "pruned" covers two very different stories: the rule yielded because its targets are
        # already determined, or nothing reads what it writes
        outranked = all(r.outsym) do s
            haskey(res.vals, s) && _precedence(res.provenance[s]) > _precedence(r.provenance)
        end
        if outranked
            verbose && printstyled(io, " - $type: $(_rule_ref(r)) yields, \
                                        its target $(r.outsym) is already determined\n")
        elseif r.optional || r.provenance === :weak_formula || r.provenance === :guess_formula
            # these are meant to be droppable, so no warning; a library full of them would
            # otherwise drown the log
            verbose && printstyled(io, " - $type: $(_rule_ref(r)) had no effect, \
                                        nothing reads $(r.outsym)\n")
        elseif error_unresolvable # the init path; `NWState` prunes such rules routinely
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
    pin && (note = strip(note * " (obs)"))
    via = isnothing(label) ? "" : "(via $label)"
    ":$s &$op $v &$note &$via"
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
