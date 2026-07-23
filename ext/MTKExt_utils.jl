
"""
    eq_type(eq::Equation)

Checks the type of the equation. Returns:
- `(:explicit_diffeq, lhs_variable)` for explicit differential equations
- `(:implicit_diffeq, nothing)` for implicit differential equations
- `(:explicit_algebraic, lhs_variable)` for explicit algebraic equations
- `(:implicit_algebraic, nothing)` for implicit algebraic equations

"""
function eq_type(eq::Equation)
    rhs_differentials = _collect_differentials(eq.rhs)
    if !isempty(rhs_differentials)
        return (:implicit_diffeq, nothing)
    end

    if isdifferential(eq.lhs)
        return (:explicit_diffeq, only(eq.lhs.args))
    end

    lhsvars = get_variables(eq.lhs)
    rhsvars = get_variables(eq.rhs)
    commonvars = intersect(lhsvars, rhsvars)
    if isempty(commonvars) && length(lhsvars) == 1 && isequal(eq.lhs, only(lhsvars))
        return (:explicit_algebraic, eq.lhs)
    else
        return (:implicit_algebraic, nothing)
    end
end


function rhs_differentials(eqs::Vector{Equation})
    diffs = Set{ST}()
    for eq in eqs
        _collect_differentials!(diffs, eq.rhs)
    end
    return diffs
end
rhs_differentials(eq::Equation) = _collect_differentials!(Set{ST}(), eq.rhs)

_collect_differentials(ex) = _collect_differentials!(Set{ST}(), ex)
_collect_differentials(eq::Equation) = _collect_differentials!(Set{ST}(), eq.lhs+eq.rhs)

function _collect_differentials!(found, ex)
    if iscall(ex)
        if operation(ex) isa Differential
            push!(found, ex)
        else
            for arg in arguments(ex)
                _collect_differentials!(found, arg)
            end
        end
    end
    return found
end

"""
    getproperty_symbolic(sys, var; might_contain_toplevel_ns=true)

Like `getproperty` but works on a greater varaity of "var"
- var can be Num or Symbolic (resolved using genname)
- strip namespace of sys if present (don't strip if `might_contain_top_level_ns=false`)
- for nested variables (foo₊bar₊baz) resolve them one by one
"""
function getproperty_symbolic(sys, var; might_contain_toplevel_ns=true)
    ns = string(getname(sys))
    varname = string(getname(var))
    # split of the toplevel namespace if necessary
    if might_contain_toplevel_ns && startswith(varname, ns*"₊")
        if getname(sys) ∈ getname.(ModelingToolkitBase.get_systems(sys))
            @warn "Namespace :$ns appears multiple times, this might lead to unexpected, since it is not clear whether the namespace should be stripped or not."
        end
        varname = replace(varname, r"^"*ns*"₊" => "")
    end
    parts = split(varname, "₊")
    # Descend a `₊`-segment ONLY when it names an actual child subsystem of the
    # current node; the first segment that is not a subsystem starts a flat variable
    # name, so the remaining segments (re-joined by `₊`) are resolved as one variable.
    # This is necessary because a flattened/transpiled leaf can carry `₊` *inside* a
    # variable name (e.g. `washout₊derivative₊x` is a single unknown, not the chain
    # washout→derivative→x). Driven by `get_systems` (the real structure), so it is
    # not a try/catch guess. First level uses namespace=false (strip the toplevel),
    # subsequent levels namespace=true (accumulate the full name) — as before.
    r = sys
    i = 1
    while i <= length(parts)
        first = (i == 1)
        if r isa System && Symbol(parts[i]) in getname.(ModelingToolkitBase.get_systems(r))
            r = getproperty(r, Symbol(parts[i]); namespace = !first)
            i += 1
        else
            r = getproperty(r, Symbol(join(parts[i:end], "₊")); namespace = !first)
            break
        end
    end
    unwrap(r)
end

function generate_massmatrix(eqs::AbstractVector{Equation})
    V = map(eqs) do eq
        type = eq_type(eq)[1]
        if type === :explicit_diffeq
            1
        elseif type === :implicit_algebraic
            0
        else
            error("Cant build mass matrix entry for $(eq) of type $type")
        end
    end
    M = Diagonal(V)
    return M==I ? I : M
end

function warn_missing_features(sys)
    cev = _collect_continuous_events(sys)
    dev = _collect_discrete_events(sys)
    if !isempty(cev) || !isempty(dev)
        @warn "MTK-Model $(sys.name) contains events. They will be ignored by NetworkDynamics.jl. Use ComponentCallbacks on a component level for that!"
    end

    if !isempty(ModelingToolkitBase.initialization_equations(sys))
        @warn "Model has explicit init equation. Those are currently ignored by NetworkDynamics.jl."
    end

    calls = filter(iscall, parameters(sys) ∪ unknowns(sys))
    any(x -> operation(x) === Base.getindex, calls) && @warn "NetworkDynamics does not support vector-variables or vector-parameters in MTK models. Detected: $(join(repr.(filter(x -> operation(x) === Base.getindex, calls)), ", "))"

    has_discretes = any(unknowns(sys)) do s
        s.metadata[Symbolics.VariableSource][1] == :discretes
    end
    if has_discretes
        @warn "Model contains variables defined via @discretes. This is not supported in NetworkDynamics. \
        For discretely changing variables use normal @parameters with callbacks on the component level."
    end
end
function _collect_continuous_events(sys)
    vcat(
        ModelingToolkitBase.get_continuous_events(sys),
        [_collect_continuous_events(sys) for sys in ModelingToolkitBase.get_systems(sys)]...
    )
end
function _collect_discrete_events(sys)
    vcat(
        ModelingToolkitBase.get_discrete_events(sys),
        [_collect_discrete_events(sys) for sys in ModelingToolkitBase.get_systems(sys)]...
    )
end

function check_metadata(exprs)
    nometadata = []
    for ex in exprs
        if ex isa Equation
            _check_metadata!(nometadata, ex.rhs)
            _check_metadata!(nometadata, ex.lhs)
        else
            _check_metadata!(nometadata, ex)
        end
    end
    return unique!(nometadata)
end
function _check_metadata!(nometadata, expr)
    vars = get_variables_deriv(expr)
    for v in vars
        isnothing(Symbolics.metadata(v)) && push!(nometadata, v)
    end
end

function fix_metadata!(invalid_eqs, sys)
    missingmetadata = check_metadata(invalid_eqs)
    if isempty(missingmetadata)
        return invalid_eqs
    end

    metadatasubs = Dict()
    allsyms = ModelingToolkitBase.all_symbols(sys)
    filter!(s->!contains(repr(s), "Initial"), allsyms)
    allnames = string.(ModelingToolkitBase.getname.(allsyms))

    for invalids in missingmetadata
        invalidname = getname(invalids)
        valid = if hasproperty(sys, getname(invalidname))
            getproperty_symbolic(sys, invalids) # works on namespaced symbols foo₊bar directly
        else
            idxs = findall(contains(string(invalidname)), allnames)
            if length(idxs) == 1
                allsyms[only(idxs)]
            else
                @warn "Could not resolve invalid symbol $invalidname, options are $(allsyms[idxs])"
            end
        end
        metadatasubs[invalids] = valid
    end

    fixedeqs = [Symbolics.substitute(eq, metadatasubs) for eq in invalid_eqs]
    if !isempty(check_metadata(fixedeqs))
        @warn "Some transformation droped metadata ($missingmetadata)! Could NOT be fixed. $(check_metadata(fixedeqs))"
    else
        @warn "Some transformation droped metadata ($missingmetadata)! Could be fixed."
    end
    invalid_eqs .= fixedeqs
end

function remove_implicit_output_fn!(eqs)
    r = SymbolicUtils.@rule implicit_output(~~x) => 0
    chain = SymbolicUtils.Chain([r])
    rewriter = SymbolicUtils.Prewalk(chain)

    for i in eachindex(eqs)
        eq = eqs[i]
        eqs[i] = rewriter(eq.lhs) ~ rewriter(eq.rhs)
    end

    eqs
end

function collect_comp_metadata(sys, key)
    md = OrderedDict{String, Any}()
    _collect_comp_metadata!(md, sys, key, nothing)
    md
end
function _collect_comp_metadata!(md, sys, key, current_namespace)
    # get namespaced model name
    ns_modelname = if isnothing(current_namespace)
        "" # toplevel needs no name
    elseif isempty(current_namespace)
        string(getname(sys))
    else
        string(current_namespace, "₊", getname(sys))
    end

    # collect from subsystems first (children before parents)
    for s in ModelingToolkitBase.get_systems(sys)
        _collect_comp_metadata!(md, s, key, ns_modelname)
    end

    # collect from current system
    thismd = ModelingToolkitBase.get_metadata(sys)
    if haskey(thismd, key)
        if haskey(md, ns_modelname) && md[ns_modelname] != thismd[key]
            _modname = isempty(ns_modelname) ? "toplevel model" : ns_modelname
            @warn "Overwriting metadata $key for $_modname while collecting"
        end
        md[ns_modelname] = thismd[key]
    end
    nothing
end


function apply_component_postprocessing!(cf)
    sys = cf.metadata[:odesystem]
    md = collect_comp_metadata(sys, ComponentPostprocessing)
    for (modelname, ppfs) in md
        if ppfs isa Union{AbstractVector, Tuple}
            for ppf in ppfs
                _check_symbol(ppf, modelname)
                ppf(cf, modelname)
            end
        else
            _check_symbol(ppfs, modelname)
            ppfs(cf, modelname)
        end
    end
end
_check_symbol(x, _) = nothing
function _check_symbol(ppf::Symbol, name)
    _name = isempty(name) ? "toplevel model" : "subsystem \"$name\""
    error(
        """
        The postprocessing function `$ppf` included in $_name could not be captured. \
        This is a limitation of the `@mtkmodel` macro. You need to make sure, that the function \
        $ppf is defined before the model definition! I.e. write

        function $ppf end
        @mtkmodel ModelUsedAs begin
            @metadata begin
                ComponentPostprocessing = $ppf
            end
        end

        or put the entire function above the model definition.
        """
    )
end

# docstring lives in utils.jl
function NetworkDynamics.set_mtk_defaults(sys::System, pairs)
    isempty(pairs) && return sys
    defs = values(pairs)
    names = keys(pairs)
    symbols = map(names) do name
        try
            getproperty_symbolic(sys, name; might_contain_toplevel_ns=false)
        catch e
            throw(ArgumentError("Could not resolve variable of name $name in system $(getname(sys))! None of the defaults have been set."))
        end
    end
    # both dicts are copies: `get_initial_conditions` hands back the live dict, and `sys`
    # must come out of here untouched — every change leaves on the returned system
    defdict = copy(ModelingToolkitBase.get_initial_conditions(sys))
    binddict = ModelingToolkitBase.SymmapT()
    for (k, v) in ModelingToolkitBase.get_bindings(sys)
        binddict[k] = v
    end

    # Mirror MTK's `collect_defaults!`, which routes a variable's `default` metadata to
    # either the initial conditions (numeric) or the bindings (symbolic) at `System`
    # construction. Values forwarded through this function arrive too late for that pass,
    # so they are classified here instead.
    bindings_changed = false
    for (s, v) in zip(symbols, defs)
        us = unwrap(s)
        if SII.symbolic_type(v) === SII.NotSymbolic()
            defdict[s] = v
            if haskey(binddict, us)
                delete!(binddict, us)
                bindings_changed = true
            end
        else
            binddict[us] = unwrap(v)
            delete!(defdict, us)
            bindings_changed = true
        end
    end

    @reset sys.initial_conditions = defdict
    bindings_changed || return sys

    @reset sys.bindings = ModelingToolkitBase.ROSymmapT(binddict)
    ## a cached bindings graph (built by `complete`) would be stale now
    @reset sys.parameter_bindings_graph = nothing
    sys
end

# docstring lives in utils.jl
function NetworkDynamics.set_initf(sys::System, pairs::Pair...)
    isempty(pairs) && return sys
    for (target, _) in pairs
        u = unwrap(target)
        vars = get_variables(u)
        if length(vars) != 1 || !isequal(only(vars), u)
            throw(ArgumentError("set_initf target $target is not a single variable."))
        end
    end
    existing = SymbolicUtils.getmetadata(sys, NetworkDynamics.SystemInitFormulas, Pair[])
    combined = vcat(existing, [unwrap(t) => unwrap(e) for (t, e) in pairs])
    SymbolicUtils.setmetadata(sys, NetworkDynamics.SystemInitFormulas, combined)
end

# docstring lives in utils.jl
function NetworkDynamics.set_guessf(sys::System, pairs::Pair...)
    isempty(pairs) && return sys
    for (target, _) in pairs
        u = unwrap(target)
        vars = get_variables(u)
        if length(vars) != 1 || !isequal(only(vars), u)
            throw(ArgumentError("set_guessf target $target is not a single variable."))
        end
    end
    existing = SymbolicUtils.getmetadata(sys, NetworkDynamics.SystemGuessFormulas, Pair[])
    combined = vcat(existing, [unwrap(t) => unwrap(e) for (t, e) in pairs])
    SymbolicUtils.setmetadata(sys, NetworkDynamics.SystemGuessFormulas, combined)
end

"""
    match_diff_states(eqs, states::Set)

Match a set of unknowns to a vector of equations.
Returns a vector of states with diff states at the correct places
and non-diff states filled in arbitrarily.
"""
function match_diff_states(eqs, states)
    @assert length(eqs) == length(states) "Number of equations and states must match to be able to match them, got $(length(eqs)) equations and $(length(states)) states"
    states_new = Vector{eltype(states)}(undef, length(eqs))
    diff_state_set = Set(states) # for fast membership test
    for (i, eq) in pairs(eqs)
        type, var = eq_type(eq)
        if type == :explicit_diffeq
            states_new[i] = pop!(diff_state_set, var)
        elseif type == :implicit_algebraic
            # do nothing
        else
            error("Got equation $eq of type $type, which is not expected at that stage (only explicit diff eqs and implicit alg. eqs). Something went wrong.")
        end
    end
    # Fill algebraic slots in the order states appear in the input vector (not hash order).
    alg_states = [s for s in states if s ∈ diff_state_set]
    ai = 0
    for i in eachindex(states_new)
        if !isassigned(states_new, i)
            ai += 1
            states_new[i] = alg_states[ai]
        end
    end
    @assert ai == length(alg_states) "Not all states could be matched, something went wrong"
    @assert length(eqs) == length(states_new)
    states_new
end

# Like `get_variables` but replaces any `D(x)` terms with the inner variable `x`.
# `get_variables` intentionally treats `D(x)` as an atomic symbol (expected behavior),
# but we often need the bare variable to check dependencies in the equation graph.
function get_variables_deriv(ex)
    set = get_variables(ex)
    for s in set
        @match s begin
            SymbolicUtils.BSImpl.Term(; f=Differential(_, 1), args=symvec) => begin
                pop!(set, s)
                push!(set, only(symvec))
            end
            _ => nothing
        end
    end
    set
end

# A value is a *symbolic constant* iff it has no free variables: a plain `Number`, or a
# `BasicSymbolic` that folds to one (e.g. an MTKv11 literal `0.0`, or `sqrt(2)`). A value that
# references other variables/parameters is NOT constant; only constants are stored as plain
# `:guess` metadata. A non-constant `:guess` is rejected by `_get_metadata` — the `guessf`
# variable option (or `set_guessf`) is the way to spell a symbolic guess.
is_symbolic_constant(x) = (u = unwrap(x); u isa Number || isempty(get_variables(u)))

# Numeric value of a symbolic constant; only meaningful when `is_symbolic_constant(x)`.
symbolic_constant_value(x) = (u = unwrap(x); u isa Number ? u : Float64(Symbolics.value(u)))

vigformula_docstring = raw"""
    VariableInitFormula
    VariableGuessFormula

Metadata types behind the `initf` / `guessf` variable options. A symbol carrying

    @variables x(t) [initf  = <expression>]     # VariableInitFormula
    @variables x(t) [guessf = <expression>]     # VariableGuessFormula

declares "at initialization, set (`initf`) / guess (`guessf`) `x` from `<expression>`". The
expression may reference any unknown, parameter or observable of the system; it is lowered to
an [`InitFormula`](@ref) / [`GuessFormula`](@ref) at compile time, with the names exactly as
written — classification and observable expansion happen at init time. Valid on unknowns and
parameters alike; on a variable which simplification turns into a non-alias observable, the
formula *pins* that observable (see the initialization docs).

`initf` is a constraint (lands in defaults, consistency-checked); `guessf` is only a hint
(lands in guesses, never checked, skipped when its inputs cannot be resolved — leaving any
scalar `guess` as the fallback). For targets belonging to a subsystem use [`set_initf`](@ref)
/ [`set_guessf`](@ref).
"""


@doc vigformula_docstring
struct VariableInitFormula end

@doc vigformula_docstring
struct VariableGuessFormula end

Symbolics.option_to_metadata_type(::Val{:initf})  = VariableInitFormula
Symbolics.option_to_metadata_type(::Val{:guessf}) = VariableGuessFormula

"""
    ParameterBoundTo

Symbolic metadata type backing `@parameters S_b [bound_to = :busbar₊S_b]`. Declares that a
parameter is a structural alias of another symbol in the *same* component. Realized as a real
MTK binding in the `VertexModel`/`EdgeModel` constructor, before compilation (see
[`resolve_bound_to`](@ref)): the bound parameter is demoted to an observable of its target and
leaves `psym`, so there is exactly one true parameter for the quantity.

The metadata value is the target's full (namespaced with `₊`) symbol name. `bound_to` is
component-local — the structural binding can only be injected before compilation, so the target
must live in the same component; cross-component defaulting is `default_from` instead. An
explicit default on a `bound_to` parameter, or an unresolvable target, is an error.
"""
struct ParameterBoundTo end

Symbolics.option_to_metadata_type(::Val{:bound_to}) = ParameterBoundTo

"""
    resolve_bound_to(sys) -> System

Turn `bound_to` parameter metadata into MTK bindings before compilation. For each parameter
`@parameters p [bound_to = :target]`, inject the binding `p = target` (via
[`set_mtk_defaults`](@ref)), so that `complete`/`mtkcompile` demotes `p` to an observable while
`target` stays the parameter. Returns `sys` unchanged if no parameter carries `bound_to`.

Errors if a `bound_to` parameter also has an explicit default (contradictory intent) or if the
target cannot be resolved in the component. See [`ParameterBoundTo`](@ref).
"""
function resolve_bound_to(sys)
    ps = parameters(sys)
    pset = Set(unwrap.(ps))
    defs = initial_conditions(sys)
    bindings = Pair{Symbol,Any}[]
    # every symbol that could carry `bound_to` metadata — parameters, states, observed. It is
    # only meaningful on a parameter (it becomes an MTK binding); anywhere else is a mistake.
    candidates = Iterators.flatten((ps, unknowns(sys), (eq.lhs for eq in observed(sys))))
    for v in candidates
        uv = unwrap(v)
        target = Symbolics.getmetadata(uv, ParameterBoundTo, nothing)
        isnothing(target) && continue
        if uv ∉ pset
            throw(ArgumentError(
                "`bound_to` on `$(getname(v))` is invalid: it only applies to parameters, but \
                 `$(getname(v))` is not a parameter of component `$(getname(sys))`."))
        end
        tname = target isa QuoteNode ? target.value : target
        if haskey(defs, v)
            throw(ArgumentError(
                "Parameter `$(getname(v))` has both an explicit default (`$(defs[v])`) and \
                 `bound_to = $tname`. These are contradictory: a `bound_to` parameter is \
                 eliminated from the parameters. Remove the default or the binding."))
        end
        tsym = try
            getproperty_symbolic(sys, tname; might_contain_toplevel_ns=false)
        catch
            throw(ArgumentError(
                "`bound_to = $tname` on parameter `$(getname(v))` could not be resolved in \
                 component `$(getname(sys))`. Available parameters: $(sort(getname.(ps)))."))
        end
        push!(bindings, getname(v) => tsym)
    end
    isempty(bindings) && return sys
    NetworkDynamics.set_mtk_defaults(sys, (; bindings...))
end

"""
    collect_initf(sys)
    collect_guessf(sys)

Collect the `initf`/`guessf` variable metadata and the [`set_initf`](@ref)/[`set_guessf`](@ref)
system metadata of `sys` and all its subsystems into a list of `target => expression` pairs,
namespaced to the level of `sys`. Deliberately a list, not a dict: a target carrying both a
variable-level and a system-level recipe must surface as two entries, so `_dedupe_resolved`
can dedupe them when identical and error/warn when they conflict — never silently prefer one.

Must be called on the **hierarchical** (pre-flattening) system: `renamespace` renames a
symbol but does not descend into the expressions stored in its metadata, so the formula of a
subsystem symbol read off a flattened system still refers to the subsystem's own symbols by
their bare names. Recursing here and namespacing each collected pair with `namespace_expr`
(which honors `ParentScope`, so expressions passed in from the parent survive untouched) is
the same strategy MTK's own `bindings`/`guesses` accessors use. `set_initf`/`set_guessf` pairs
are written in the local names of the system they were attached to, so the same recursion
namespaces them correctly too.
"""
collect_initf(sys)  = _collect_formula_metadata(sys, VariableInitFormula,  NetworkDynamics.SystemInitFormulas)
collect_guessf(sys) = _collect_formula_metadata(sys, VariableGuessFormula, NetworkDynamics.SystemGuessFormulas)

function _collect_formula_metadata(sys, VarMetaType, SysMetaKey)
    pairs = Pair{ST,Any}[]
    # a subsystem variable referenced in this level's equations shows up in this level's
    # unknowns too, still carrying its metadata — but with the expression in the *subsystem's*
    # local names. Only the recursion below sees it in the right namespace, so it must be
    # skipped here.
    subvars = Set{ST}()
    for subsys in ModelingToolkitBase.get_systems(sys)
        for v in vcat(ModelingToolkitBase.get_unknowns(subsys), ModelingToolkitBase.get_ps(subsys))
            push!(subvars, unwrap(ModelingToolkitBase.namespace_expr(v, subsys)))
        end
    end
    for v in vcat(ModelingToolkitBase.get_unknowns(sys), ModelingToolkitBase.get_ps(sys))
        u = unwrap(v)
        u ∈ subvars && continue
        SymbolicUtils.hasmetadata(u, VarMetaType) || continue
        push!(pairs, u => unwrap(SymbolicUtils.getmetadata(u, VarMetaType)))
    end
    for (target, expr) in SymbolicUtils.getmetadata(sys, SysMetaKey, Pair[])
        push!(pairs, unwrap(target) => unwrap(expr))
    end
    for subsys in ModelingToolkitBase.get_systems(sys)
        for (target, expr) in _collect_formula_metadata(subsys, VarMetaType, SysMetaKey)
            push!(pairs, ModelingToolkitBase.namespace_expr(target, subsys) =>
                         ModelingToolkitBase.namespace_expr(expr, subsys))
        end
    end
    pairs
end

"""
    initf_to_initformulas(pairs)
    guessf_to_guessformulas(pairs)

Turn the `target => expression` pairs collected by
[`collect_initf`](@ref)/[`collect_guessf`](@ref) into InitFormulas/GuessFormulas. Targets may
be unknowns, parameters, inputs — or observables, in which case the formula *pins* the
observable as an init-time dataflow node.

An InitFormula is a constraint (not a hint), so conflicting definitions for the same raw
target are a genuine over-determination and raise an error. A GuessFormula is only a
convergence hint, so conflicting definitions are deduped with a warning, never fatal.
"""
initf_to_initformulas(pairs)   = _metadata_to_formulas(pairs, InitFormula;  fail=:error, kind="initf")
guessf_to_guessformulas(pairs) = _metadata_to_formulas(pairs, GuessFormula; fail=:warn,  kind="guessf")

function _metadata_to_formulas(pairs, ::Type{FT}; fail::Symbol, kind::String) where {FT}
    (isnothing(pairs) || isempty(pairs)) && return nothing

    resolved = Any[]
    for (target, expr) in pairs
        r = _resolve_formula(target, expr; kind)
        r === nothing && continue
        push!(resolved, r)
    end

    resolved = _dedupe_resolved(resolved; fail, kind=string(nameof(FT)))
    isempty(resolved) ? nothing : [_build_formula(FT, r) for r in resolved]
end

"""
    assert_no_state_bindings(sys)

Throw if `sys` binds any unknown to an expression, i.e. if it was built with
`bindings = [x => expr]` or, equivalently, `@variables x(t) = expr` for a symbolic `expr`.

Parameter bindings are fine and are left alone: for a parameter, a symbolic default is a
genuine runtime dependency which MTK lowers to an observed equation, removing the parameter
from the compiled model. Only for an unknown is the intent ambiguous, and there it is spelled
`initf`.

Runs on the **hierarchical** system, so it reports only bindings the user actually wrote
(`mtkcompile` emits bindings of its own). Note this rules out `bound_parameters` for telling
the two kinds apart — that one needs a completed system — hence the `isparameter` check.
"""
function assert_no_state_bindings(sys)
    bindings = ModelingToolkitBase.bindings(sys)
    isempty(bindings) && return nothing

    offenders = [target => expr for (target, expr) in bindings
                 if !ismissing(unwrap_const(expr)) &&
                    !ModelingToolkitBase.isparameter(unwrap(target))]
    isempty(offenders) && return nothing

    list = join(("  $target = $expr" for (target, expr) in offenders), "\n")
    rewrite = join(("  @variables $target [initf = $expr]" for (target, expr) in offenders), "\n")
    throw(ArgumentError(
        """
        System :$(getname(sys)) binds unknown(s) to an expression:
        $list
        A binding on an unknown (`bindings = [x => expr]`, or `@variables x(t) = expr` with a \
        symbolic `expr`) does not say when the expression should hold. Declare it as an explicit \
        initialization equation instead:
        $rewrite
        which sets the target once, during initialization, and leaves the dynamics alone.

        Note this does not apply to *parameters*: a symbolic default on a parameter is a runtime \
        dependency which MTK lowers to an observed equation, shadowing the parameter away. That \
        keeps working and is often what you want.
        """))
end

# ── shared helpers for _metadata_to_formulas (initf and guessf) ──
#
# Both an InitFormula and a GuessFormula are "set `target := f(inputs)`" objects, so the
# resolution, conflict-deduplication and function-building are identical; only the conflict
# policy (`fail`: error for initf, warn for guessf) differs between the two.

# Shape one `lhs => rhs` pair into a `(target, rhs, inputs)` form, or skip it (with a
# warning) when it structurally cannot become a formula. Deliberately *raw*: target and
# inputs keep the names the user wrote, so classification and observable expansion happen at
# init time in `normalize` — one resolution path, shared with hand-attached formulas.
# Existence of the raw names is checked at attach time (`add_initformula_lenient!` etc.).
function _resolve_formula(lhs_sym, rhs_expr; kind)
    lhs_vars = get_variables(lhs_sym)
    if length(lhs_vars) != 1 || !isequal(only(lhs_vars), unwrap(lhs_sym))
        @warn "$kind target $lhs_sym is not a single variable. Skip."
        return nothing
    end
    target = getname(lhs_sym)

    input_symbolic = collect(get_variables(rhs_expr))
    input_names    = Symbol[getname(s) for s in input_symbolic]
    # no raw self-dependency (the formula constructors throw on it); a dependency hidden
    # behind an observable is only detectable at init time, where `normalize` reports it
    if target ∈ input_names
        @warn "$kind for $lhs_sym depends on its own target. Skip."
        return nothing
    end

    (; src=lhs_sym, target, rhs=unwrap(rhs_expr), input_symbolic, input_names)
end

"""
    add_initformula_lenient!(c, formula)
    add_guessformula_lenient!(c, formula)

Attach an MTK-lowered formula to a compiled component, demoting the compatibility check to a
warn-and-skip. Lowering is best effort: a recipe whose target or inputs did not survive
simplification (or which references the independent variable directly) must not fail
compilation — the model works without it, the initialization just knows less. Only the
compatibility check is caught; the attach itself runs unchecked, having already been
validated.
"""
function add_initformula_lenient!(c, formula)
    try
        assert_initformula_compat(c, formula)
    catch e
        e isa ArgumentError || rethrow()
        @warn "Skipping an InitFormula which does not fit the compiled component :$(c.name): $(e.msg)" formula
        return nothing
    end
    add_initformula!(c, formula; check=false)  # already validated
    nothing
end
function add_guessformula_lenient!(c, formula)
    try
        assert_guessformula_compat(c, formula)
    catch e
        e isa ArgumentError || rethrow()
        @warn "Skipping a GuessFormula which does not fit the compiled component :$(c.name): $(e.msg)" formula
        return nothing
    end
    add_guessformula!(c, formula; check=false)  # already validated
    nothing
end

# Collapse entries with the same raw target. Identical definitions (equal rhs) are always
# safe to dedupe silently. Genuinely *conflicting* definitions (same target, differing rhs)
# are a warning for guesses (`fail=:warn`, only a hint) but an error for initf
# (`fail=:error`, two constraints forcing one state). Comparison is on the symbolic rhs, not
# the prettyprint. Note this is syntactic and name-level: two formulas targeting different
# members of one alias class both survive here, and the init-time duplicate-writer check
# reports them once normalization has collapsed the class.
function _dedupe_resolved(resolved; fail::Symbol, kind)
    by_target = OrderedDict{Symbol,Vector{Any}}()
    for r in resolved
        push!(get!(() -> Any[], by_target, r.target), r)
    end
    kept = Any[]
    for (target, group) in by_target
        chosen = first(sort(group; by = g -> repr(g.rhs)))
        push!(kept, chosen)
        length(group) == 1 && continue
        if all(g -> isequal(g.rhs, chosen.rhs), group)
            @debug "$kind: $(length(group)) identical definitions for $target; keeping one."
        else
            msg = "$kind: conflicting definitions target $target (differing right-hand \
                   sides). Keeping $(repr(chosen.rhs)); dropping the rest."
            fail === :error ? error(msg) : @warn msg
        end
    end
    kept
end

# Build the actual Init/GuessFormula from a resolved entry. Identical for both formula types
# apart from the type name. The prettyprint mirrors the macro form a hand-written formula
# prints as (see `_macro_source_string`), even though these come from `initf`/guess metadata
# rather than an `@initformula` call — it reads nicer and keeps the two origins consistent.
function _build_formula(::Type{FT}, r) where {FT}
    macroname = "@" * lowercase(string(nameof(FT)))
    f = Symbolics.build_function([r.rhs], r.input_symbolic; expression=Val(false))[2]

    # spell the inputs as the `:sym` the macro form uses, in place of the symbolic variables
    rhsstring = repr(r.rhs)
    for input in r.input_symbolic
        rhsstring = replace(rhsstring, repr(input) => repr(getname(input)))
    end
    # `repr` prints a numeric coefficient times a variable as juxtaposition (`3x`); with the
    # variable now a `:sym` that would read as a range (`3:x`), so restore the explicit `*`.
    rhsstring = replace(rhsstring, r"(?<=[0-9.]):(?=[A-Za-z_])" => " * :")
    prettyprint = "$macroname begin\n    $(repr(r.target)) = $(rhsstring)\nend"
    FT(f, [r.target], r.input_names, prettyprint)
end

"""
    extract_aliasmap(c::ComponentModel, obseqs::Vector{Equation})

Extracts the `AliasMap` of a compiled component from its observed equations, i.e.
every observable which is a pure `factor * symbol` alias.

Note that `pick_best_alias_names` already consolidates *identity* alias groups onto a single
representative and re-inserts direct `alias ~ main` observations; what remains for us are the
scaled/sign-flipped aliases (`get_alias` matches no coefficient) plus chains through them.

A chain either reaches a settable symbol — the usual case — or bottoms out on an observable
whose own equation is not a pure alias (a sum, say). That terminal observable is a valid
canonical too: it is where `generate_obs_expansion` bottoms out anyway, so recording it
unifies the names of a class whose shared value has no storage slot at all. See
[`AliasMap`](@ref) for what that buys and what it costs.

Anything that is not a pure scaled alias — affine (`x + 1`), sums, nonlinear terms, symbolic
coefficients — stays an ordinary observable.
"""
function extract_aliasmap(c::NetworkDynamics.ComponentModel, obseqs)
    settable = settable_symbols(c)
    obs = NetworkDynamics.obssym(c)

    # one-step links only; resolved transitively below
    steps = Dict{Symbol,Tuple{Float64,Symbol}}()
    for eq in obseqs
        m = _match_scaled_var(eq.rhs)
        isnothing(m) && continue
        steps[getname(eq.lhs)] = m
    end

    am = AliasMap()
    for alias in keys(steps)
        # a settable symbol must never be recorded as an alias of another settable symbol;
        # keep both un-aliased instead (`assert_aliasmap_compat` would reject the entry)
        if alias ∈ settable
            @debug "Not aliasing :$alias, it is a settable symbol of the component."
            continue
        end
        factor, canonical = _resolve_alias(alias, steps, settable, Symbol[])
        # `a ~ 2*t` matches as a scaled alias of the independent variable, which is neither
        # settable nor an observable — no canonical, hence no entry
        if canonical ∉ settable && canonical ∉ obs
            @debug "Not aliasing :$alias, its root :$canonical is neither settable nor an observable."
            continue
        end
        am[alias] = (factor, canonical)
    end
    am
end

# Match `ex` against `factor * var` with a numeric, nonzero `factor` and a single variable.
# Returns `(factor, varname)` or `nothing`. Accepts exactly the shape `linear_expansion`
# reports as `factor*v + offset` with a single variable, numeric nonzero factor and zero
# offset; everything else (multi-variable, nonlinear, symbolic/affine, `0*x`) is rejected by
# the guards below. `linear_expansion` returns MTKv11 symbolic literals, so factor and offset
# go through `unwrap_const` and must land on a `Number` to count as numeric.
# Differentials need no handling: `generate_io_function` rejects rhs differentials upfront.
function _match_scaled_var(ex)
    vars = get_variables(ex)
    length(vars) == 1 || return nothing
    v = only(vars)

    _factor, _offset, islinear = Symbolics.linear_expansion(ex, v)
    islinear || return nothing
    factor = unwrap_const(_factor)
    offset = unwrap_const(_offset)
    (factor isa Number && offset isa Number) || return nothing
    iszero(offset) || return nothing
    (isfinite(factor) && !iszero(factor)) || return nothing

    (Float64(factor), getname(v))
end

# Follow one-step links until the chain ends, multiplying factors along the way
# (`obs2 ~ -obs1`, `obs1 ~ 2x` ⇒ `(-2.0, :x)`). It ends on the first settable symbol, or else
# on the first symbol without a further link — the terminal observable, which is the
# canonical of a class that has no settable member. Only ever called for an `s` which has a
# link, so it always resolves.
function _resolve_alias(s, steps, settable, visiting)
    if s ∈ visiting
        error("Cyclic alias chain detected at :$s via $(join(visiting, " → ")). \
               Observed equations are topologically sorted, this should never happen.")
    end
    factor, target = steps[s]
    (target ∈ settable || !haskey(steps, target)) && return (factor, target)

    push!(visiting, s)
    rest = _resolve_alias(target, steps, settable, visiting)
    pop!(visiting)
    (factor * rest[1], rest[2])
end

"""
    generate_obs_expansion(cf::ComponentModel, syms::Vector{Symbol}; stop_at) -> (roots, f)

MTK implementation of the core stub; see `NetworkDynamics.generate_obs_expansion` for the
contract. Expansion is a plain `fixpoint_sub` against the stored `:observed` equations:
they are acyclic, and the output-defining equations are absent from them, so substituting
to a fixpoint necessarily bottoms out on settable symbols (plus the independent variable).
Equations whose lhs name is in `stop_at` are excluded from the substitution set, so the
fixpoint additionally bottoms out on those symbols.

The independent variable is split off from the roots and passed to the closure separately —
callers source root values from the defaults/guesses dicts, where a `t` has no business
being. Symbols with no observed equation never enter the symbolic part at all; they are
copied through by index, which spares us reconstructing a symbolic variable for them.
"""
function NetworkDynamics.generate_obs_expansion(cf::NetworkDynamics.ComponentModel, syms::Vector{Symbol};
                                                stop_at=Set{Symbol}())
    if !NetworkDynamics.has_metadata(cf, :observed)
        throw(ArgumentError("Cannot expand $syms: component :$(cf.name) has no `:observed` \
                             metadata. Only components compiled from ModelingToolkit carry \
                             the symbolic observed equations needed for expansion."))
    end
    obseqs = NetworkDynamics.get_metadata(cf, :observed)
    obseqs = filter(eq -> getname(eq.lhs) ∉ stop_at, obseqs)
    subs = OrderedDict(eq.lhs => eq.rhs for eq in obseqs)
    byname = Dict(getname(eq.lhs) => eq.lhs for eq in obseqs)

    obsidx = findall(s -> haskey(byname, s), syms)
    plainidx = findall(s -> !haskey(byname, s), syms)
    exprs = [fixpoint_sub(subs[byname[syms[i]]], subs) for i in obsidx]

    iv = only(independent_variables(NetworkDynamics.get_metadata(cf, :odesystem_simplified)))
    # `get_variables` hands back a set, so collect before flattening
    symroots = unique(reduce(vcat, [collect(get_variables(ex)) for ex in exprs]; init=[]))
    filter!(v -> !isequal(v, iv), symroots)
    rootnames = getname.(symroots)

    _warn_unsettable_roots(cf, rootnames, stop_at)

    # plain symbols are their own roots; `unique` because a formula may well read both a
    # settable symbol and an alias of it
    roots = unique(vcat(rootnames, syms[plainidx]))
    pos = Dict(s => i for (i, s) in enumerate(roots))
    symrootpos = [pos[n] for n in rootnames]
    plainpos = [pos[syms[i]] for i in plainidx]

    g = isempty(exprs) ? nothing : build_function(exprs, symroots, iv; expression=Val(false))[1]
    n = length(syms)
    expand = function (rootvals, t)
        out = Vector{Float64}(undef, n)
        isnothing(g) || (out[obsidx] .= g(view(rootvals, symrootpos), t))
        for (k, i) in enumerate(plainidx)
            out[i] = rootvals[plainpos[k]]
        end
        out
    end
    (roots, expand)
end

# Mirrors the `obs_deps ⊆ params ∪ inputs ∪ unknowns` warning in `generate_io_function`: a
# root that is not settable cannot be sourced from the defaults dict, so the formula reading
# it will skip. Warn rather than throw — the skip is already diagnosed downstream. Symbols
# the expansion was asked to stop at are deliberate roots, not accidents — never warn on them.
function _warn_unsettable_roots(cf, rootnames, stop_at)
    unsettable = setdiff(rootnames, settable_symbols(cf), stop_at)
    isempty(unsettable) && return nothing
    @warn "Observable expansion for :$(cf.name) bottomed out on non-settable symbol(s) \
           $(collect(unsettable)). Formulas reading them will be skipped."
end

function NetworkDynamics.multiline_repr(eqs::Vector{Equation}; prefix="")
    lines = map(eqs) do eq
        prefix * repr(eq.lhs) * " &~ " * repr(eq.rhs)
    end
    join(NetworkDynamics.align_strings(lines), "\n")
end

function _compare_mtkcompile(VEModel, args, kwargs)
    m1 = VEModel(args...; mtkcompile=true, kwargs...)
    m2 = VEModel(args...; mtkcompile=false, kwargs...)

    # bench for .5 second
    maxiter = 5
    t1min = typemax(Float64)
    iter = 0
    start = time()
    while time() - start < 1 && iter < maxiter
        t = @elapsed VEModel(args...; mtkcompile=true, kwargs..., verbose=false)
        t1min = min(t1min, t)
        iter += 1
    end
    t2min = typemax(Float64)
    iter = 0
    start = time()
    while time() - start < 1 && iter < maxiter
        t = @elapsed VEModel(args...; mtkcompile=false, kwargs..., verbose=false)
        t2min = min(t2min, t)
        iter += 1
    end

    printstyled("Timings:", color=:blue, bold=true)
    print("\n  with mtk:        ", NetworkDynamics.str_significant(t1min; sigdigits=4), " seconds")
    print("\n  witout mtk:      ", NetworkDynamics.str_significant(t2min; sigdigits=4), " seconds")
    print("\n  speedup without: ")
    printstyled(NetworkDynamics.str_significant(t2min/t1min; sigdigits=3)*"x", color=(t2min<t1min) ? :green : :red, bold=true)
    println()

    printstyled("States:", color=:blue, bold=true)
    if NetworkDynamics.insym(m1) == NetworkDynamics.insym(m2)
        printstyled("\n  Insym match! ", color=:green, bold=true)
    else
        printstyled("\n  Insym dont match! ", color=:red, bold=true)
    end
    if NetworkDynamics.outsym(m1) == NetworkDynamics.outsym(m2)
        printstyled("\n  Outsym match! ", color=:green, bold=true)
    else
        printstyled("\n  Outsym dont match! ", color=:red, bold=true)
    end
    samesym = Set(NetworkDynamics.sym(m1)) == Set(NetworkDynamics.sym(m2))
    if samesym
        printstyled("\n  Syms match! ", color=:green, bold=true)
    else
        printstyled("\n  Syms dont match! ", color=:red, bold=true)
        dim1 = NetworkDynamics.dim(m1)
        dim2 = NetworkDynamics.dim(m2)
        if dim1 == dim2
            println("\n  Both have $(dim1) states")
        else
            print("\n    with MTK:    ")
            printstyled("$(dim1) states", color=(dim2>dim1) ? :green : :red)
            print("\n    without MTK: ")
            printstyled("$(dim2) states", color=(dim1>dim2) ? :green : :red)
        end
        identical = collect(NetworkDynamics.sym(m1) ∩ NetworkDynamics.sym(m2))
        print("\n    Identical syms:  ", inline_repr(identical))
        m1_specific = collect(setdiff(NetworkDynamics.sym(m1), identical))
        print("\n    Only with mtk:   ", inline_repr(m1_specific))
        m2_specific = collect(setdiff(NetworkDynamics.sym(m2), identical))
        print("\n    Only witout mtk: ", inline_repr(m2_specific))
    end
    samepsyms =  Set(NetworkDynamics.psym(m1)) == Set(NetworkDynamics.psym(m2))
    if samepsyms
        printstyled("\n  Psyms match! ", color=:green, bold=true)
    else
        printstyled("\n  Psyms dont match! ", color=:red, bold=true)
    end
    !isnothing(m1.extin) || !isnothing(m2.extin) && error("Comparison for extin not supported.")

    # Numerical comparison: seed the model with more states from rand inputs,
    # then reconstruct the smaller model's inputs using m_big's state + observed values.
    printstyled("\nComparison of f & g outputs with random inputs: ", color=:blue, bold=true)
    if !samepsyms
        printstyled("\n  Skipped: psyms don't match", color=:yellow, bold=true)
    else
        m_big, m_small = NetworkDynamics.dim(m1) >= NetworkDynamics.dim(m2) ? (m1, m2) : (m2, m1)
        psym_perm = [findfirst(==(p), NetworkDynamics.psym(m_small)) for p in NetworkDynamics.psym(m_big)]

        # Seed from m_small: its states are all "real" (no algebraically lifted extras),
        # so the result is consistent without any post-hoc fixup.
        outs_small, du_small, u_small, ins, p_small, t = NetworkDynamics.rand_inputs_fg(m_small)
        NetworkDynamics.compfg(m_small)(outs_small, du_small, u_small, ins, p_small, t)

        # collect m_small's observed and output values for filling m_big
        obsout_small = zeros(length(NetworkDynamics.obssym(m_small)))
        if !isempty(obsout_small) && !isnothing(m_small.obsf)
            m_small.obsf(obsout_small, u_small, ins..., p_small, t)
        end
        outs_small_flat  = reduce(vcat, outs_small)
        outsym_small_flat = reduce(vcat, NetworkDynamics.outsym(m_small))

        # for each state of m_big: look it up in m_small.sym, obssym, then outsym
        u_big_raw = map(NetworkDynamics.sym(m_big)) do s
            idx = findfirst(==(s), NetworkDynamics.sym(m_small))
            !isnothing(idx) && return u_small[idx]
            idx = findfirst(==(s), NetworkDynamics.obssym(m_small))
            !isnothing(idx) && return obsout_small[idx]
            idx = findfirst(==(s), outsym_small_flat)
            !isnothing(idx) && return outs_small_flat[idx]
            return nothing
        end

        n_missing = count(isnothing, u_big_raw)
        if n_missing > 0 || any(isnothing, psym_perm)
            missing_syms = NetworkDynamics.sym(m_big)[findall(isnothing, u_big_raw)]
            printstyled("\n  Skipped: $n_missing state(s) of larger model could not be seeded: $(inline_repr(missing_syms))",
                        color=:yellow, bold=true)
        else
            u_big  = Float64[x for x in u_big_raw]
            p_big  = p_small[psym_perm]
            outs_big = copy.(outs_small)
            du_big   = zeros(NetworkDynamics.dim(m_big))

            NetworkDynamics.compfg(m_big)(outs_big, du_big, u_big, ins, p_big, t)

            _mm_diag(m) = begin
                mm = m.mass_matrix
                mm isa UniformScaling ? fill(Int(mm.λ), NetworkDynamics.dim(m)) : LinearAlgebra.diag(mm)
            end
            sym_big   = NetworkDynamics.sym(m_big)
            sym_small = NetworkDynamics.sym(m_small)
            diff_idx_big = findall(!=(0), _mm_diag(m_big))
            alg_idx_big  = findall(==(0), _mm_diag(m_big))
            diff_idx_s   = findall(!=(0), _mm_diag(m_small))
            alg_idx_s    = findall(==(0), _mm_diag(m_small))

            # named diff-state pairs: same symbol name in both models
            named_diff_pairs = Tuple{Int,Int}[]
            matched_big_diff = Set{Int}()
            unmatched_small_diff = Int[]
            for j in diff_idx_s
                idx = findfirst(==(sym_small[j]), sym_big)
                if !isnothing(idx) && idx ∈ diff_idx_big && idx ∉ matched_big_diff
                    push!(named_diff_pairs, (idx, j))
                    push!(matched_big_diff, idx)
                else
                    push!(unmatched_small_diff, j)
                end
            end

            # greedy matching for remaining alg + unmatched diff states
            remaining_big   = [filter(i -> i ∉ matched_big_diff, diff_idx_big); alg_idx_big]
            remaining_small = [unmatched_small_diff; alg_idx_s]
            greedy_matched = 0
            greedy_total   = length(remaining_small)
            if !isempty(remaining_small) && !isempty(remaining_big)
                cost = [min(abs(du_big[i] - du_small[j]), abs(du_big[i] + du_small[j]))
                        for i in remaining_big, j in remaining_small]
                used_big   = falses(length(remaining_big))
                used_small = falses(greedy_total)
                for _ in 1:greedy_total
                    best = Inf; bi = bj = 0
                    for i in eachindex(remaining_big), j in 1:greedy_total
                        !used_big[i] && !used_small[j] && cost[i,j] < best &&
                            (best = cost[i,j]; bi = i; bj = j)
                    end
                    bi == 0 && break
                    used_big[bi] = used_small[bj] = true
                    best < 1e-6 && (greedy_matched += 1)
                end
            end

            maxodiff = maximum(abs.(reduce(vcat, outs_big) .- reduce(vcat, outs_small)), init=0.0)
            maxdiff_diff = isempty(named_diff_pairs) ? 0.0 :
                maximum(abs(du_big[i] - du_small[j]) for (i, j) in named_diff_pairs)

            alg_syms_big   = Set(sym_big[alg_idx_big])
            alg_syms_small = Set(sym_small[alg_idx_s])
            alg_sets_differ = alg_syms_big != alg_syms_small
            if alg_sets_differ
                printstyled("\n  (!) Alg state sets differ — output & remaining-state comparison",
                            color=:yellow)
                printstyled("\n      may be unreliable (models not on a shared constraint manifold)",
                            color=:yellow)
            end
            if maxodiff < 1e-6
                printstyled("\n  Outputs match!         ", color=:green, bold=true)
            else
                c = alg_sets_differ ? :yellow : :red
                printstyled("\n  Outputs dont match!    ", color=c, bold=true)
            end
            print("max diff: ", NetworkDynamics.str_significant(maxodiff; sigdigits=3))
            if isempty(named_diff_pairs)
                printstyled("\n  No named diff states   ", color=:yellow, bold=true)
            elseif maxdiff_diff < 1e-6
                printstyled("\n  Diff states match!     ", color=:green, bold=true)
            else
                printstyled("\n  Diff states dont match!", color=:red, bold=true)
            end
            !isempty(named_diff_pairs) &&
                print("max diff: ", NetworkDynamics.str_significant(maxdiff_diff; sigdigits=3))
            if greedy_total == 0
                # nothing to report
            elseif greedy_matched == greedy_total
                printstyled("\n  Remaining states match! ($greedy_matched/$greedy_total)", color=:green, bold=true)
            else
                c = alg_sets_differ ? :yellow : (greedy_matched == 0 ? :red : :yellow)
                printstyled("\n  Remaining states: $greedy_matched/$greedy_total matched", color=c, bold=true)
            end
        end
    end
    println()

    # return based on default
    return_mtkcompile = MTKCOMPILE_DEFAULT[] === :compare ? false : MTKCOMPILE_DEFAULT[]
    if return_mtkcompile
        return m1 # return the one with mtkcompile
    else
        return m2 # return the one without mtkcompile
    end
end
