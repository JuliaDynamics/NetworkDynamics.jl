
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
function NetworkDynamics.set_mtk_defaults!(sys::System, pairs)
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
    defdict = ModelingToolkitBase.get_initial_conditions(sys)
    for (s, v) in zip(symbols, defs)
        defdict[s] = v
    end
    sys
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
# references other variables/parameters is NOT constant — for guesses/bindings those belong to
# the GuessFormula/InitFormula path, while constants are stored as plain `:guess` metadata.
# Shared by `_get_metadata` (keep numeric guesses) and `guesses_to_guessformulas!` (skip them).
is_symbolic_constant(x) = (u = unwrap(x); u isa Number || isempty(get_variables(u)))

# Numeric value of a symbolic constant; only meaningful when `is_symbolic_constant(x)`.
symbolic_constant_value(x) = (u = unwrap(x); u isa Number ? u : Float64(Symbolics.value(u)))

"""
    guesses_to_guessformulas!(sys; obs_subs=Dict())

Extracts symbolic guesses from a (simplified) system and returns matching GuessFormulas.
The symbolic guesses are removed from the system's guesses dict in-place (hence `!`),
so they won't be silently dropped by metadata extraction.
Constant (numeric) guesses are left in place — they are handled as `:guess` metadata.
"""
function guesses_to_guessformulas!(sys; obs_subs=Dict())
    # `get_guesses` (raw backing dict), not `guesses` (a freshly-constructed, recursively
    # namespaced copy): `sys` is already the flattened/simplified system here, so the two are
    # identical — but we `delete!` symbolic guesses below and need the *mutable* backing dict
    # for that to stick. (`_get_metadata` uses `guesses(sys)` because it runs on the original,
    # possibly hierarchical, system where only the recursive accessor namespaces correctly.)
    @assert isempty(ModelingToolkitBase.get_systems(sys)) "Should be called on flattend/simplified system!"
    allguesses = ModelingToolkitBase.get_guesses(sys)
    isempty(allguesses) && return nothing

    # A guess may target any guessable symbol (unknown or parameter — params can be free
    # during init), so targets and inputs share the same valid set here.
    valid_inputs = valid_targets = Set(getname.(vcat(unknowns(sys), parameters(sys))))

    resolved = Any[]
    for (lhs_sym, rhs_expr) in collect(allguesses)
        # skip constant guesses — leave them for :guess metadata
        is_symbolic_constant(rhs_expr) && continue
        # This is a symbolic guess: it can never be a numeric :guess metadata entry, so
        # remove it from the system dict unconditionally (whether or not it becomes a
        # GuessFormula below) to avoid it being mishandled later.
        delete!(allguesses.dict, lhs_sym)

        r = _resolve_formula(lhs_sym, rhs_expr; obs_subs, valid_inputs, valid_targets, kind="Guess")
        r === nothing && continue
        push!(resolved, r)
    end

    # A GuessFormula is only a convergence hint, so conflicting targets (alias-merged
    # states carrying different guess expressions) are deduped with a warning, never fatal.
    resolved = _dedupe_resolved(resolved; fail=:warn, kind="GuessFormula")
    guessformulas = Set(_build_formula(GuessFormula, r) for r in resolved)
    isempty(guessformulas) ? nothing : guessformulas
end

"""
    bindings_to_initformulas(sys; states, obs_subs)

Extracts the `bindings` from a system and returns matching InitFormulas.

This will **ignore** parameter bindings! Parameter bindings will become observed
equations in an earlier step.

An InitFormula is a constraint (not a hint), so conflicting targets (two bindings forcing
the same state to *different* values after alias substitution) are a genuine
over-determination and raise an error.
"""
function bindings_to_initformulas(sys; states::Vector{ST}, obs_subs)
    bindings = ModelingToolkitBase.bindings(sys)
    isempty(bindings) && return nothing

    valid_inputs  = Set(getname.(vcat(unknowns(sys), parameters(sys))))
    valid_targets = Set(getname.(states))

    resolved = Any[]
    for (_lhs, _rhs) in bindings
        ismissing(unwrap_const(_rhs)) && continue # FIXME somehow, mtkcompile tends to spit out missing bindings?
        _lhs ∈ ModelingToolkitBase.bound_parameters(sys) && continue # skip parameter bindings, they will become observed equations
        r = _resolve_formula(_lhs, _rhs; obs_subs, valid_inputs, valid_targets, kind="Binding")
        r === nothing && continue
        push!(resolved, r)
    end

    resolved = _dedupe_resolved(resolved; fail=:error, kind="InitFormula")
    isempty(resolved) ? nothing : Set(_build_formula(InitFormula, r) for r in resolved)
end

# ── shared helpers for guesses_to_guessformulas! and bindings_to_initformulas ──
#
# Both an InitFormula and a GuessFormula are "set `target := f(inputs)`" objects, so the
# resolution, validity guards, conflict-deduplication and function-building are identical;
# only the upstream filtering (numeric guess vs. parameter binding) and the conflict policy
# (`fail`) live in the two callers above.

# Resolve one `lhs => rhs` pair against the observed-substitution map into a canonical
# `(target, rhs, inputs)` form, or skip it (with a warning) when it cannot become a
# well-formed formula. `obs_subs` collapses every alias-group member onto its surviving
# "main" symbol, so two entries on different members of one alias group resolve to the *same*
# target here — exactly the collision `_dedupe_resolved` then handles.
function _resolve_formula(lhs_sym, rhs_expr; obs_subs, valid_inputs, valid_targets, kind)
    # (1) target must resolve to a single surviving unknown; an eliminated variable expands
    #     to a compound expression after alias substitution (`getname` would then throw).
    lhs_sub  = Symbolics.substitute(lhs_sym, obs_subs)
    lhs_vars = get_variables(lhs_sub)
    if length(lhs_vars) != 1 || !isequal(only(lhs_vars), lhs_sub)
        @warn "$kind for $lhs_sym expands to $lhs_sub after alias substitution, which is not \
               a single unknown (likely eliminated). Skip."
        return nothing
    end
    # (2) it must not resolve to an observable (those are computed, not set).
    if lhs_sub ∈ keys(obs_subs)
        @warn "$kind for $lhs_sym resolves to observable $lhs_sub, which is computed rather \
               than set. Skip."
        return nothing
    end
    target = getname(lhs_sub)
    # (3) and it must be a settable target of the system (a surviving unknown/state).
    if target ∉ valid_targets
        @warn "$kind for $lhs_sym resolves to $target, which is not a settable target of the \
               system. Skip."
        return nothing
    end

    rhs_sub        = fixpoint_sub(rhs_expr, obs_subs)
    input_symbolic = collect(get_variables(rhs_sub))
    input_names    = getname.(input_symbolic)
    # (4) no self-dependency (the formula constructors forbid it).
    if target ∈ input_names
        @warn "$kind for $lhs_sym depends on its own resolved target $target after \
               substitution. Skip."
        return nothing
    end
    # (5) every input must be a symbol the system actually exposes (a surviving unknown or
    #     parameter); a rhs referencing anything else is a malformed formula that would be
    #     rejected fatally downstream (`assert_*formula_compat`), so skip it with a warning.
    #     This is a generic well-formedness guard, NOT a namespace workaround: guesses are
    #     always read through `guesses`/`get_guesses` on the enclosing (flattened) system,
    #     which namespace correctly, so a well-formed guess never lands here.
    badins = setdiff(Set(input_names), valid_inputs)
    if !isempty(badins)
        @warn "$kind for $lhs_sym references symbol(s) not in the system \
               ($(join(badins, ", "))). Skip."
        return nothing
    end

    (; src=lhs_sym, target, rhs=rhs_sub, input_symbolic, input_names)
end

# Collapse entries that resolved to the same target. Identical definitions (equal rhs after
# obs substitution, e.g. alias-merged states carrying the same expression) are always safe to
# dedupe silently. Genuinely *conflicting* definitions (same target, differing rhs) are a
# warning for guesses (`fail=:warn`, only a hint) but an error for bindings (`fail=:error`,
# two constraints forcing one state). Comparison is on the symbolic rhs, not the prettyprint.
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
            @debug "$kind: $(length(group)) identical definitions for $target after obs \
                    substitution; keeping one."
        else
            msg = "$kind: conflicting definitions target $target after obs substitution \
                   (differing right-hand sides). Keeping $(repr(chosen.rhs)); dropping the rest."
            fail === :error ? error(msg) : @warn msg
        end
    end
    kept
end

# Build the actual Init/GuessFormula from a resolved entry. Identical for both formula types
# apart from the type name baked into the prettyprint block.
function _build_formula(::Type{FT}, r) where {FT}
    label = string(nameof(FT))
    f = Symbolics.build_function([r.rhs], r.input_symbolic; expression=Val(false))[2]

    rhsstring = repr(r.rhs)
    for input in r.input_symbolic
        rhsstring = replace(rhsstring, repr(input) => "u[" * repr(getname(input)) * "]")
    end
    prettyprint = """
    $label([$(repr(r.target))], [$(join(repr.(r.input_names), ", "))]) do out, u
        out[$(repr(r.target))] = $(rhsstring)
    end"""
    FT(f, [r.target], r.input_names, prettyprint)
end

"""
    extract_aliasmap(c::ComponentModel, obseqs::Vector{Equation})

Extracts the `AliasMap` of a compiled component from its observed equations, i.e.
every observable which is a pure `factor * settable_symbol` alias.

Note that `pick_best_alias_names` already consolidates *identity* alias groups onto a single
representative and re-inserts direct `alias ~ main` observations; what remains for us are the
scaled/sign-flipped aliases (`get_alias` matches no coefficient) plus chains through them.

Anything that is not a pure scaled alias — affine (`x + 1`), sums, nonlinear terms, symbolic
coefficients — as well as chains whose root is not settable, stays an ordinary observable.
"""
function extract_aliasmap(c::NetworkDynamics.ComponentModel, obseqs)
    settable = settable_symbols(c)

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
        resolved = _resolve_alias(alias, steps, settable, Symbol[])
        isnothing(resolved) && continue
        # never an identity entry: the root is settable, `alias` is not
        am[alias] = resolved
    end
    am
end

# Match `ex` against `factor * var` with a numeric, nonzero `factor` and a single variable.
# Returns `(factor, varname)` or `nothing`.
#
# `linear_expansion(ex, v)` decomposes `ex` into `(factor, offset, islinear)` such that
# `ex == factor*v + offset`, which is precisely the shape we accept. The guards reject, in
# order: multi-variable rhs (`x + y`; also `k*x`, since `get_variables` counts parameters),
# nonlinear rhs (`x^2`, `sin(x)` — not linear), a symbolic factor or offset (`k*x` again),
# affine rhs (`x + 1` — nonzero offset) and a degenerate `0*x`.
#
# `linear_expansion` hands back MTKv11 symbolic literals, so both parts go through
# `unwrap_const`; requiring the result to be a `Number` keeps this to numeric literals, and a
# variable-free but unfolded coefficient is conservatively left as an ordinary observable.
#
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

# Follow one-step links until the root is settable, multiplying factors along the way
# (`obs2 ~ -obs1`, `obs1 ~ 2x` ⇒ `(-2.0, :x)`). Returns `nothing` when the chain dead-ends on
# a non-settable symbol.
function _resolve_alias(s, steps, settable, visiting)
    haskey(steps, s) || return nothing
    if s ∈ visiting
        error("Cyclic alias chain detected at :$s via $(join(visiting, " → ")). \
               Observed equations are topologically sorted, this should never happen.")
    end
    factor, target = steps[s]
    target ∈ settable && return (factor, target)

    push!(visiting, s)
    rest = _resolve_alias(target, steps, settable, visiting)
    pop!(visiting)
    isnothing(rest) && return nothing
    (factor * rest[1], rest[2])
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
