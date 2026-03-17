
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

"""
    lhs_var(eq::Equation)

Returns the variable on the lhs of the equation for equations.
"""
lhs_var(eq::Equation) = eq_type(eq)[2]

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
    r = getproperty(sys, Symbol(parts[1]); namespace=false)
    for part in parts[2:end]
        r = getproperty(r, Symbol(part); namespace=true)
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
            # TODO: upstream fixed, can be removed
            # https://github.com/SciML/ModelingToolkit.jl/issues/3016
            # getproperty(sys, getname(invalidname); namespace=false)
            getproperty_symbolic(sys, invalids) # like getproperty but works on namespaced symbols foo₊bar directly
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
    _unsortedstates = Set(states) # creates copy
    for (i, eq) in pairs(eqs)
        type, var = eq_type(eq)
        if type == :explicit_diffeq
            states_new[i] = pop!(_unsortedstates, var)
        elseif type == :implicit_algebraic
            # do nothing
        else
            error("Got equation $eq of type $type, which is not expected at that stage (only explicit diff eqs and implicit alg. eqs). Something went wrong.")
        end
    end
    for i in eachindex(states_new)
        if !isassigned(states_new, i)
            states_new[i] = pop!(_unsortedstates)
        end
    end
    @assert isempty(_unsortedstates) "Not all states could be matched, something went wrong"
    @assert length(eqs) == length(states_new)
    states_new
end

# WORKAOROUND: get_variables does not descend into Differential anymore
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

"""
    guesses_to_guessformulas!(sys; alias_substitutions=Dict())

Extracts symbolic guesses from a (simplified) system and returns matching GuessFormulas.
The symbolic guesses are removed from the system's guesses dict in-place (hence `!`),
so they won't be silently dropped by metadata extraction.
Constant (numeric) guesses are left in place — they are handled as `:guess` metadata.
"""
function guesses_to_guessformulas!(sys; obs_subs=Dict())
    allguesses = ModelingToolkitBase.get_guesses(sys)
    isempty(allguesses) && return nothing

    guessformulas = Set{GuessFormula}()

    for (lhs_sym, rhs_expr) in collect(allguesses)
        # skip constant guesses — leave them for :guess metadata
        isempty(get_variables(rhs_expr)) && continue

        # remove from system guesses dict
        delete!(allguesses.dict, lhs_sym)

        lhs_sub  = Symbolics.substitute(lhs_sym,  obs_subs)
        rhs_sub  = fixpoint_sub(rhs_expr, obs_subs)

        target = [getname(lhs_sub)]
        input_symbolic = collect(get_variables(rhs_sub))
        input_names = getname.(input_symbolic)
        f = Symbolics.build_function([rhs_sub], input_symbolic; expression=Val(false))[2]

        rhsstring = repr(rhs_sub)
        for input in input_symbolic
            rhsstring = replace(rhsstring, repr(input) => "u["*repr(getname(input))*"]")
        end
        prettyprint = """
        GuessFormula([$(join(repr.(target), ", "))], [$(join(repr.(input_names), ", "))]) do out, u
            out[$(repr(only(target)))] = $(rhsstring)
        end"""
        push!(guessformulas, GuessFormula(f, target, input_names, prettyprint))
    end

    _warn_duplicate_formula_targets(guessformulas, "GuessFormula")
    isempty(guessformulas) ? nothing : guessformulas
end

"""
    bindings_to_initformulas(sys)

Extractes the `bindings` from a system and returns matchin InitFormulas.
"""
function bindings_to_initformulas(sys; obs_subs=Dict())
    bindings = ModelingToolkitBase.bindings(sys)
    isempty(bindings) && return nothing

    initformulas = Set{InitFormula}()

    for (lhs, rhs) in bindings
        lhs = Symbolics.substitute(lhs, obs_subs)
        rhs = fixpoint_sub(rhs, obs_subs)
        target = [getname(lhs)]
        input_symbolic = collect(get_variables(rhs))
        input_names = getname.(input_symbolic)
        f = Symbolics.build_function([rhs], input_symbolic; expression=Val(false))[2]

        rhsstring = repr(rhs)
        for input in input_symbolic
            rhsstring = replace(rhsstring, repr(input) => "u["*repr(getname(input))*"]")
        end
        prettyprint  = """
        InitFormula([$(join(repr.(target), ", "))], [$(join(repr.(input_names), ", "))]) do out, u
            out[$(repr(only(target)))] = $(rhsstring)
        end"""
        push!(initformulas, InitFormula(f, target, input_names, prettyprint))
    end
    _warn_duplicate_formula_targets(initformulas, "InitFormula")
    return initformulas
end

function _warn_duplicate_formula_targets(formulas, kind)
    seen = Set{Symbol}()
    duplicates = Symbol[]
    for f in formulas, t in f.outsym
        t ∈ seen ? push!(duplicates, t) : push!(seen, t)
    end
    if !isempty(unique!(duplicates))
        @warn "Multiple $kind formulas target the same symbol(s) $duplicates after obs substitution. \
               This is not supported yet and may lead to unclear behavior."
    end
end

function NetworkDynamics.multiline_repr(eqs::Vector{Equation}; prefix="")
    lines = map(eqs) do eq
        prefix * repr(eq.lhs) * " &~ " * repr(eq.rhs)
    end
    join(NetworkDynamics.align_strings(lines), "\n")
end
