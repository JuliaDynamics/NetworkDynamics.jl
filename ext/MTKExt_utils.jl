
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

This will **ingore** parameter bindings! Parameter bindings will become observed
equations in an earlier step.
"""
function bindings_to_initformulas(sys; states::Vector{ST}, obs_subs)
    bindings = ModelingToolkitBase.bindings(sys)
    isempty(bindings) && return nothing

    initformulas = Set{InitFormula}()
    for (_lhs, _rhs) in bindings
        _lhs ∈ ModelingToolkitBase.bound_parameters(sys) && continue # skip parameter bindings, they will become observed equations
        lhs = Symbolics.substitute(_lhs, obs_subs) # apply alias transformation
        rhs = fixpoint_sub(_rhs, obs_subs)

        try
            if lhs ∉ Set(states)
                @warn "Binding $_lhs <= $_rhs targets a non-state variable after expansion in known symbols: $lhs <= $rhs. This most likely means that $_lhs was solved and is not an unknown anymore. Skip formula."
                continue
            end
        catch
            @warn "Could not handle binding with $_lhs <= $_rhs expanded to $lhs <= $rhs. Skip!"
            continue
        end

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

    # try find match between sym and psyms
    if (samesym || all(s -> s ∈ m1.obssym, setdiff(m2.sym, m1.sym))) && samepsyms
        sym_perm = [findfirst(==(s), NetworkDynamics.sym(m1)) for s in NetworkDynamics.sym(m2)]
        psym_perm = [findfirst(==(p), NetworkDynamics.psym(m1)) for p in NetworkDynamics.psym(m2)]

        outs1, du1, u1, ins, p1, t = NetworkDynamics.rand_inputs_fg(m1)
        outs2 = copy.(outs1)
        du2 = copy(du1)
        u2 = map(enumerate(sym_perm)) do (i,p)
            if !isnothing(p)
                u1[p]
            else
                sym = NetworkDynamics.sym(m2)[i]
                obsidx = findfirst(==(sym), NetworkDynamics.obssym(m1))
                obsout = zeros(length(m1.obssym))
                m2.obsf(obsout, u1, ins..., p1, t)
                obsout[obsidx]
            end
        end

        p2 = p1[psym_perm]
        NetworkDynamics.compfg(m1)(outs1, du1, u1, ins, p1, t)
        NetworkDynamics.compfg(m2)(outs2, du2, u2, ins, p2, t)
        # split state indices into diff (mass=1) and alg (mass=0)
        _mm_diag2(m) = begin
            mm = m.mass_matrix; s = NetworkDynamics.sym(m)
            mm isa UniformScaling ? fill(Int(mm.λ), length(s)) : LinearAlgebra.diag(mm)
        end
        diff_idx1 = findall(!=(0), _mm_diag2(m1))
        alg_idx1  = findall(==(0), _mm_diag2(m1))
        diff_idx2 = findall(!=(0), _mm_diag2(m2))
        alg_idx2  = findall(==(0), _mm_diag2(m2))

        maxodiff = maximum(abs.(reduce(vcat, outs1) .- reduce(vcat, outs2)), init=0.0)
        # diff states: sym_perm maps m2 positions → m1; compare directly
        maxdiff_diff = maximum(abs.(du1[sym_perm[diff_idx2]] .- du2[diff_idx2]), init=0.0)
        # alg states: order is arbitrary — find best matching by min-cost bijection
        # cost(i,j) = min(|du1_alg[i] - du2_alg[j]|, |du1_alg[i] + du2_alg[j]|)
        du1_alg = du1[alg_idx1]; du2_alg = du2[alg_idx2]
        n_alg = length(du1_alg)
        alg_cost = [min(abs(du1_alg[i] - du2_alg[j]), abs(du1_alg[i] + du2_alg[j]))
                    for i in 1:n_alg, j in 1:n_alg]
        alg_matched = Pair{Int,Int}[]
        used1_alg = falses(n_alg); used2_alg = falses(n_alg)
        for _ in 1:n_alg
            best = Inf; bi = bj = 0
            for i in 1:n_alg, j in 1:n_alg
                !used1_alg[i] && !used2_alg[j] && alg_cost[i,j] < best &&
                    (best = alg_cost[i,j]; bi = i; bj = j)
            end
            push!(alg_matched, bi => bj)
            used1_alg[bi] = used2_alg[bj] = true
        end
        maxdiff_alg = isempty(alg_matched) ? 0.0 :
            maximum(alg_cost[p.first, p.second] for p in alg_matched)

        printstyled("\nComparison of f & g outputs with random inputs: ", color=:blue, bold=true)
        if maxodiff < 1e-6
            printstyled("\n  Outputs match!         ", color=:green, bold=true)
        else
            printstyled("\n  Outputs dont match!    ", color=:red, bold=true)
        end
        print("max diff: ", NetworkDynamics.str_significant(maxodiff; sigdigits=3))
        if maxdiff_diff < 1e-6
            printstyled("\n  Diff states match!     ", color=:green, bold=true)
        else
            printstyled("\n  Diff states dont match!", color=:red, bold=true)
        end
        print("max diff: ", NetworkDynamics.str_significant(maxdiff_diff; sigdigits=3))
        if maxdiff_alg < 1e-6
            printstyled("\n  Alg states match!      ", color=:green, bold=true)
        else
            printstyled("\n  Alg states dont match! ", color=:red, bold=true)
            alg_syms1 = NetworkDynamics.sym(m1)[alg_idx1]
            alg_syms2 = NetworkDynamics.sym(m2)[alg_idx2]
            for (i, j) in alg_matched
                c = alg_cost[i, j]
                c >= 1e-6 && print("\n    $(alg_syms1[i]) <-> $(alg_syms2[j]): diff=",
                    NetworkDynamics.str_significant(c; sigdigits=3))
            end
        end
        print("max diff: ", NetworkDynamics.str_significant(maxdiff_alg; sigdigits=3))
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
