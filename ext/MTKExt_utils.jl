
"""
    eq_type(eq::Equation)

Checks the type of the equation. Returns:
- `(:explicit_diffeq, lhs_variable)` for explicit differential equations
- `(:explicit_algebraic, lhs_variable)` for explicit algebraic equations
- `(:implicit_algebraic, lhs_variable)` for implicit algebraic equations

"""
function eq_type(eq::Equation)
    rhs_differentials = _collect_differentials(eq.rhs)
    if !isempty(rhs_differentials)
        throw(ArgumentError("Got equation $eq which contains differentials on the rhs, this is currently not supported!"))
    end
    @match eq.lhs begin
        SymbolicUtils.BSImpl.Const(; val) => begin
            if val != 0
                throw(ArgumentError("Got equation $eq which is a non-zero constant on the lhs, this is currently not supported!"))
            end
            (:implicit_algebraic, nothing)
        end
        SymbolicUtils.BSImpl.Term(; f=Differential(t,1), args=symvec) => begin
            @argcheck length(symvec) == 1 "Diff. eq $eq has more than one variable in lhs!"
            (:explicit_diffeq, only(symvec))
        end
        SymbolicUtils.BSImpl.Term(; f=SymbolicUtils.BSImpl.Sym(), args=args) => begin
            # classic sym call x(t)
            @argcheck length(args) == 1 "Can't parse equation $eq, lhs is a Sym call with more than one argument!"
            rhs_vars = get_variables_fix(eq.rhs)
            if eq.lhs ∈ rhs_vars
                (:implicit_algebraic, eq.lhs)
            else
                (:explicit_algebraic, eq.rhs)
            end
        end
        _ => throw(ArgumentError("Can't determine eq type of $eq."))
    end
end

"""
    lhs_var(eq::Equation)

Returns the variable on the lhs of the equation for equations.
"""
lhs_var(eq::Equation) = eq_type(eq)[2]

function rhs_differentials(eqs::Vector{Equation})
    diffs = Set{BasicSymbolic}()
    for eq in eqs
        _collect_differentials!(diffs, eq.rhs)
    end
    return diffs
end
rhs_differentials(eq::Equation) = _collect_differentials!(Set{BasicSymbolic}(), eq.rhs)

_collect_differentials(ex) = _collect_differentials!(Set{BasicSymbolic}(), ex)

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
    vars = get_variables_fix(expr)
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
    defdict = ModelingToolkit.get_defaults(sys)
    for (s, v) in zip(symbols, defs)
        defdict[s] = v
    end
    sys
end

# WORKAOROUND: get_variables does not descend into Differential anymore
function get_variables_fix(ex)
    set = get_variables(ex)
    new = map(collect(set)) do v
        @match v begin
            SymbolicUtils.BSImpl.Term(; f=Differential(t,1), args=symvec) => begin
                only(symvec)
            end
            _ => v
        end
    end
    typeof(set)(new)
end

####
#### linear algebraic equation reduction
####
function reduce_linear_algebraic(eqs, obseqs, states; outputs=[], verbose)
    types = [eq_type(eq)[1] for eq in eqs]
    try
        @assert !any(==(:explicit_algebraic), types) "Can't have explicit algebraic at that stage"

        alg_mask = [t == :implicit_algebraic for t in types]
        any(alg_mask) || return eqs, obseqs, states
        alg_idx = findall(alg_mask)

        if verbose
            @info "Trying to reduce $(length(alg_idx)) implicit algebraic equations" eqs[alg_idx]
        end

        coeff = _build_coeff_mat(eqs[alg_idx], states[alg_idx])
        outset = Set(outputs)
        is_output = [s ∈ outset for s in states[alg_idx]]
        matches = _match_equations_to_states(coeff, states[alg_idx], is_output)

        if verbose
            @info "Found $(length(matches)) matches in those equations" state_eqs_matching = Any[states[alg_idx[m[2]]] => eqs[alg_idx[m[1]]] for m in matches]
        end

        isempty(matches) && return eqs, obseqs, states

        # decompose matched pairs into SCCs and solve per-SCC
        sccs = _matched_sccs(coeff, matches)
        verbose &&@info "Decomposed matches into groups $(length(sccs)) groups to solve:"

        solved_eq_local = Int[]
        solved_st_local = Int[]
        newobs = Equation[]
        for (gi, scc) in enumerate(sccs)
            scc_rows = [matches[i][1] for i in scc]
            scc_cols = [matches[i][2] for i in scc]
            scc_eqs = eqs[alg_idx[scc_rows]]
            scc_sts = states[alg_idx[scc_cols]]

            try
                sol = Symbolics.symbolic_linear_solve(scc_eqs, scc_sts)
                append!(newobs, scc_sts .~ sol)
                append!(solved_eq_local, scc_rows)
                append!(solved_st_local, scc_cols)
                verbose && @info "Solved group $gi: $scc" scc_eqs scc_sts .~ sol
            catch
                # SCC not jointly linear (nonlinear cross-terms) → leave as constraints
                verbose && @info "Failed to solve group $gi for $scc_sts" scc_eqs
            end
        end

        isempty(newobs) && return eqs, obseqs, states

        _insert_sorted!(obseqs, newobs)
        keep_eq_idx = sort(union(findall(.!alg_mask), alg_idx[setdiff(1:length(alg_idx), solved_eq_local)]))
        keep_st_idx = sort(union(findall(.!alg_mask), alg_idx[setdiff(1:length(alg_idx), solved_st_local)]))
        neweqs = eqs[keep_eq_idx]
        newstates = states[keep_st_idx]
        verbose && @info "Left with equations after reduction:" newstates .=> neweqs
        return neweqs, obseqs, newstates
    catch e
        @warn "Reduction failed, fall back" exception=(e, catch_backtrace())
        eqs, obseqs, states
    end
end
function _insert_sorted!(obseqs, newobs)
    # newobs must be in topological order (dependencies before dependents).
    # Insert sequentially so each equation finds its dependencies already in obseqs.
    # We insert each equation just after the last equation in obseqs whose LHS appears
    # on the RHS of eq (i.e., after all its obs dependencies are defined).
    # Note: existing obs may reference `eq.lhs` as a state variable (not as an obs),
    # so we do not check first_dependant here — that check is unsound when states
    # and obs share the same symbol.
    for eq in newobs
        obssym = [e.lhs for e in obseqs]
        vars = get_variables_fix(eq.rhs)
        idx = findlast(sym -> sym ∈ vars, obssym)
        last_dependency = isnothing(idx) ? 0 : idx
        insert!(obseqs, last_dependency+1, eq)
    end
    nothing
end
function _build_coeff_mat(lineqs, linstates)
    expanders = [Symbolics.LinearExpander(s) for s in linstates]
    coeff = Matrix{Any}(nothing, length(lineqs), length(linstates))
    for (j, ex) in enumerate(expanders)
        for (i, eq) in enumerate(lineqs)
            a, b, lin = ex(eq)
            if lin
                coeff[i, j] = a
            end
        end
    end
    coeff
end

function _is_zero_coeff(c)
    # isequal(unwrap_const(c), 0)
    c isa Number && return iszero(c)
    if c isa BasicSymbolic
        uc = unwrap_const(c)
        return uc isa Number && iszero(uc)
    end
    return false
end

function _is_constant_coeff(c, state_vars)
    c isa Number && return true
    cvars = get_variables(c)
    return isempty(intersect(cvars, state_vars))
end

function _maximum_bipartite_matching(adj::AbstractMatrix{Bool}; col_match::Union{Nothing,Vector{Int}}=nothing)
    nrow, ncol = size(adj)
    cm = isnothing(col_match) ? zeros(Int, ncol) : copy(col_match)
    matched_rows = Set(c for c in cm if c > 0)

    function try_augment!(row, visited)
        for col in 1:ncol
            adj[row, col] && !visited[col] || continue
            visited[col] = true
            if cm[col] == 0 || try_augment!(cm[col], visited)
                cm[col] = row
                return true
            end
        end
        return false
    end

    for row in 1:nrow
        row in matched_rows && continue
        try_augment!(row, falses(ncol))
    end
    return cm  # cm[j] = row matched to col j, 0 if unmatched
end

function _match_equations_to_states(coeff, state_vars, is_output)
    n_eq, n_st = size(coeff)
    adj_const = [!isnothing(coeff[i,j]) && !_is_zero_coeff(coeff[i,j]) && _is_constant_coeff(coeff[i,j], state_vars) for i in 1:n_eq, j in 1:n_st]
    adj_full  = [!isnothing(coeff[i,j]) && !_is_zero_coeff(coeff[i,j]) for i in 1:n_eq, j in 1:n_st]

    # mask out output columns, match non-outputs first, then extend to outputs
    not_out = .!is_output
    cm = _maximum_bipartite_matching(adj_const .* not_out')
    cm = _maximum_bipartite_matching(adj_full .* not_out'; col_match=cm)
    # now extend matching to include output columns
    cm = _maximum_bipartite_matching(adj_const; col_match=cm)
    cm = _maximum_bipartite_matching(adj_full; col_match=cm)

    return [(cm[j], j) for j in 1:n_st if cm[j] > 0]
end

"""
Decompose matched pairs into SCCs of the dependency digraph.
Returns SCCs in topological order (dependencies before dependents).
Each SCC is a vector of indices into `matches`.
"""
function _matched_sccs(coeff, matches)
    n = length(matches)
    n == 0 && return Vector{Int}[]

    # map matched col → index in matches
    col_to_idx = Dict(c => i for (i, (_, c)) in enumerate(matches))

    # build adjacency list: match i depends on match j if equation of i
    # has a nonzero coefficient for the state of j
    adj = [Int[] for _ in 1:n]
    for (i, (r, _)) in enumerate(matches)
        for (j, (_, c)) in enumerate(matches)
            i == j && continue
            if !isnothing(coeff[r, c]) && !_is_zero_coeff(coeff[r, c])
                push!(adj[i], j)
            end
        end
    end

    # Tarjan's SCC algorithm
    index_counter = Ref(0)
    stack = Int[]
    on_stack = falses(n)
    indices = fill(-1, n)
    lowlinks = fill(-1, n)
    sccs = Vector{Int}[]

    function strongconnect(v)
        indices[v] = index_counter[]
        lowlinks[v] = index_counter[]
        index_counter[] += 1
        push!(stack, v)
        on_stack[v] = true

        for w in adj[v]
            if indices[w] == -1
                strongconnect(w)
                lowlinks[v] = min(lowlinks[v], lowlinks[w])
            elseif on_stack[w]
                lowlinks[v] = min(lowlinks[v], indices[w])
            end
        end

        if lowlinks[v] == indices[v]
            scc = Int[]
            while true
                w = pop!(stack)
                on_stack[w] = false
                push!(scc, w)
                w == v && break
            end
            push!(sccs, scc)
        end
    end

    for v in 1:n
        indices[v] == -1 && strongconnect(v)
    end

    # Tarjan's yields SCCs in topological order (dependencies before dependents); no reversal needed
    return sccs
end
