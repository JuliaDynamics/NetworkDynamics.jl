"""
    reduce_linear_algebraic

This function was implemented to cover the shortcomings of MTKBase vs full MTK.
It implements a subset of the simplification pipeline previously handled by MTK, which is now
split of into the GPL licensed MTK.

This method is not as powerful but may recover some of whats lost by MTKBase.
It shouldn't to a lot if the system was allready full reduced using MTK.
"""
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

    return sccs
end
