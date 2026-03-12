function expand_alg_equations(eqs, obseqs, states, inputs; verbose)
    alg_idx = findall(eqs) do eq
        type = eq_type(eq)
        type[1] == :explicit_algebraic && error("Can't have explicit algebraic at that stage")
        type[1] == :implicit_algebraic
    end
    isempty(alg_idx) && return eqs

    obs_subs = OrderedDict(eq.lhs => eq.rhs for eq in obseqs)
    if !isempty(obs_subs)
        for i in alg_idx
            eqs[i] = fixpoint_sub(eqs[i], obs_subs)
        end
    end
    eqs
end

function remove_aliases(eqs, obseqs, states, outputs; verbose)
    # 1. Find the alias pairs
    states_set = Set(states)
    alias_pairs = Tuple{ST,ST}[]
    eq_alias_idxs = Int[]
    for (i, eq) in enumerate(eqs)
        alias = get_alias(eq)
        isnothing(alias) && continue
        push!(alias_pairs, alias)
        push!(eq_alias_idxs, i)
    end
    obs_alias_idx = Int[]
    for (i, eq) in enumerate(obseqs)
        alias = get_alias(eq)
        isnothing(alias) && continue
        push!(alias_pairs, alias)
        push!(obs_alias_idx, i)
    end

    # 2. Group them into alias groups
    groups = _alias_connected_components(alias_pairs)

    if verbose
        str = "Found $(length(groups)) alias groups"
        str *= length(groups) > 0 ? ":" : "."
    end

    # 3. Pick main per group, build substitution dict and alias observations
    outset = Set(outputs)
    sortf(s) = s ∈ outset ? typemin(Int) : count('₊', string(s))
    alias_subs = Dict()
    new_alias_obs = Equation[]
    removed_states = Set{ST}()
    for group in groups
        sorted = sort!(collect(group), by=sortf)
        main = first(sorted)
        for r in @view(sorted[2:end])
            alias_subs[r] = main
            push!(new_alias_obs, r ~ main)
            push!(removed_states, r)
        end
        if verbose
            str *= "\n  Group: $main replaces $(inline_repr(sorted[2:end]))"
        end
    end
    verbose && @info str

    # 4. insert obs, reduce states and reduce 
    eqs_new = let
        eqs_new = copy(eqs)
        eqs_new = deleteat!(eqs_new, eq_alias_idxs)
        eqs_new = Symbolics.substitute_in_deriv.(eqs_new, Ref(alias_subs))
    end

    obseqs_new = let
        obseqs_new = copy(obseqs)
        obseqs_new = deleteat!(obseqs_new, obs_alias_idx)
        obseqs_new = Symbolics.substitute.(obseqs_new, Ref(alias_subs))
        _insert_sorted!(obseqs_new, new_alias_obs)
        obseqs_new
    end

    # reducing the state vector is tricky!
    # essentially we need a primitive matching to match diff states where they belong
    # and then fill the rest with non-diff states
    unmatched_states = setdiff(states, removed_states)
    states_new = match_diff_states(eqs_new, unmatched_states)

    eqs_new, obseqs_new, states_new
end
"""
    _alias_connected_components(pairs)

Given alias pairs `[(a,b), (c,d), ...]`, returns groups of transitively connected
variables as a `Vector{Vector}`. E.g. pairs `(a,b), (b,c)` → `[[a, b, c]]`.
"""
function _alias_connected_components(pairs)
    adj = Dict{ST, Vector{ST}}()
    for (a, b) in pairs
        push!(get!(adj, a, ST[]), b)
        push!(get!(adj, b, ST[]), a)
    end
    visited = Set{ST}()
    groups = Vector{Vector{ST}}()
    for v in keys(adj)
        v ∈ visited && continue
        group = ST[]
        stack = ST[v]
        while !isempty(stack)
            u = pop!(stack)
            u ∈ visited && continue
            push!(visited, u)
            push!(group, u)
            for w in adj[u]
                w ∈ visited || push!(stack, w)
            end
        end
        push!(groups, group)
    end
    return groups
end

"""
    reduce_linear_algebraic

This function was implemented to cover the shortcomings of MTKBase vs full MTK.
It implements a subset of the simplification pipeline previously handled by MTK, which is now
split off into the GPL-licensed MTK.

This method is not as powerful but may recover some of what's lost by MTKBase.
It shouldn't do much if the system was already fully reduced using MTK.
"""
function reduce_linear_algebraic(eqs, obseqs, states; outputs=[], ff_inputs=Set(), verbose)
    try
        outset = Set(outputs)
        alg_idx = findall(eqs) do eq
            type = eq_type(eq)
            type[1] == :implicit_algebraic
        end
        # Sort so non-output states come first: the bipartite matching tries columns in
        # order, so this naturally prefers solving for non-outputs (internal states).
        sort!(alg_idx, by = i -> states[i] ∈ outset ? 1 : 0)

        if verbose
            str = "Trying to reduce $(length(alg_idx)) implicit algebraic equations"
            str *= "\n" * multiline_repr(eqs[alg_idx], prefix="  ")
            @info str
        end

        alg_eqs = eqs[alg_idx]

        coeff = _build_coeff_mat(alg_eqs, states[alg_idx]; outset, ff_inputs)

        matches = _match_equations_to_states(coeff)

        if verbose
            str =  "Found $(length(matches)) matches in those equations"
            states_eqs_matching = OrderedDict(states[alg_idx[m[2]]] => eqs[alg_idx[m[1]]] for m in matches)
            str *= "\n" * multiline_repr(states_eqs_matching, prefix="  ")
            @info str
        end

        isempty(matches) && return eqs, obseqs, states

        # decompose matched pairs into SCCs and solve per-SCC
        sccs = _matched_sccs(coeff, matches)
        verbose && @info "Decomposed matches into groups $(length(sccs)) groups to solve:"

        solved_eq_local = Int[]
        solved_st_local = Int[]
        newobs = Equation[]
        for (gi, scc) in enumerate(sccs)
            scc_rows = [matches[i][1] for i in scc]
            scc_cols = [matches[i][2] for i in scc]
            # use expanded equations so all state dependencies are visible to the solver
            scc_eqs = alg_eqs[scc_rows]
            scc_sts = states[alg_idx[scc_cols]]

            try
                sol = _symbolic_linear_solve_clean(scc_eqs, scc_sts)
                append!(newobs, scc_sts .~ sol)
                append!(solved_eq_local, scc_rows)
                append!(solved_st_local, scc_cols)
                if verbose
                    str = "Solved group $gi for $(inline_repr(scc_sts))"
                    str *= "\n" * multiline_repr(scc_eqs, prefix="  ")
                    str *= "\nSolution:"
                    soldict = OrderedDict(scc_sts .=> sol)
                    str *= "\n" * multiline_repr(soldict, prefix="  ")
                    @info str
                end
            catch
                if verbose
                    str = "Failed to solve group $gi for $(inline_repr(scc_sts))"
                    str *= "\n" * multiline_repr(scc_eqs, prefix="  ")
                    @info str
                end
            end
        end

        isempty(newobs) && return eqs, obseqs, states

        _insert_sorted!(obseqs, newobs)
        neweqs    = eqs[setdiff(eachindex(eqs), alg_idx[solved_eq_local])]
        newstates = states[setdiff(eachindex(states), alg_idx[solved_st_local])]
        if verbose
            str = "Left with equations after reduction:\n"
            str *= multiline_repr(newstates .=> neweqs, prefix="  ")
            @info str
        end
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
# Solve a small linear system using Cramer's rule for n≤2, falling back to
# symbolic_linear_solve for larger systems. Cramer's rule puts det(A) as the
# single common denominator, avoiding spurious intermediate divisions (e.g.
# dividing by one component of a complex current when the natural denominator
# is the magnitude squared).
function _symbolic_linear_solve_clean(scc_eqs, scc_sts)
    n = length(scc_sts)
    if n > 2
        return Symbolics.symbolic_linear_solve(scc_eqs, scc_sts)
    end
    # extract coefficient matrix A and RHS b from  0 ~ expr  (i.e. expr = 0)
    exprs = [eq.rhs - eq.lhs for eq in scc_eqs]
    A = Symbolics.jacobian(exprs, collect(scc_sts))
    zero_subs = Dict(s => 0 for s in scc_sts)
    b = [-Symbolics.substitute(e, zero_subs) for e in exprs]
    # Cramer's rule
    sol = if n == 1
        [b[1] / A[1,1]]
    else  # n == 2
        det = A[1,1]*A[2,2] - A[1,2]*A[2,1]
        [(A[2,2]*b[1] - A[1,2]*b[2]) / det,
         (A[1,1]*b[2] - A[2,1]*b[1]) / det]
    end
    Symbolics.simplify.(sol, expand=true)
end

"""
Dependency type of an equation w.r.t. a state variable:
  :none      — state does not appear in the equation
  :linear    — equation is linear in the state (a*x + rest, a may involve other states)
  :nonlinear — equation is nonlinear in the state (e.g. x^2, sin(x))
              also used to block FF matching: output states in equations that depend on
              ff_inputs are marked :nonlinear to prevent solving them (which would
              create a feedforward path later undone with spurious divisions).
"""
function _build_coeff_mat(lineqs, linstates; outset=Set(), ff_inputs=Set())
    coeff = Matrix{Symbol}(undef, length(lineqs), length(linstates))
    # Precompute per-equation FF-input flag (only depends on eq, not state)
    has_ff_input = isempty(ff_inputs) ? falses(length(lineqs)) :
                   [!isempty(Set(get_variables_fix(eq)) ∩ ff_inputs) for eq in lineqs]
    for (j, s) in enumerate(linstates)
        is_output = s ∈ outset
        ex = Symbolics.LinearExpander(s)

        for (i, eq) in enumerate(lineqs)
            a, b, lin = ex(eq)
            # mark output in eq with ff_input dependency as nonlinear to block FF matching
            if !lin || (is_output && has_ff_input[i])
                coeff[i, j] = :nonlinear
            elseif isequal(unwrap_const(a), 0)
                coeff[i, j] = :none
            else
                coeff[i, j] = :linear
            end
        end
    end
    coeff
end

function _maximum_bipartite_matching(adj::AbstractMatrix{Bool})
    nrow, ncol = size(adj)
    cm = zeros(Int, ncol)

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
        try_augment!(row, falses(ncol))
    end
    return cm  # cm[j] = row matched to col j, 0 if unmatched
end

function _match_equations_to_states(coeff)
    n_st = size(coeff, 2)
    adj = coeff .== :linear
    cm = _maximum_bipartite_matching(adj)
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
            if coeff[r, c] === :linear
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

# Returns (diff_term, coeff) if lhs is of the form coeff * D(x), otherwise nothing.
# Handles: numeric scalar, symbolic scalar, product of scalars (a*b*D(x)), negative coefficients.
function get_scaled_diff(lhs)
    lhs = SymbolicUtils.unwrap(lhs)
    @match lhs begin
        SymbolicUtils.BSImpl.AddMul(; coeff, dict) => begin
            diff_keys = [(k, v) for (k, v) in dict
                         if iscall(k) && operation(k) isa Differential && isequal(v, 1)]
            length(diff_keys) == 1 || return nothing
            diff_term = diff_keys[1][1]
            other_keys = [(k, v) for (k, v) in dict if !isequal(k, diff_term)]
            # guard: no other differential factors (e.g. D(x)*D(y))
            any(kv -> iscall(kv[1]) && operation(kv[1]) isa Differential, other_keys) && return nothing
            coeff_expr = coeff * prod(k^v for (k, v) in other_keys; init=1)
            (diff_term, coeff_expr)
        end
        _ => nothing
    end
end

function get_alias(eq, t=ModelingToolkitBase.t_nounits)
    vars = get_variables_fix(eq)
    length(vars) == 2 || return nothing
    a, b = vars
    # parameters are not aliases — an equation like `x ~ V` (V a parameter) defines x,
    # it does not alias two unknowns
    any(ModelingToolkitBase.isparameter, (a, b)) && return nothing
    match = if isequal(unwrap_const(eq.lhs), 0)
        #zero lhs
        isequal(eq.rhs, a-b) || isequal(eq.rhs, b-a)
    else # nonzero rhs
        isequal(eq, a ~ b) || isequal(eq, b ~ a)
    end
    if match
        return (a, b)
    else
        return nothing
    end
end
