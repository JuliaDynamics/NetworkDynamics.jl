"""
    remove_aliases(eqs, obseqs, states, outputs; verbose)

Detects aliasgroups in eqs and obseqs, finds a "main" state to represent each
aliasgroup and uses that everywhere.
"""
function remove_aliases(eqs, obseqs, states, outputs; verbose)
    ST = SymbolicUtils.BasicSymbolicImpl.var"typeof(BasicSymbolicImpl)"{SymbolicUtils.SymReal}
    aliasgroups = Set{Set{ST}}()

    _collect_and_remove_aliases!(aliasgroups, eqs, states)
    @assert length(states) == length(eqs) "Expected one state per equation after removing aliases, but got $(length(states)) states and $(length(eqs)) equations"
    _collect_and_remove_aliases!(aliasgroups, obseqs)

    alias_subs = Dict{ST, ST}()
    alias_obs = Equation[]

    outset = Set(outputs)
    sortf = function(s)
        if s ∈ outset
            return typemin(Int)
        else
            st = string(s)
            count('₊', st)
        end
    end
    if verbose
        infostr = "Found $(length(aliasgroups)) alias groups:"
    end
    for group in aliasgroups
        candidates = collect(group)
        sort!(candidates, by=sortf)
        main = first(candidates)
        rest = view(candidates, 2:length(candidates))

        if verbose
            infostr *= "\n  $main => " * join(collect(rest), ", ")
        end

        for r in rest
            alias_subs[r] = main
            push!(alias_obs, r ~ main)
        end
    end

    verbose && @info infostr

    _eqs = Symbolics.substitute(eqs, alias_subs)
    _obseqs = Symbolics.substitute(obseqs, alias_subs)
    _insert_sorted!(_obseqs, alias_obs)
    _states = Symbolics.substitute(states, alias_subs)

    _eqs, _obseqs, _states
end
function _collect_and_remove_aliases!(aliasgroups, eqs::Vector{Equation}, states=nothing)
    rm_eq = Int[]
    if !isnothing(states)
        algebraic = [eq_type(eq)[1] == :implicit_algebraic for eq in eqs]
        rm_st = Int[]
        claimed = Set()
    end
    for (i, eq) in enumerate(eqs)
        alias = get_alias(eq)
        alias isa Tuple || continue
        insert_alias!(aliasgroups, alias)
        push!(rm_eq, i)
        if !isnothing(states)
            a, b = alias
            idx = findfirst(j -> algebraic[j] &&
                                 (isequal(states[j], a) || isequal(states[j], b)) &&
                                 states[j] ∉ claimed,
                            eachindex(states))
            if !isnothing(idx)
                push!(rm_st, idx)
                push!(claimed, states[idx])
            end
        end
    end
    deleteat!(eqs, rm_eq)
    if !isnothing(states)
        deleteat!(states, sort!(rm_st))
    end
end
function insert_alias!(aliasgroups, tup)
    a, b = tup
    for group in aliasgroups
        if a ∈ group || b ∈ group
            push!(group, a)
            push!(group, b)
            return
        end
    end
    push!(aliasgroups, Set([a, b]))
end

"""
    reduce_linear_algebraic

This function was implemented to cover the shortcomings of MTKBase vs full MTK.
It implements a subset of the simplification pipeline previously handled by MTK, which is now
split of into the GPL licensed MTK.

This method is not as powerful but may recover some of whats lost by MTKBase.
It shouldn't to a lot if the system was allready full reduced using MTK.
"""
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
    Symbolics.simplify.(sol)
end

function reduce_linear_algebraic(eqs, obseqs, states; outputs=[], ff_inputs=Set(), verbose)
    try
        outset = Set(outputs)
        alg_idx = findall(eqs) do eq
            type = eq_type(eq)
            type[1] == :explicit_algebraic && error("Can't have explicit algebraic at that stage")
            type[1] == :implicit_algebraic
        end
        # Sort so non-output states come first: the bipartite matching tries columns in
        # order, so this naturally prefers solving for non-outputs (internal states).
        sort!(alg_idx, by = i -> states[i] ∈ outset ? 1 : 0)

        if verbose
            @info "Trying to reduce $(length(alg_idx)) implicit algebraic equations" eqs[alg_idx]
        end

        # Expand obs that depend on algebraic states (so coeff matrix sees true state
        # dependencies) and obs that depend on ff_inputs (so FF blocking sees indirect
        # input dependencies through observed variables).
        _transparent_deps = ff_inputs ∪ states[alg_idx]
        obs_subs = _obs_subs_for_states(obseqs, _transparent_deps)
        alg_eqs_expanded = isempty(obs_subs) ? eqs[alg_idx] : Symbolics.substitute(eqs[alg_idx], obs_subs)

        coeff = _build_coeff_mat(alg_eqs_expanded, states[alg_idx]; outset, ff_inputs)

        matches = _match_equations_to_states(coeff)

        if verbose
            @info "Found $(length(matches)) matches in those equations" state_eqs_matching = Any[states[alg_idx[m[2]]] => eqs[alg_idx[m[1]]] for m in matches]
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
            scc_eqs = alg_eqs_expanded[scc_rows]
            scc_sts = states[alg_idx[scc_cols]]

            try
                sol = _symbolic_linear_solve_clean(scc_eqs, scc_sts)
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
        neweqs    = eqs[setdiff(eachindex(eqs), alg_idx[solved_eq_local])]
        newstates = states[setdiff(eachindex(states), alg_idx[solved_st_local])]
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

"""
Build a substitution dict from obseqs for all observed variables that transitively
depend on any of the given states. This is used to expand algebraic equations so
that all state dependencies are explicit (not hidden behind observed variables).
"""
function _obs_subs_for_states(obseqs, states::Set)
    # obs are in topological order: dependencies before dependents.
    # We iterate forward, accumulating which obs symbols depend on states,
    # storing fully-expanded RHS values so a single substitute call suffices.
    subs = Dict()  # lhs => fully-expanded rhs (in terms of states only)
    for eq in obseqs
        vars = get_variables_fix(eq.rhs)
        if any(v -> v ∈ states || haskey(subs, v), vars)
            expanded_rhs = Symbolics.substitute(eq.rhs, subs)
            subs[eq.lhs] = expanded_rhs
        end
    end
    subs
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
