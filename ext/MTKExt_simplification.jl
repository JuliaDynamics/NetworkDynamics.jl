
"""
    pick_best_alias_names(eqs, obseqs, states, outputs, inputs; verbose)

Post-processing pass that consolidates alias chains produced by `reduce_equations`.

After reduction, the obs list may contain pure alias equations of the form `a ~ b` that
connect several names for the same quantity (e.g. `busbar₊u_r ~ busbar₊terminal₊u_r ~
gen₊terminal₊u_r`).  This pass:
1. Collects all such alias equations from `obseqs`.
2. Groups transitively connected aliases into clusters.
3. Picks one *main* representative per cluster — preferring differential states, then
   input variables, then output variables, then variables with fewer `₊` separators
   (less deeply nested).
4. Applies the substitution `nonmain → main` throughout `eqs`, `obseqs`, and `states`,
   and reinserts canonical `nonmain ~ main` alias observations.
"""
function pick_best_alias_names(eqs, obseqs, states, outputs, inputs; verbose)
    # 1. Find the alias pairs
    alias_pairs = Tuple{ST,ST}[]
    alias_idx = Int[]
    for (i, eq) in enumerate(obseqs)
        alias = get_alias(eq)
        isnothing(alias) && continue
        push!(alias_pairs, alias)
        push!(alias_idx, i)
    end

    # 2. Group them into alias groups
    groups = _alias_connected_components(alias_pairs)

    if verbose
        str = "Found $(length(groups)) alias groups"
        str *= length(groups) > 0 ? ":" : "."
    end

    # 3. Pick main per group, build substitution dict and alias observations
    diffstateset = Set{ST}()
    for eq in eqs
        type, var = eq_type(eq)
        type == :explicit_diffeq && push!(diffstateset, var)
    end
    inset = Set{ST}(inputs)
    outset = Set{ST}(outputs)
    sortf = function(s)
        prio = if s ∈ diffstateset
            0
        elseif s ∈ inset
            1
        elseif s ∈ outset
            2
        else
            3
        end
        (prio, count('₊', string(s)))
    end
    alias_subs = Dict()
    new_alias_obs = Equation[]
    for group in groups
        sorted = sort!(collect(group), by=sortf)
        main = first(sorted)
        for r in @view(sorted[2:end])
            alias_subs[r] = main
            push!(new_alias_obs, r ~ main)
        end
        if verbose
            str *= "\n  $main replaces $(inline_repr(sorted[2:end]))"
        end
    end
    verbose && @info str

    # 4. Apply substitution to eqs, obs (minus old aliases), and states; reinsert canonical aliases.
    eqs_new = Symbolics.substitute.(eqs, Ref(alias_subs))

    obseqs_new = let
        obseqs_new = copy(obseqs)
        obseqs_new = deleteat!(obseqs_new, alias_idx)
        obseqs_new = Symbolics.substitute.(obseqs_new, Ref(alias_subs))
        _insert_sorted!(obseqs_new, new_alias_obs)
        obseqs_new
    end

    states_new = Symbolics.substitute.(states, Ref(alias_subs))
    allunique(states_new) || error("Alias elimination resulted in duplicate state names: $(inline_repr(states_new))! This should never happen.")

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
    reduce_equations(eqs, obseqs, states; outputs=[], ff_inputs=Set(), verbose)

Reduce implicit algebraic equations by solving them symbolically and moving solutions
into the observation list.

Iteratively calls `_match_and_solve` until no further reductions are possible. Each round
uses a Hungarian-algorithm bipartite match to pair equations with states, then calls
`_solve_plan` to decide which matched pairs to actually solve.

Feed-forward detection
-----------------------
`_solve_plan` identifies all simple dependency paths from output-state matches (output-sinks)
to `ff_input`-dependent matches (FF-sources).  Solving any match along such a path would
introduce algebraic feed-forward from an input to an output through the observation chain.

A greedy minimum hitting set selects the smallest set of matches to forbid (keep as residual
constraints), with preference for `:linear_state` matches (coefficient involves another state)
to avoid numerically risky denominators. The remaining solvable matches are returned as SCCs
in topological order.
"""
function reduce_equations(eqs::Vector{Equation}, obseqs::Vector{Equation}, states::Vector{ST}; outset::Set{ST}, ff_inputs::Set{ST}, verbose)
    # scalarize equations if necesary
    if any(eq->Symbolics.isarraysymbolic(eq.lhs), eqs)
        _eqs = Equation[]
        for eq in eqs
            if Symbolics.isarraysymbolic(eq.lhs)
                append!(_eqs, Symbolics.scalarize(eq))
            else
                push!(_eqs, eq)
            end
        end
        eqs = _eqs
    end

    if length(states) != length(eqs)
        if length(states) > length(eqs)
            superflous = setdiff(Set{ST}(states), get_variables_deriv(eqs))
            if length(superflous) == length(states) - length(eqs)
                @warn "Some provided states do not appear in the equations and will be ignored: " * inline_repr(superflous)
                states = filter(s -> s ∉ superflous, states)
            else
                throw(ArgumentError("Expected same number of states and equations. Got $(length(states)) states and $(length(eqs)) equations. Please report issue and try building component with `mtkcompile=true` in the meantime."))
            end
        else
            throw(ArgumentError("Expected same number of states and equations. Got $(length(states)) states and $(length(eqs)) equations. Please report issue and try building component with `mtkcompile=true` in the meantime."))
        end
    end

    obs_unsrtd = copy(obseqs) # we'll be inserting into this, so make a copy to avoid mutating the input
    # match states include D(x) instead of x
    diffstates_in_eqs = filter(isdifferential, Symbolics.get_variables(eqs))
    diffstates_inner_map = Dict(only(s.args) => s for s in diffstates_in_eqs)
    match_states = map(s -> get(diffstates_inner_map, s, s), states)
    match_states_set = Set{ST}(match_states)
    # ordered vector: match_states first, then any extra entries from states (the non-diff originals)
    # this keeps a deterministic order while including both D(x) and x
    all_states = vcat(match_states, filter(s -> s ∉ match_states_set, states))

    # check that we don't have any obs to expand
    if !isempty(Set(obs.lhs for obs in obs_unsrtd) ∩ match_states_set)
        error("Observation LHS cannot be a state variable! Found: $(Set(obs.lhs for obs in obs_unsrtd) ∩ match_states_set)")
    end

    # normalize alg functions which are not yet 0 ~ expr
    eqs = map(eqs) do eq
        lhs = unwrap_const(eq.lhs)
        lhszero = lhs isa Number && isequal(lhs, 0)
        if !lhszero && !isdifferential(eq.lhs)
            0 ~ eq.lhs - eq.rhs
        else
            eq
        end
    end

    # iteratively match states and move solved equations to obs_unsrtd
    while !isempty(eqs)
        before = hash((eqs, obs_unsrtd, match_states))
        eqs, obs_unsrtd, match_states = _match_and_solve(eqs, obs_unsrtd, match_states, all_states; outset, ff_inputs, verbose)
        before == hash((eqs, obs_unsrtd, match_states)) && break
    end

    # last set: bring solved diffeqs back to eqs
    rmidx = Int[]
    for (i, eq) in pairs(obs_unsrtd)
        if isdifferential(eq.lhs)
            pushfirst!(eqs, eq.lhs ~ eq.rhs)
            pushfirst!(match_states, only(eq.lhs.args))
            push!(rmidx, i)
        end
    end
    deleteat!(obs_unsrtd, rmidx)

    if verbose
        str = "Left with $(length(eqs)) equations after reduction:\n"
        str *= multiline_repr(match_states .=> eqs, prefix="  ")
        @info str
    end
    nodiff(sts) = map(s -> isdifferential(s) ? only(s.args) : s, sts)
    return eqs, _topological_sort(obs_unsrtd), nodiff(match_states)
end

function _match_and_solve(eqs, obs_unsrtd, match_states::Vector, all_states::Vector; outset, ff_inputs, verbose)
    non_extended_states = length(match_states)
    coeff, has_input, match_states, extra_match_states = _build_coeff_mat(eqs, obs_unsrtd, match_states, all_states; ff_inputs)

    if verbose
        str = "Trying to match $(length(eqs)) equations to $(length(match_states)) states for reduction"
        if non_extended_states < length(match_states)
            str *= " (including $(length(match_states) - non_extended_states) extended diff states)"
        end
        str *= "\nStates: " * inline_repr(match_states)
        str *= "\n" * multiline_repr(eqs, prefix="  ")
        @info str
    end

    solvable_matches, bk_matches = _match_equations_to_states(coeff)
    @assert length(solvable_matches) + length(bk_matches) == length(match_states)

    # Reorder states so that sorted_sts[i] is the state matched to equation i.
    # Also permute coeff columns consistently so coeff_sorted[i,j] still means
    # "equation i depends on sorted_sts[j]".
    state_perm = [-1 for i in 1:length(match_states)]
    for (eqi, si) in Iterators.flatten((solvable_matches, bk_matches))
        state_perm[eqi] = si
    end
    sorted_match_sts = match_states[state_perm]
    coeff_sorted = coeff[:, state_perm]

    # Solvable indices: rows whose matched state has a linear coefficient.
    # After reordering, match i is always at diagonal position i.
    solvable_eq_idx = sort!([eqi for (eqi, _) in solvable_matches])

    if verbose
        str = "Found $(length(solvable_matches)) solvable matches for $(length(eqs)) eqs"
        states_eqs_matching = OrderedDict(sorted_match_sts[i] => eqs[i] for i in solvable_eq_idx)
        str *= "\n" * multiline_repr(states_eqs_matching, prefix="  ")
    end

    # assert that all differential states have a linear match
    diffidx = findall(isdifferential, sorted_match_sts)
    linear_idx_set = Set(solvable_eq_idx)
    if !(Set(diffidx) ⊆ linear_idx_set)
        unmatched = sorted_match_sts[setdiff(diffidx, linear_idx_set)]
        verbose && @info str
        error("Could not match all differential states to equations! Unmatched differentials:\n" *
              multiline_repr(unmatched, prefix="  "))
    end

    # similar to has_input property for equations we need is_output for states which influence output states
    output_like = symbolic_dependencies(obs_unsrtd, outset)
    is_output = [s ∈ output_like for s in sorted_match_sts]
    sccs = _solve_plan(solvable_eq_idx, coeff_sorted, sorted_match_sts, has_input, is_output; verbose)

    solved_eq_idx = Int[]

    for (gi, scc) in enumerate(sccs)
        # scc is a Vector of indices into linear_matches; linear_matches[k] = (i,i)
        scc_idx   = [solvable_eq_idx[k] for k in scc]
        scc_sts_i = sorted_match_sts[scc_idx]
        expander = selective_expander(obs_unsrtd, scc_sts_i; expand_diffs=true)
        scc_eqs_i = expander(eqs[scc_idx])

        try
            sol = _symbolic_linear_solve_clean(scc_eqs_i, scc_sts_i)
            append!(obs_unsrtd, scc_sts_i .~ sol)
            append!(solved_eq_idx, scc_idx)
            if verbose
                if length(scc) == 1
                    str *= "\nSolved group $gi for $(only(scc_sts_i)) (trivial)"
                else
                    str *= "\nSolved group $gi for $(inline_repr(scc_sts_i))"
                    str *= "\n" * multiline_repr(OrderedDict(scc_sts_i .=> sol), prefix="  ")
                end
            end
        catch err
            if verbose
                str *= "\n\nFailed to solve group $gi for $(inline_repr(scc_sts_i)): $err"
            end
        end
    end
    verbose && @info str

    sort!(solved_eq_idx)
    deleteat!(eqs, solved_eq_idx)
    deleteat!(sorted_match_sts, solved_eq_idx)

    ####
    #### Handle solved extra states (conflicting solutions for D(x) vs x)
    ####
    solved_extra_states = setdiff(extra_match_states, sorted_match_sts)
    sorted_match_sts = filter(s -> s ∉ Set{ST}(extra_match_states), sorted_match_sts)
    if !isempty(solved_extra_states)
        if verbose
            str = ""
        end
        for s in solved_extra_states
            # find conclicting obs
            diffob_idx = findfirst(obs -> isdifferential(obs.lhs) && isequal(only(obs.lhs.args), s), obs_unsrtd)
            algob_idx = findfirst(obs -> isequal(obs.lhs, s), obs_unsrtd)
            if isnothing(diffob_idx) || isnothing(algob_idx)
                error("Expected to find both diff and alg obs for solved extended state $s, but found:\n" *
                      (isnothing(diffob_idx) ? "  No diff obs\n" : "  Diff obs: $(obs_unsrtd[diffob_idx])\n") *
                      (isnothing(algob_idx) ? "  No alg obs\n" : "  Alg obs: $(obs_unsrtd[algob_idx])\n"))
            end
            # we get rid of the dif obs
            diffeq = obs_unsrtd[diffob_idx]
            algeq = obs_unsrtd[algob_idx]
            deleteat!(obs_unsrtd, diffob_idx)
            diffop = operation(only(diffeq.lhs))

            neweq = 0 ~ Symbolics.expand_derivatives(diffop(algeq.rhs)) - diffeq.rhs
            push!(eqs, neweq)

            if verbose
                str *= "Solving exteded state lead to conflicting solutions for $s"
                str *= "\n - Diff obs: $diffeq"
                str *= "\n - Alg obs:  $algeq"
                str *= "\n=> New eq:   $neweq"
                @info str
            end
        end
    end

    ####
    #### Substitute known differentials
    ####
    required_diffs = filter(isdifferential, Symbolics.get_variables(eqs))
    if !isempty(required_diffs)
        diffop = operation(first(required_diffs))
        known_diffs = Dict{ST,ST}(obs.lhs => obs.rhs for obs in obs_unsrtd if isdifferential(obs.lhs))
        if !(required_diffs ⊆ keys(known_diffs))
            for obseq in obs_unsrtd
                isdifferential(obseq.lhs) && continue
                known_diffs[diffop(obseq.lhs)] = Symbolics.expand_derivatives(diffop(obseq.rhs))
            end
        end
        if !(required_diffs ⊆ keys(known_diffs))
            throw(RHSDifferentialsError([repr(only(d.args)) for d in setdiff(remaining_diffs, keys(known_diffs))]))
        elseif verbose
            @info "Substitute known differentials: " * multiline_repr(known_diffs)
        end
        for (i, eq) in pairs(eqs)
            eqs[i] = fixpoint_sub(eq, known_diffs)
        end
    end

    return eqs, obs_unsrtd, sorted_match_sts
end

"""
    _solve_plan(solvable_eq_idx, coeff, sorted_sts, has_input, is_output) -> sccs

Determine which solvable matched pairs should actually be solved and return them as SCCs
in dependency-first topological order.

Two dependency graphs are built: `g_solvable` (solvable dependencies between matched
pairs) and `g_unsolvable` (nonlinear dependencies). An initial SCC decomposition of
`g_solvable` is computed, then combined with BFS over `g_unsolvable` to find all paths
from output-sink matches (states that influence outputs) to FF-source matches (equations
that depend on `ff_inputs`).

A greedy minimum hitting set selects the smallest set of matches to forbid along those
paths. Candidates are ranked by a three-level cost:
1. `:linear_state` diagonal coefficient (denominator would involve another state → numerically risky to solve).
2. Output state (keeping an output as a residual is preferable to solving it algebraically along an FF path).
3. Number of `:unsolvable` (nonlinear) appearances in other equations (more appearances → riskier to solve).
SCCs are then recomputed on the non-forbidden subgraph and returned in dependency-first
topological order.
"""
function _solve_plan(solvable_eq_idx, coeff, sorted_sts, has_input, is_output; verbose)
    n = length(solvable_eq_idx)
    n == 0 && return Vector{Vector{Int}}[]

    g_solvable   = Graphs.SimpleDiGraph(n)
    g_unsolvable = Graphs.SimpleDiGraph(n)
    for (i, r) in enumerate(solvable_eq_idx)
        for (j, c) in enumerate(solvable_eq_idx)
            i == j && continue
            type = coeff[r, c]
            if type !== :none && type !== :unsolvable # any connection
                Graphs.add_edge!(g_solvable, i, j)
            elseif type == :unsolvable
                Graphs.add_edge!(g_unsolvable, i, j)
            end
        end
    end

    sccs = Graphs.strongly_connected_components_tarjan(g_solvable)

    ff_sources = Set{Int}(i for (i, r) in enumerate(solvable_eq_idx) if has_input[r])
    out_sinks  = Set{Int}(i for (i, c) in enumerate(solvable_eq_idx) if is_output[c])

    forbidden_connections = Dict{Int, Set{Int}}() # src → [dst1, dst2, ...]
    for src in out_sinks
        forbidden_connections[src] = copy(ff_sources)
    end

    # edges go from dependent -> dependency
    # the SCC is ordered from dependencies to dependents (i.e. just backward links)
    # so we need to find **forward** links via g_unsolvable: walk backward through SCCs
    # and check if any unsolvable dependency lands in an already-visited (upstream) SCC
    upstream = Set{Int}()
    for scc in Iterators.reverse(sccs)
        for dependency in Graphs.BFSIterator(g_unsolvable, scc; neighbors_type=Graphs.outneighbors)
            if dependency in upstream
                out = get!(forbidden_connections, dependency) do
                    Set{Int}()
                end
                # println("Equation $(first(scc)) depends non-linear on eq. $(dependency), which is solved after")
                push!(out, first(scc))
            end
        end
        union!(upstream, scc)
    end

    # return if no need to break cycles
    isempty(forbidden_connections) && return sccs

    pathes_to_break = Set{Vector{Int}}()
    for (src, dsts) in forbidden_connections
        for p in Graphs.all_simple_paths(g_solvable, src, dsts)
            push!(pathes_to_break, p)
        end
    end

    if verbose
        str = "Detected algebraic pathes which needs to be cut (avoid feed-forward from input to output or prevent nonlinear loop):"
        for p in pathes_to_break
            str *= "\n  " * inline_repr(sorted_sts[p])
        end
    end

    # return if no cycles detected
    isempty(pathes_to_break) && return sccs

    # Classify each match by diagonal coefficient type.
    # Lower cost = better candidate to keep as residual (tearing variable).
    match_cost = map(solvable_eq_idx, is_output) do r, io
        # prefer to break linear_state matches to avoid division by potential 0
        prio1 = coeff[r, r] == :linear_state ? 0 : 1
        # second: prefer states which are outputs
        prio2 = io ? 0 : 1
        # third: prefer states which don't appear nonlinear in other equations
        prio3 = count(c -> c == :unsolvable, view(coeff, :, r))
        (prio1, prio2, prio3)
    end

    # Greedy minimum hitting set: build inverted index node -> path indices
    node_to_paths = Dict{Int, Set{Int}}()
    for (i, path) in enumerate(pathes_to_break)
        for node in path
            push!(get!(Set{Int}, node_to_paths, node), i)
        end
    end

    remaining_paths = Set(1:length(pathes_to_break))
    chosen_nodes = Int[]

    while !isempty(remaining_paths)
        best_node = argmin(keys(node_to_paths)) do i
            hits = length(intersect(node_to_paths[i], remaining_paths))
            cost1, cost2, cost3 = match_cost[i]
            (-hits, cost1, cost2, cost3, i)  # i as tiebreaker: deterministic across Julia versions
        end
        push!(chosen_nodes, best_node)
        setdiff!(remaining_paths, node_to_paths[best_node])
    end

    # Reverse greedy cleanup: remove redundant chosen nodes
    for k in length(chosen_nodes):-1:1
        node = chosen_nodes[k]
        others = [n for n in chosen_nodes if n != node]
        isempty(others) && continue
        other_coverage = union((node_to_paths[n] for n in others)...)
        if issubset(node_to_paths[node], other_coverage)
            deleteat!(chosen_nodes, k)
        end
    end

    forbidden_matches = [i in Set(chosen_nodes) for i in 1:n]

    if verbose
        broken = sorted_sts[solvable_eq_idx[forbidden_matches]]
        str *= "\n => break by keep $(inline_repr(broken)) in equations."
    end

    non_forbidden = [i for i in 1:n if !forbidden_matches[i]]

    g_solvable_subset = g_solvable[non_forbidden]
    _sccs = Graphs.strongly_connected_components_tarjan(g_solvable_subset)

    # map back to original indices
    sccs = [[non_forbidden[k] for k in scc] for scc in _sccs]

    if verbose
        str *= "\n - Decomposed matches into $(length(sccs)) solvable groups"
        @info str
    end

    return sccs
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
        _insert_sorted!(obseqs, eq)
    end
    nothing
end
function _insert_sorted!(obseqs, eq::Equation)
    obssym = [e.lhs for e in obseqs]
    vars = get_variables_deriv(eq.rhs)
    idx = findlast(sym -> sym ∈ vars, obssym)
    last_dependency = isnothing(idx) ? 0 : idx
    insert!(obseqs, last_dependency+1, eq)
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
end

"""
Dependency type of an equation w.r.t. a state variable:
  :none         — state does not appear in the equation
  :explicit     — linear in the state; coefficient is exactly ±1 (unit coefficient, no division needed)
  :linear_const — linear in the state; coefficient is free of states (constant/parameter-only)
  :linear_state — linear in the state; coefficient involves other states (risky denominator)
  :unsolvable   — unsolvable in the state (e.g. x^2, sin(x))
"""
function _build_coeff_mat(lineqs, obseqs, match_states::Vector, all_states::Vector; ff_inputs=Set())
    @assert length(lineqs) == length(match_states)
    all_states_set = Set{ST}(all_states)
    needs_extension = false
    # Use all_states (ordered Vector) to derive extended states — this ensures that
    # even after differential states are solved and removed from match_states, their
    # inner algebraic states remain available for the extension codepath.
    match_states_set = Set{ST}(match_states)
    extended_states = ST[only(s.args) for s in all_states if isdifferential(s) && only(s.args) ∉ match_states_set]
    extended_states_set = Set{ST}(extended_states)

    N = length(match_states) + length(extended_states)
    # we build coeff with more rows ("fake" equations which do not exist)
    coeff = fill(:none, N, N)

    extended_match_states = vcat(match_states, extended_states)
    columns = Dict{ST,Int}(s => i for (i, s) in enumerate(extended_match_states))
    linexpanders = Dict{ST,LinearExpander}()
    # we need to exapnd the eq in everything to uncover both ff and linear_state dependencies
    expand_obs = selective_expander(obseqs, Any)
    has_input = falses(N)

    function get_coeff_type(eq, s)
        linex = get!(linexpanders, s) do
            LinearExpander(s)
        end
        a, b, lin = linex(eq)
        if !lin
            return :unsolvable
        elseif isequal(unwrap_const(a), 0)
            return :none
        else
            a_val = unwrap_const(a)
            if a_val isa Number && (isequal(a_val, 1) || isequal(a_val, -1))
                return :explicit
            else
                coeff_vars = get_variables_deriv(a_val)
                return isdisjoint(coeff_vars, all_states_set) ? :linear_const : :linear_state
            end
        end
    end

    for (i, eq) in enumerate(lineqs)
        eq = expand_obs(eq)
        allsyms = get_variables(eq)
        has_input[i] = !isdisjoint(allsyms, ff_inputs)
        syms = Set{ST}(extended_match_states ∩ allsyms)
        if any(isdifferential, syms)
            # for differentials we only allow solving for diff states
            for s in syms
                if isdifferential(s)
                    coeff[i, columns[s]] = get_coeff_type(eq, s)
                else
                    coeff[i, columns[s]] = :unsolvable
                end
            end
        elseif allsyms ⊆ extended_states_set
            # this is a *last resort*, we onlye extend when we have no choice to match a different state
            for s in syms
                t = get_coeff_type(eq, s)
                if t ∈ (:explicit, :linear_const, :linear_state)
                    needs_extension = true
                end
                coeff[i, columns[s]] = t
            end
        else # normal equation
            for s in syms
                coeff[i, columns[s]] = get_coeff_type(eq, s)
            end
        end
    end
    if needs_extension
        # get rid of "unused" extended states
        mask = map(1:N) do i
            i <= length(match_states) || any(c -> c !== :none, view(coeff, :, i))
        end
        coeff = coeff[mask, mask]
        extended_match_states = extended_match_states[mask]
        # make sure that "fake" equations get taken by extended states if possible
        coeff[length(lineqs)+1:end, length(match_states)+1:end] .= :nonlinear
    else
        # don't keep the extended stuff
        coeff = coeff[1:length(lineqs), 1:length(match_states)]
    end

    if needs_extension
        # Use filter on the ordered Vector (not Set intersection) so extra_match_states
        # preserves the deterministic order of extended_match_states.
        coeff, has_input, extended_match_states, filter(s -> s ∈ extended_states_set, extended_match_states)
    else
        coeff, has_input, match_states, ST[]
    end
end

"""
    _match_equations_to_states(coeff) -> (solvable, bookkeeping)

Single-pass min-cost maximum-cardinality bipartite matching via the Hungarian
algorithm.  The cost matrix encodes a strict five-level priority hierarchy:

| edge type       | cost       | meaning                                   |
|-----------------|------------|-------------------------------------------|
| `:explicit`     | 0          | coefficient is ±1; trivial solve          |
| `:linear_const` | 1          | safe to solve; constant coefficient       |
| `:linear_state` | 2          | solvable but coefficient involves states  |
| `:unsolvable`   | 2(n+1)+1   | bookkeeping only; tracks dependency       |
| `:none`         | C_nl²      | last resort; state absent from eq         |

where `n = size(coeff, 1)` (always square).  The sentinel gaps guarantee that
the Hungarian solution maximises solvable-match cardinality first, then
prefers `:explicit` over `:linear_const` over `:linear_state`, then prefers
`:unsolvable` over `:none` for bookkeeping — all in one call.

Both return values are `Vector{Tuple{Int,Int}}` of `(row, col)` pairs.
Solvable pairs have `cost ≤ 2`; bookkeeping pairs have `cost > 2`.
"""
function _match_equations_to_states(coeff)
    n = size(coeff, 1)  # always square
    C_nl   = 2*(n + 1) + 1   # > n*2 = max total solvable cost
    C_none = C_nl^2           # > n*C_nl = max total unsolvable cost

    cost = Matrix{Int}(undef, n, n)
    for i in 1:n, j in 1:n
        cost[i, j] = @match coeff[i, j] begin
            :explicit     => 0
            :linear_const => 1
            :linear_state => 2
            :unsolvable    => C_nl
            _             => C_none
        end
    end

    assignment, _ = hungarian(cost)

    solvable    = Tuple{Int,Int}[]
    bookkeeping = Tuple{Int,Int}[]
    for i in 1:n
        j = assignment[i]
        push!(cost[i, j] <= 2 ? solvable : bookkeeping, (i, j))
    end
    assignment, total_cost = hungarian(cost)
    return solvable, bookkeeping
end



function get_alias(eq)
    vars = get_variables_deriv(eq)
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

function simplify_with_mtkcompile(_sys, allinputs, alloutputs; verbose)
    missing_inputs = Set{ST}()
    sys = if ModelingToolkitBase.iscomplete(_sys)
        deepcopy(_sys)
    else
        _openinputs = setdiff(allinputs, Set(parameters(_sys)))
        all_eq_vars = mapreduce(get_variables_deriv, union, full_equations(_sys), init=Set{ST}())
        if !(_openinputs ⊆ all_eq_vars)
            missing_inputs = setdiff(_openinputs, all_eq_vars)
            verbose && @warn "The specified inputs ($missing_inputs) do not appear in the equations of the system!"
            _openinputs = setdiff(_openinputs, missing_inputs)
        end

        implicit_outputs = setdiff(alloutputs, all_eq_vars)
        if !isempty(implicit_outputs)
            throw(
                ArgumentError("The outputs $(getname.(implicit_outputs)) do not appear in the equations of the system! \
                    Try to to make them explicit using the keyword `assume_io_coupling` or the more manual `implicit_output`\n" *
                    NetworkDynamics.implicit_output_docstring)
            )
        end

        verbose && @info "Simplifying system with inputs $(inline_repr(_openinputs)) and outputs $(inline_repr(alloutputs))"
        try
            mtkcompile(_sys; inputs=_openinputs, outputs=alloutputs, simplify=false)
        catch e
            if e isa ModelingToolkitBase.ExtraEquationsSystemException
                msg = "The system could not be compiled because of extra equations! \
                    Sometimes, this can be related to fully implicit output equations. \
                    Check `@doc implicit_output` for more information."
                throw(ArgumentError(msg))
            end
            rethrow(e)
        end
    end

    # extract the main equations and observed equations
    eqs::Vector{Equation} = equations(sys)
    obseqs_sorted::Vector{Equation} = copy(observed(sys))
    fix_metadata!(eqs, sys);
    fix_metadata!(obseqs_sorted, sys);

    allparams = parameters(sys) # contains inputs!
    # mtkcompile/complete calls remove_bound_parameters_from_ps which removes params
    # lets push those bindings to the obseqs
    bps = ModelingToolkitBase.bound_parameters(sys)
    if !isempty(bps)
        bindings = ModelingToolkitBase.bindings(sys)
        newobs = [bp ~ bindings[bp] for bp in bps]
        _insert_sorted!(obseqs_sorted, newobs)
    end

    @argcheck allinputs ⊆ Set(allparams) ∪ missing_inputs
    params = setdiff(allparams, Set(allinputs))

    # create states vector in correct ordering
    states = match_diff_states(eqs, unknowns(sys))

    sys, eqs, obseqs_sorted, states, params
end

function simplify_without_mtkcompile(_sys, allinputs, alloutputs; verbose, ff_to_constraint)
    sys = ModelingToolkitBase.complete(_sys)

    eqs::Vector{Equation} = equations(sys)
    obseqs_sorted::Vector{Equation} = copy(observed(sys))

    allparams = parameters(sys) # contains inputs!
    # mtkcompile/complete calls remove_bound_parameters_from_ps which removes params
    # lets push those bindings to the obseqs
    bps = ModelingToolkitBase.bound_parameters(sys)
    if !isempty(bps)
        bindings = ModelingToolkitBase.bindings(sys)
        newobs = [bp ~ bindings[bp] for bp in bps]
        newobssorted = _topological_sort(newobs) # ensure new obs are in correct order
        append!(obseqs_sorted, newobssorted)
    end

    # maybe inputs have been given as parameters
    params = setdiff(allparams, Set{ST}(allinputs))

    # don't consider inputs states; preserve unknowns(sys) order for deterministic matching
    inputs_set = Set{ST}(allinputs)
    states = filter(s -> s ∉ inputs_set, unknowns(sys))

    # check if we need to block input ff
    ff_inputs = ff_to_constraint ? inputs_set : Set{ST}()

    # reduce equations (solves for differentials and algebraic variables)
    eqs, obseqs, states = reduce_equations(
        eqs, obseqs_sorted, states;
        outset=Set{ST}(alloutputs), ff_inputs, verbose
    )

    sys, eqs, obseqs, states, params
end

isdifferential(s) = iscall(unwrap(s)) && operation(unwrap(s)) isa Differential

function selective_expander(obseqs, targets; expand_diffs=false)
    if isempty(obseqs) || (targets != Any && isempty(targets))
        return identity
    end

    if targets == Any
        obs_subset = _topological_sort(obseqs)
    else
        sinks = [o.lhs for o in obseqs if !isdisjoint(get_variables(o.rhs), targets)]
        g, map = _dependency_graph(obseqs)
        sinks_i = [map[s] for s in sinks]
        necessary = collect(Int, Graphs.BFSIterator(g, sinks_i; neighbors_type=Graphs.inneighbors))

        obs_subset = obseqs[necessary]
        perm = reverse(Graphs.topological_sort(g[necessary]))
        permute!(obs_subset, perm)
    end

    substitutions = OrderedDict{ST,ST}()
    for obs in obs_subset
        substitutions[obs.lhs] = substitute(obs.rhs, substitutions)
    end

    if expand_diffs
        diff_subs = Dict{ST,ST}(o.lhs => o.rhs for o in obseqs if isdifferential(o.lhs))
        # first: expand diff in diff (non sorted)
        for (lhs, rhs) in diff_subs
            diff_subs[lhs] = fixpoint_sub(rhs, diff_subs)
        end
        # second: add them to selective substitutions and expand selective (sorted)
        for (lhs, rhs) in diff_subs
            substitutions[lhs] = substitute(rhs, substitutions)
        end
    end

    (term) -> begin
        substitute(term, substitutions)
    end
end

"""
    symbolic_dependencies(obseqs, states)

Given a list of states, return a set of superset of the provided states that also
includes any "observed" state from obseqs which transitively influences the provided state.
"""
function symbolic_dependencies(obseqs, states::Set)
    isempty(obseqs) && return Set{ST}(states)
    g, s_to_i_map = _dependency_graph(obseqs)
    srcs = Int[s_to_i_map[s] for s in states if haskey(s_to_i_map, s)]
    transitive = Set{ST}(states)
    for i in Graphs.BFSIterator(g, srcs)
        union!(transitive, get_variables(obseqs[i]))
    end
    transitive
end

function _topological_sort(obseqs)
    isempty(obseqs) && return obseqs
    g, s_to_i_map = _dependency_graph(obseqs)
    sorted_indices = Graphs.topological_sort(g)
    obseqs[reverse(sorted_indices)]
end

function _dependency_graph(obseqs)
    s_to_i_map = Dict{ST, Int}(eq.lhs => i for (i, eq) in enumerate(obseqs))
    g = Graphs.SimpleDiGraph(length(s_to_i_map))
    for (i, obs) in enumerate(obseqs)
        vars = get_variables(obs.rhs)
        for v in vars
            if haskey(s_to_i_map, v)
                Graphs.add_edge!(g, i, s_to_i_map[v])
            end
        end
    end
    (g, s_to_i_map)
end
