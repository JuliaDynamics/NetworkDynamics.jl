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

function observed_outputs_to_states(eqs, obseqs, states, outputs; verbose)
    outset = Set(outputs)
    out_obs = Int[]
    for (i, eq) in pairs(obseqs)
        if eq.lhs ∈ outset
            push!(out_obs, i)
        end
    end
    isempty(out_obs) && return eqs, obseqs, states

    if verbose
        str = "Promote observed outputs to states:\n"
        str *= multiline_repr(obseqs[out_obs], prefix="  ")
        @info str
    end
    new_obseqs = copy(obseqs)
    deleteat!(new_obseqs, out_obs)
    new_eqs = copy(eqs)
    new_states = copy(states)
    for i in out_obs
        eq = obseqs[i]
        push!(new_eqs, 0 ~ eq.lhs - eq.rhs)
        push!(new_states, eq.lhs)
    end
    new_eqs, new_obseqs, new_states
end

"""
    pick_best_alias_names(eqs, obseqs, states, outputs; verbose)

Post-processing pass that consolidates alias chains produced by `reduce_linear_algebraic`.

After reduction, the obs list may contain pure alias equations of the form `a ~ b` that
connect several names for the same quantity (e.g. `busbar₊u_r ~ busbar₊terminal₊u_r ~
gen₊terminal₊u_r`).  This pass:
1. Collects all such alias equations from `obseqs`.
2. Groups transitively connected aliases into clusters.
3. Picks one *main* representative per cluster — preferring differential states, then output variables, then variables
   with fewer `₊` separators (less deeply nested).
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

State preference
----------------
Columns (states) are sorted so that deeply-nested states (many `₊` in their name) and
non-output states come first in the bipartite matching.  This causes the DFS to prefer
eliminating internal/nested variables, keeping output states and short-named states in the
remaining equation system.

Feed-forward detection
-----------------------
After building SCCs of the matched-pair dependency graph, `_solve_plan` identifies which
SCCs lie on a directed path from an `ff_input`-dependent SCC (FF-source) to an output-state
SCC (output-sink) in the SCC meta-DAG.  Solving such an SCC would introduce algebraic
feed-forward from an input to an output through the observation chain.

For each FF-source → output-sink path, the algorithm walks from the output end toward the
source and forbids the first `:linear_state` SCC it encounters (one whose matched coefficient
involves another state, making division by that state necessary).  If no `:linear_state` SCC
exists on a path, the output-sink SCC itself is forbidden.  This keeps the forbidden set as
small as possible while preferring to cut at numerically risky solve points.
"""
function reduce_equations(eqs::Vector{Equation}, obseqs::Vector{Equation}, states::Vector{ST}; outset::Set{ST}, ff_inputs::Set{ST}, verbose)
    obs_unsrtd = copy(obseqs) # we'll be inserting into this, so make a copy to avoid mutating the input
    # match states include D(x) instead of x
    diffstates_in_eqs = filter(isdifferential, Symbolics.get_variables(eqs))
    diffstates_inner_map = Dict(only(s.args) => s for s in diffstates_in_eqs)
    match_states = map(s -> get(diffstates_inner_map, s, s), states)
    match_states_set = Set{ST}(match_states)
    all_states = match_states_set ∪ states # includs both D(x) and x

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
            push!(match_states, only(eq.lhs.args))
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

function _match_and_solve(eqs, obs_unsrtd, match_states, all_states; outset, ff_inputs, verbose)
    @assert length(match_states) == length(eqs)
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

    # similarlar to has_input property for equations we need is_output for states which influce output states
    output_like = symbolic_dependencies(obs_unsrtd, outset)
    is_output = [s ∈ output_like for s in sorted_match_sts]
    sccs, forbidden_matches = _solve_plan(solvable_eq_idx, coeff_sorted, sorted_match_sts, has_input, is_output)

    if verbose
        n_forbidden = sum(forbidden_matches)
        if n_forbidden > 0
            forbidden_sts = [sorted_match_sts[solvable_eq_idx[k]] for k in findall(forbidden_matches)]
            str *= "\n - Skipping $(inline_repr(forbidden_sts)) — on input→output feed-forward path"
        end
        str *= "\n - Decomposed matches into $(length(sccs)) solvable groups"
    end

    solved_eq_idx = Int[]

    for (gi, scc) in enumerate(sccs)
        # scc is a Vector of indices into linear_matches; linear_matches[k] = (i,i)
        scc_idx   = [solvable_eq_idx[k] for k in scc]
        scc_sts_i = sorted_match_sts[scc_idx]
        expander = selective_expander(obs_unsrtd, scc_sts_i)
        scc_eqs_i = expander(eqs[scc_idx])

        try
            sol = _symbolic_linear_solve_clean(scc_eqs_i, scc_sts_i)
            append!(obs_unsrtd, scc_sts_i .~ sol)
            append!(solved_eq_idx, scc_idx)
            if verbose
                str *= "\n\nSolved group $gi for $(inline_repr(scc_sts_i))"
                str *= "\n" * multiline_repr(OrderedDict(scc_sts_i .=> sol), prefix="  ")
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
    sorted_match_sts = filter(s -> s ∉ extra_match_states, sorted_match_sts)
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
    _solve_plan(matches, coeff, sorted_sts, has_input, is_output) -> (sccs, forbidden_matches)

Decompose matched pairs into SCCs of the dependency digraph (topological order,
dependencies before dependents) and return a `BitVector` of length `n` marking
individual matches that must not be solved because they lie on a feed-forward path
from an `ff_input`-dependent equation to an output state.

FF detection is performed on the **raw** match-level dependency graph before SCCs are
formed.  In the raw graph, edge `i→j` means "equation of match `i` depends on state of
match `j`", so a path `out_sink → … → ff_source` (following dependency edges from the
output end) is exactly the algebraic FF chain that must be cut.

For each such path the algorithm walks from the output-sink end toward the ff-source and
forbids the first `:linear_state` match encountered (one whose matched coefficient
involves another state, making division by that state necessary).  If no `:linear_state`
match exists on a path the output-sink match itself is forbidden.

After forbidden matches are identified, SCCs are computed only on the surviving
non-forbidden subgraph, so cycles that existed solely through a forbidden match are
already eliminated.  The returned `sccs` therefore contain only solvable matches and
are already in dependency-first topological order.
"""
function _solve_plan(solvable_eq_idx, coeff, sorted_sts, has_input, is_output)
    n = length(solvable_eq_idx)
    n == 0 && return (Vector{Vector{Int}}[], falses(0))

    # Step 1: Dependency graph on matched pairs.
    # Edge i→j means "equation of match i has a nonzero coeff for state of match j",
    # i.e. match i depends on match j and must be solved after it.
    # We include :unsolvable dependencies too, so that equations like
    # `x ~ atan(a, b)` are correctly ordered after `a` and `b` in the SCC output.
    g = Graphs.SimpleDiGraph(n)
    for (i, r) in enumerate(solvable_eq_idx)
        for (j, c) in enumerate(solvable_eq_idx)
            i == j && continue
            if coeff[r, c] !== :none
                Graphs.add_edge!(g, i, j)
            end
        end
    end

    # Early exit when FF detection is unnecessary.
    if all(!, has_input)
        sccs = Graphs.strongly_connected_components_tarjan(g)
        return (sccs, falses(n))
    end

    # Step 2: Classify each match by the coefficient type of its diagonal entry.
    match_type = fill(:linear_const, n)
    for (mi, rc) in enumerate(solvable_eq_idx)
        if coeff[rc, rc] === :linear_state
            match_type[mi] = :linear_state
        end
    end

    # Step 3: Identify FF-source matches (equation directly references ff_inputs)
    # and output-sink matches (matched state is an output).
    ff_sources = [i for (i, r) in enumerate(solvable_eq_idx) if has_input[r]]
    out_sinks  = [i for (i, c) in enumerate(solvable_eq_idx) if is_output[c]]

    # Step 4: Forbid individual matches to break FF paths.
    # In the dependency graph g (i→j = "i depends on j"), a path from out_sink to ff_source
    # means out_sink transitively depends on ff_source — the algebraic FF chain.
    # For each such path walk from the output-sink end and forbid the first :linear_state
    # match; fall back to the output-sink match itself if none is found.
    forbidden_matches = falses(n)
    for src in ff_sources
        for dst in out_sinks
            if src == dst
                forbidden_matches[src] = true
                continue
            end
            for path in Graphs.all_simple_paths(g, dst, src)
                # path[1] == dst (output-sink), path[end] == src (ff-source).
                forbid_idx = dst   # fallback: forbid the output-sink match
                for mi in path
                    if match_type[mi] === :linear_state
                        forbid_idx = mi
                        break
                    end
                end
                forbidden_matches[forbid_idx] = true
            end
        end
    end

    # Step 5: SCC decomposition on non-forbidden matches only.
    # Cycles that ran through a forbidden match are already broken.
    non_forbidden = [i for i in 1:n if !forbidden_matches[i]]
    isempty(non_forbidden) && return (Vector{Vector{Int}}[], forbidden_matches)

    new_idx = zeros(Int, n)
    for (k, i) in enumerate(non_forbidden)
        new_idx[i] = k
    end
    g_sub = Graphs.SimpleDiGraph(length(non_forbidden))
    for i in non_forbidden
        for j in Graphs.outneighbors(g, i)
            forbidden_matches[j] && continue
            Graphs.add_edge!(g_sub, new_idx[i], new_idx[j])
        end
    end

    # For our edge convention (i→j = "i depends on j"), strongly_connected_components_tarjan
    # already returns SCCs with dependencies before dependents — no reversal needed.
    sccs_sub = Graphs.strongly_connected_components_tarjan(g_sub)

    # Remap subgraph SCC indices back to original match indices.
    sccs = [[non_forbidden[k] for k in scc] for scc in sccs_sub]
    return (sccs, forbidden_matches)
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
function _build_coeff_mat(lineqs, obseqs, match_states, all_states; ff_inputs=Set())
    @assert length(lineqs) == length(match_states)
    needs_extension = false
    extended_states = ST[only(s.args) for s in all_states if isdifferential(s)]
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
                return isdisjoint(coeff_vars, all_states) ? :linear_const : :linear_state
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
        coeff, has_input, extended_match_states, extended_states_set ∩ extended_match_states 
    else
        coeff, has_input, match_states, Set{ST}()
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

function simplify_with_mtkcompile(_sys, allinputs, alloutputs)
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
    # whose value is bound to another param (e.g. Sn = S_b) from ps, but those symbols
    # still appear in equations. Add them back so build_function can reference them and
    # InitFormulas targeting them remain valid.
    append!(allparams, ModelingToolkitBase.bound_parameters(sys))
    @argcheck allinputs ⊆ Set(allparams) ∪ missing_inputs
    params = setdiff(allparams, Set(allinputs))

    # create states vector in correct ordering
    states = match_diff_states(eqs, unknowns(sys))
 
    # pick best alias names
    eqs, obseqs, states = pick_best_alias_names(eqs, obseqs, states, alloutputs; verbose)

    sys, eqs, obseqs_sorted, states, params
end

function simplify_without_mtkcompile(_sys, allinputs, alloutputs; verbose, ff_to_constraint)
    sys = ModelingToolkitBase.complete(_sys)

    eqs::Vector{Equation} = equations(sys)
    obseqs_sorted::Vector{Equation} = copy(observed(sys))

    allparams = parameters(sys) # contains inputs!
    append!(allparams, ModelingToolkitBase.bound_parameters(sys))

    # maybe inputs have been given as parameters
    params = setdiff(allparams, Set{ST}(allinputs))

    # don't consider inputs states
    states = collect(setdiff(Set{ST}(unknowns(sys)), Set{ST}(allinputs)))

    # check if we need to block input ff
    ff_inputs = ff_to_constraint ? Set{ST}(allinputs) : Set{ST}()

    # reduce equations (solves for differentials and algebraic variables)
    eqs, obseqs, states = reduce_equations(
        eqs, obseqs_sorted, states;
        outset=Set{ST}(alloutputs), ff_inputs, verbose
    )


    _sys, eqs, obseqs, states, params
end

isdifferential(s) = iscall(unwrap(s)) && operation(unwrap(s)) isa Differential 

function selective_expander(obseqs, targets)
    obssubs = Dict{ST, ST}(obs.lhs => obs.rhs for obs in obseqs)
    (term) -> begin
        fixpoint_sub(term, obssubs)
    end
end

"""
    symbolic_dependencies(obseqs, states)

Given a list of states, return a set of superset of the provided states that also
includes any "observed" state from obseqs which transitively influences the provided state.
"""
function symbolic_dependencies(obseqs, states)
    isempty(obseqs) && return Set{ST}(states)
    g, s_to_i_map = _dependecy_graph(obseqs)
    srcs = Int[s_to_i_map[s] for s in states if haskey(s_to_i_map, s)]
    transitive = Set{ST}(states)
    for i in Graphs.BFSIterator(g, srcs)
        union!(transitive, get_variables(obseqs[i]))
    end
    transitive
end

function _topological_sort(obseqs)
    isempty(obseqs) && return obseqs
    g, s_to_i_map = _dependecy_graph(obseqs)
    sorted_indices = Graphs.topological_sort(g)
    obseqs[reverse(sorted_indices)]
end

function _dependecy_graph(obseqs)
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
