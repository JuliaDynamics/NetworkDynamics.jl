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
function pick_best_alias_names(eqs, obseqs, states, outputs; verbose)
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
        type = eq_type(eq)
        type[1] == :explicit_diffeq && push!(diffstateset, type[2])
    end
    outset = Set(outputs)
    sortf = function(s)
        prio = if s ∈ diffstateset
            0
        elseif s ∈ outset
            1
        else
            2
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
    eqs_new = Symbolics.substitute_in_deriv.(eqs, Ref(alias_subs))

    obseqs_new = let
        obseqs_new = copy(obseqs)
        obseqs_new = deleteat!(obseqs_new, alias_idx)
        obseqs_new = Symbolics.substitute.(obseqs_new, Ref(alias_subs))
        _insert_sorted!(obseqs_new, new_alias_obs)
        obseqs_new
    end

    states_new = Symbolics.substitute.(states, Ref(alias_subs))
    allunique(states_new) || error("Alias elimination resulted in duplicate state names: $(states_new)! This should never happen.")

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
    state_set = Set(states)

    if verbose
        str = "Trying to match $(length(eqs)) equations to $(length(states)) states for reduction"
        str *= "\nStates: " * inline_repr(states)
        str *= "\n" * multiline_repr(eqs, prefix="  ")
        @info str
    end

    @assert length(states) == length(eqs)

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

    # Expand state_set to include inner differential variables (e.g. θ for D(θ)).
    # These are genuine dynamic variables, so coefficients like cos(θ) or sin(θ) are
    # state-dependent (:linear_state) rather than constant (:linear_const).
    diff_inner_set = Set(only(s.args) for s in states if isdifferential(s))
    coeff_state_set = state_set ∪ diff_inner_set

    coeff, has_ff = _build_coeff_mat(eqs, states, coeff_state_set; ff_inputs)
    solvable_matches, bk_matches = _match_equations_to_states(coeff)
    @assert length(solvable_matches) + length(bk_matches) == length(states)

    # Reorder states so that sorted_sts[i] is the state matched to equation i.
    # Also permute coeff columns consistently so coeff_sorted[i,j] still means
    # "equation i depends on sorted_sts[j]".
    state_perm = [-1 for i in 1:length(states)]
    for (eqi, si) in Iterators.flatten((solvable_matches, bk_matches))
        state_perm[eqi] = si
    end
    sorted_sts = states[state_perm]
    coeff_sorted = coeff[:, state_perm]

    # Solvable indices: rows whose matched state has a linear coefficient.
    # After reordering, match i is always at diagonal position i.
    linear_idx = sort!([eqi for (eqi, _) in solvable_matches])

    if verbose
        str = "Found $(length(linear_idx)) linear matches out of $(length(states))"
        states_eqs_matching = OrderedDict(sorted_sts[i] => eqs[i] for i in linear_idx)
        str *= "\n" * multiline_repr(states_eqs_matching, prefix="  ")
        @info str
    end

    # assert that all differential states have a linear match
    diffidx = findall(isdifferential, sorted_sts)
    linear_idx_set = Set(linear_idx)
    if !(Set(diffidx) ⊆ linear_idx_set)
        unmatched = sorted_sts[setdiff(diffidx, linear_idx_set)]
        error("Could not match all differential states to equations! Unmatched differentials:\n" *
              multiline_repr(unmatched, prefix="  "))
    end

    # _solve_plan expects matches as (row, col) pairs; since row==col after reordering, zip with self.
    linear_matches = [(i, i) for i in linear_idx]
    sccs, forbidden_matches = _solve_plan(linear_matches, coeff_sorted, sorted_sts, has_ff, outset)

    if verbose
        n_forbidden = sum(forbidden_matches)
        if n_forbidden > 0
            forbidden_sts = [sorted_sts[linear_idx[k]] for k in findall(forbidden_matches)]
            @info "Skipping $(inline_repr(forbidden_sts)) — on input→output feed-forward path"
        end
        @info "Decomposed matches into $(length(sccs)) solvable groups"
    end

    solved_local = Int[]
    newobs = Equation[]
    primary_diff_eq = Dict{ST,ST}()

    for (gi, scc) in enumerate(sccs)
        # scc is a Vector of indices into linear_matches; linear_matches[k] = (i,i)
        scc_idx   = [linear_idx[k] for k in scc]
        scc_eqs_i = eqs[scc_idx]
        scc_sts_i = sorted_sts[scc_idx]

        diffstate = any(isdifferential, scc_sts_i)
        if !diffstate
            try
                sol = _symbolic_linear_solve_clean(scc_eqs_i, scc_sts_i)
                append!(newobs, scc_sts_i .~ sol)
                append!(solved_local, scc_idx)
                if verbose
                    str = "Solved group $gi for $(inline_repr(scc_sts_i))"
                    str *= "\nSolution:"
                    str *= "\n" * multiline_repr(OrderedDict(scc_sts_i .=> sol), prefix="  ")
                    @info str
                end
            catch err
                verbose && @info "Failed to solve group $gi for $(inline_repr(scc_sts_i)): $err"
            end
        else
            @assert all(isdifferential, scc_sts_i) "Got mixed diff/non-diff SCC, unhandled."
            sol = _symbolic_linear_solve_clean(scc_eqs_i, scc_sts_i)
            eqs[scc_idx] .= scc_sts_i .~ sol
            for (s, e) in zip(scc_sts_i, sol)
                primary_diff_eq[s] = e
            end
            if verbose
                str = "Solved group $gi for diff states $(inline_repr(scc_sts_i))"
                str *= "\nSolution:"
                str *= "\n" * multiline_repr(OrderedDict(scc_sts_i .=> sol), prefix="  ")
                @info str
            end
        end
    end

    isempty(newobs) && return eqs, obseqs, sorted_sts

    _insert_sorted!(obseqs, newobs)
    sort!(solved_local)
    deleteat!(eqs, solved_local)
    deleteat!(sorted_sts, solved_local)
    if verbose
        str = "Left with equations after reduction:\n"
        str *= multiline_repr(sorted_sts .=> eqs, prefix="  ")
        @info str
    end
    return eqs, obseqs, sorted_sts
end

"""
    _solve_plan(matches, coeff, sorted_sts, has_ff, outset) -> (sccs, forbidden_matches)

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
function _solve_plan(matches, coeff, sorted_sts, has_ff, outset)
    n = length(matches)
    n == 0 && return (Vector{Vector{Int}}[], falses(0))

    # Step 1: Dependency graph on matched pairs.
    # Edge i→j means "equation of match i has a nonzero coeff for state of match j",
    # i.e. match i depends on match j and must be solved after it.
    # We include :nonlinear dependencies too, so that equations like
    # `x ~ atan(a, b)` are correctly ordered after `a` and `b` in the SCC output.
    g = Graphs.SimpleDiGraph(n)
    for (i, (r, _)) in enumerate(matches)
        for (j, (_, c)) in enumerate(matches)
            i == j && continue
            if coeff[r, c] !== :none
                Graphs.add_edge!(g, i, j)
            end
        end
    end

    # Early exit when FF detection is unnecessary.
    if isempty(outset) || all(!, has_ff)
        sccs = Graphs.strongly_connected_components_tarjan(g)
        return (sccs, falses(n))
    end

    # Step 2: Classify each match by the coefficient type of its diagonal entry.
    match_type = fill(:linear_const, n)
    for (mi, (r, c)) in enumerate(matches)
        if coeff[r, c] === :linear_state
            match_type[mi] = :linear_state
        end
    end

    # Step 3: Identify FF-source matches (equation directly references ff_inputs)
    # and output-sink matches (matched state is an output).
    ff_sources = [i for (i, (r, _)) in enumerate(matches) if has_ff[r]]
    out_sinks  = [i for (i, (_, c)) in enumerate(matches) if sorted_sts[c] ∈ outset]

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
        obssym = [e.lhs for e in obseqs]
        vars = get_variables_deriv(eq.rhs)
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
end

"""
Dependency type of an equation w.r.t. a state variable:
  :none         — state does not appear in the equation
  :explicit     — linear in the state; coefficient is exactly ±1 (unit coefficient, no division needed)
  :linear_const — linear in the state; coefficient is free of states (constant/parameter-only)
  :linear_state — linear in the state; coefficient involves other states (risky denominator)
  :nonlinear    — nonlinear in the state (e.g. x^2, sin(x))

Returns `(coeff, has_ff)` where `has_ff[i]` is true if equation `i` directly references
an ff_input variable.  FF blocking is done in `_solve_plan`, not here.
"""
function _build_coeff_mat(lineqs, linstates, state_set; ff_inputs=Set())
    coeff = Matrix{Symbol}(undef, length(lineqs), length(linstates))
    has_ff = isempty(ff_inputs) ? falses(length(lineqs)) :
             Bool[!isempty(Set(get_variables_deriv(eq)) ∩ ff_inputs) for eq in lineqs]
    for (j, s) in enumerate(linstates)
        ex = Symbolics.LinearExpander(s)
        for (i, eq) in enumerate(lineqs)
            a, b, lin = ex(eq)
            if !lin
                coeff[i, j] = :nonlinear
            elseif isequal(unwrap_const(a), 0)
                coeff[i, j] = :none
            else
                a_val = unwrap_const(a)
                if a_val isa Number && (isequal(a_val, 1) || isequal(a_val, -1))
                    coeff[i, j] = :explicit
                else
                    coeff_vars = get_variables_deriv(a_val)
                    coeff[i, j] = isempty(coeff_vars ∩ state_set) ? :linear_const : :linear_state
                end
            end
        end
    end
    coeff, has_ff
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
| `:nonlinear`    | 2(n+1)+1   | bookkeeping only; tracks dependency       |
| `:none`         | C_nl²      | last resort; state absent from eq         |

where `n = size(coeff, 1)` (always square).  The sentinel gaps guarantee that
the Hungarian solution maximises solvable-match cardinality first, then
prefers `:explicit` over `:linear_const` over `:linear_state`, then prefers
`:nonlinear` over `:none` for bookkeeping — all in one call.

Both return values are `Vector{Tuple{Int,Int}}` of `(row, col)` pairs.
Solvable pairs have `cost ≤ 2`; bookkeeping pairs have `cost > 2`.
"""
function _match_equations_to_states(coeff)
    n = size(coeff, 1)  # always square
    C_nl   = 2*(n + 1) + 1   # > n*2 = max total solvable cost
    C_none = C_nl^2           # > n*C_nl = max total nonlinear cost

    cost = Matrix{Int}(undef, n, n)
    for i in 1:n, j in 1:n
        cost[i, j] = @match coeff[i, j] begin
            :explicit     => 0
            :linear_const => 1
            :linear_state => 2
            :nonlinear    => C_nl
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

    # states of the system are
    diffstates = filter(isdifferential, Set{ST}(get_variables(eqs)))
    diffstates_inner = [only(s.args) for s in diffstates]

    wo_inputs = setdiff(unknowns(sys), allinputs)
    just_algebraic = setdiff(wo_inputs, diffstates_inner)
    states = vcat(collect(diffstates), just_algebraic)

    # check if we need to block input ff
    ff_inputs = ff_to_constraint ? Set{ST}(allinputs) : Set{ST}()

    # reduce quations (solves for differentials and algebraic variables)
    eqs, obseqs, states = reduce_equations(
        eqs, obseqs_sorted, states;
        outset=Set{ST}(alloutputs), ff_inputs, verbose
    )

    states_nodiff = map(states) do s
        isdifferential(s) ? only(s.args) : s
    end

    _sys, eqs, obseqs, states_nodiff, params
end

isdifferential(s) = iscall(unwrap(s)) && operation(unwrap(s)) isa Differential 
