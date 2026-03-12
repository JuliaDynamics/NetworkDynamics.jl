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
3. Picks one *main* representative per cluster — preferring output variables, then variables
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
    outset = Set(outputs)
    sortf(s) = s ∈ outset ? typemin(Int) : count('₊', string(s))
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
    reduce_linear_algebraic(eqs, obseqs, states; outputs=[], ff_inputs=Set(), verbose)

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
function reduce_linear_algebraic(eqs, obseqs, states; outputs=[], ff_inputs=Set(), verbose)
    try
        outset    = Set(outputs)
        state_set = Set(states)

        alg_idx = findall(eqs) do eq
            eq_type(eq)[1] == :implicit_algebraic
        end

        isempty(alg_idx) && return eqs, obseqs, states

        if verbose
            str = "Trying to reduce $(length(alg_idx)) implicit algebraic equations"
            str *= "\n" * multiline_repr(eqs[alg_idx], prefix="  ")
            @info str
        end

        # Sort columns (states) by preference for solving:
        # deeply-nested (many ₊) and non-output states come first so the DFS matching
        # prefers to eliminate them, leaving outputs and short-named states in the system.
        sortkey(i) = (states[i] ∈ outset ? 1 : 0, -count('₊', string(states[i])))
        sort!(alg_idx, by=sortkey)

        alg_eqs = eqs[alg_idx]
        alg_sts = states[alg_idx]

        coeff, has_ff = _build_coeff_mat(alg_eqs, alg_sts, state_set; ff_inputs)

        matches = _match_equations_to_states(coeff)

        if verbose
            str = "Found $(length(matches)) matches in those equations"
            states_eqs_matching = OrderedDict(alg_sts[m[2]] => alg_eqs[m[1]] for m in matches)
            str *= "\n" * multiline_repr(states_eqs_matching, prefix="  ")
            @info str
        end

        isempty(matches) && return eqs, obseqs, states

        sccs, forbidden = _solve_plan(matches, coeff, alg_sts, has_ff, outset)
        verbose && @info "Decomposed matches into $(length(sccs)) groups to solve"

        solved_eq_local = Int[]
        solved_st_local = Int[]
        newobs = Equation[]

        for (gi, scc) in enumerate(sccs)
            scc_rows  = [matches[i][1] for i in scc]
            scc_cols  = [matches[i][2] for i in scc]
            scc_eqs_i = alg_eqs[scc_rows]
            scc_sts_i = alg_sts[scc_cols]

            if forbidden[gi]
                verbose && @info "Group $gi: skipping $(inline_repr(scc_sts_i)) — on input→output path"
                continue
            end

            try
                sol = _symbolic_linear_solve_clean(scc_eqs_i, scc_sts_i)
                append!(newobs, scc_sts_i .~ sol)
                append!(solved_eq_local, scc_rows)
                append!(solved_st_local, scc_cols)
                if verbose
                    str = "Solved group $gi for $(inline_repr(scc_sts_i))"
                    str *= "\nSolution:"
                    str *= "\n" * multiline_repr(OrderedDict(scc_sts_i .=> sol), prefix="  ")
                    @info str
                end
            catch err
                verbose && @info "Failed to solve group $gi for $(inline_repr(scc_sts_i)): $err"
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

"""
    _solve_plan(matches, coeff, alg_sts, has_ff, outset) -> (sccs, forbidden)

Decompose matched pairs into SCCs of the dependency digraph (topological order,
dependencies before dependents) and return a `BitVector` marking SCCs that must
not be solved because they lie on a feed-forward path from an `ff_input`-dependent
equation to an output state.

For each feed-forward source → output-sink path in the SCC meta-DAG, the SCC closest
to the output that is classified as `:linear_state` is forbidden.  If no such SCC
exists on a path, the output-sink SCC itself is forbidden.
"""
function _solve_plan(matches, coeff, alg_sts, has_ff, outset)
    n = length(matches)
    n == 0 && return (Vector{Int}[], BitVector())

    # Step 1: Dependency graph on matched pairs.
    # Edge i→j means "equation of match i has a nonzero coeff for state of match j",
    # i.e. match i depends on match j and must be solved after it.
    g = Graphs.SimpleDiGraph(n)
    for (i, (r, _)) in enumerate(matches)
        for (j, (_, c)) in enumerate(matches)
            i == j && continue
            if coeff[r, c] === :linear_const || coeff[r, c] === :linear_state
                Graphs.add_edge!(g, i, j)
            end
        end
    end

    # Step 2: SCC decomposition in dependency-first order.
    # For our edge convention (i→j = "i depends on j"), strongly_connected_components_tarjan
    # already returns SCCs with dependencies before dependents — no reversal needed.
    sccs = Graphs.strongly_connected_components_tarjan(g)

    # Early exit when FF detection is unnecessary.
    if isempty(outset) || all(!, has_ff)
        return (sccs, falses(length(sccs)))
    end

    n_sccs = length(sccs)

    # Step 3: Classify each SCC.
    # :linear_state if any of its matched pairs has coeff === :linear_state; else :linear_const.
    scc_type = fill(:linear_const, n_sccs)
    for (si, scc) in enumerate(sccs)
        for mi in scc
            r, c = matches[mi]
            if coeff[r, c] === :linear_state
                scc_type[si] = :linear_state
                break
            end
        end
    end

    # Map match index → SCC index.
    match_to_scc = Dict(mi => si for (si, scc) in enumerate(sccs) for mi in scc)

    # Step 4: Meta-DAG of SCCs.
    # For each match-level edge i→j (i depends on j), add meta edge scc(j)→scc(i).
    # FF flows from the dependency scc(j) toward the dependent scc(i).
    meta_g = Graphs.SimpleDiGraph(n_sccs)
    for e in Graphs.edges(g)
        si = match_to_scc[e.src]   # dependent
        sj = match_to_scc[e.dst]   # dependency
        si == sj && continue
        Graphs.add_edge!(meta_g, sj, si)
    end

    # Step 5: Identify FF-source and output-sink SCCs.
    ff_sources = Int[]
    out_sinks  = Int[]
    for (si, scc) in enumerate(sccs)
        rows = [matches[mi][1] for mi in scc]
        cols = [matches[mi][2] for mi in scc]
        any(has_ff[r] for r in rows)            && push!(ff_sources, si)
        any(alg_sts[c] ∈ outset for c in cols)  && push!(out_sinks,  si)
    end

    # Step 6: Forbid SCCs that would introduce algebraic FF.
    # For each (ff_source, output_sink) pair, walk every path from output toward the
    # input and forbid the first :linear_state SCC encountered; fall back to the sink.
    forbidden = falses(n_sccs)
    for src in ff_sources
        for dst in out_sinks
            if src == dst
                forbidden[src] = true
                continue
            end
            for path in Graphs.all_simple_paths(meta_g, src, dst)
                # Walk path from output end (dst) toward ff-source (src).
                forbid_idx = dst   # fallback: forbid the output sink
                for scc_idx in reverse(path)
                    if scc_type[scc_idx] === :linear_state
                        forbid_idx = scc_idx
                        break
                    end
                end
                forbidden[forbid_idx] = true
            end
        end
    end

    return (sccs, forbidden)
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
    Symbolics.simplify.(sol, expand=true)
end

"""
Dependency type of an equation w.r.t. a state variable:
  :none         — state does not appear in the equation
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
                coeff_vars = Set(get_variables_deriv(unwrap_const(a)))
                coeff[i, j] = isempty(coeff_vars ∩ state_set) ? :linear_const : :linear_state
            end
        end
    end
    coeff, has_ff
end

"""
    _maximum_bipartite_matching(biadj) -> cm

Maximum bipartite matching via DFS augmenting paths (Kuhn's algorithm).
`biadj[i,j]` is true if row `i` (equation) can be matched to column `j` (state).
Returns `cm` where `cm[j]` is the row matched to column `j`, or 0 if unmatched.
"""
function _maximum_bipartite_matching(biadj::AbstractMatrix{Bool})
    nrow, ncol = size(biadj)
    cm = zeros(Int, ncol)

    function try_augment!(row, visited)
        for col in 1:ncol
            biadj[row, col] && !visited[col] || continue
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
    # Both :linear_const and :linear_state are matchable; FF blocking is in _solve_plan.
    biadj = @. (coeff == :linear_const) | (coeff == :linear_state)
    cm = _maximum_bipartite_matching(biadj)
    return [(cm[j], j) for j in 1:n_st if cm[j] > 0]
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
