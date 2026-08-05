####
#### Rule resolution: one graph for observed equations, output equations and init formulas
####
# Everything here is an algebraic rule `out = f(in…)`: observed equations, output equations and
# user formulas all look the same. They go into one bucket and the execution order falls out of
# which symbols are already known.
#
# This file is the executor. It knows nothing about components or metadata, and it decides
# nothing — which rules exist and what counts as an error is up to the caller.

"""
    ResolutionRule(f, outsym, sym; provenance, optional, label=nothing)

One rule `outsym = f(sym…)` of the resolution graph. The symbol names are canonical and the
payload is a plain `f(out, u)` on vectors, which is the shape `build_function` gives us anyway.

- `provenance` — the rank the values it writes get, see [`_precedence`](@ref).
- `optional` — whether it is fine for this rule to never fire.
- `label` — short name for diagnostics.

`optional` and a low `provenance` are two different escape hatches: the first lets a rule not run
at all, the second lets it run and lose.

A rule fires as a unit once all its inputs are known, but its outputs are written one by one, so
a collision on one of them leaves the others alone.
"""
struct ResolutionRule
    f::Function # potentialy FunctionWrapper
    outsym::Vector{Symbol}
    sym::Vector{Symbol}
    provenance::Symbol
    optional::Bool
    label::Union{Nothing,String}

    function ResolutionRule(f, outsym, sym; provenance, optional, label=nothing)
        self_deps = intersect(sym, outsym)
        if !isempty(self_deps)
            throw(ArgumentError("ResolutionRule cannot depend on its own output: $self_deps"))
        end
        _precedence(provenance) # reject an unknown provenance here rather than mid-walk
        new(f, collect(Symbol, outsym), collect(Symbol, sym), provenance, optional, label)
    end
end

"""
    _precedence(provenance::Symbol) -> Int

Precedence of a value's provenance, higher wins:

    :strong_formula > :provided > :guess_formula > :guess > :derived > :weak_formula

Formulas rank above the values they are meant to refine. Values a rule computed rank below
anything the user wrote down, guesses included, but above weak formulas, which only exist to
step aside for something better.
"""
function _precedence(provenance::Symbol)
    provenance === :strong_formula ? 6 :
    provenance === :provided       ? 5 :
    provenance === :guess_formula  ? 4 :
    provenance === :guess          ? 3 :
    provenance === :derived        ? 2 :
    provenance === :weak_formula   ? 1 :
    throw(ArgumentError("Unknown value provenance :$provenance"))
end

# prefer the label, `rule for [:out…]` reads badly mid-sentence
_rule_ref(r::ResolutionRule) = isnothing(r.label) ? "rule for $(r.outsym)" : "`$(r.label)`"

# obs/output rules are exactly the `:derived` ones. Not `!optional`, which only coincides.
_is_user_rule(r::ResolutionRule) = r.provenance !== :derived

"""
    ResolutionRule(f::InitFormula, am::AliasMap)
    ResolutionRule(f::GuessFormula, am::AliasMap)

Wrap a user formula as a rule. The symbol lists are canonicalized through `am`, while the user's
closure keeps being called under the names they originally wrote.
"""
ResolutionRule(f::InitFormula, am::AliasMap) =
    _wrap_formula(f, am; provenance=f.weak ? :weak_formula : :strong_formula)
ResolutionRule(f::GuessFormula, am::AliasMap) =
    _wrap_formula(f, am; provenance=:guess_formula)

function _wrap_formula(f, am; provenance)
    outsym = _canonical_names(f.outsym, am, f)
    sym = _canonical_names(f.sym, am, f)
    origout, origin = f.outsym, f.sym
    payload = function (out::Vector{Float64}, u::Vector{Float64})
        f.f(SymbolicView(out, origout), SymbolicView(u, origin))
        nothing
    end
    ResolutionRule(payload, outsym, sym; provenance, optional=false, label=f.label)
end

# two members of one alias class are the same variable, so listing both would write it twice
function _canonical_names(syms, am, f)
    names = [canonicalize(am, s) for s in syms]
    if !allunique(names)
        dupes = unique(n for n in names if count(isequal(n), names) > 1)
        throw(ArgumentError("$(typeof(f).name.name) uses $syms, which collapse onto $dupes — \
                             members of one alias class are one variable."))
    end
    names
end

"""
    ResolutionResult

What [`resolve_rules`](@ref) did. It only records things; erroring or warning is up to the caller.

- `vals` — the resolved values, usable as a `defaults`/`guesses` dict.
- `provenance` — where each value came from.
- `writer` — which rule wrote each value. Diagnostics need this, since several rules can share
  one provenance.
- `fired` / `pruned` — per rule.
- `unfired` — rules that never ran because their inputs stayed unknown.
- `nonfinite` — values that came out `NaN`/`Inf`, and who wrote them.
- `yields` — a rule lost a symbol to a higher-ranked one.
- `conflicts` — two rules of the same rank disagreed about a symbol.
"""
struct ResolutionResult
    vals::Dict{Symbol,Float64}
    provenance::Dict{Symbol,Symbol}
    writer::Dict{Symbol,Int}
    fired::BitVector
    pruned::BitVector
    unfired::Vector{@NamedTuple{rule::Int, unknown::Vector{Symbol}}}
    nonfinite::Vector{@NamedTuple{sym::Symbol, rule::Int}}
    conflicts::Vector{@NamedTuple{sym::Symbol, held::Float64, offered::Float64, rule::Int}}
    yields::Vector{@NamedTuple{sym::Symbol, offered::Float64, rule::Int}}
end

function ResolutionResult(vals, provenance, nrules::Int)
    ResolutionResult(
        vals,
        provenance,
        Dict{Symbol,Int}(),
        falses(nrules),
        falses(nrules),
        @NamedTuple{rule::Int, unknown::Vector{Symbol}}[],
        @NamedTuple{sym::Symbol, rule::Int}[],
        @NamedTuple{sym::Symbol, held::Float64, offered::Float64, rule::Int}[],
        @NamedTuple{sym::Symbol, offered::Float64, rule::Int}[]
    )
end

"""
    resolve_rules(vals, rules; targets, seedprov, maxpasses) -> ResolutionResult

Resolve everything that follows from `vals` by firing rules whose inputs are known. `vals` itself
is not modified.

The walk:

1. Seed the known values, each at the rank `seedprov` gives it.
2. Drop the rules that cannot contribute to `targets`.
3. Order the rules with Tarjan and walk them, re-scanning cyclic groups until they settle.
4. Write each output, unless something of higher rank is already sitting there. Equal rank is a
   consistency check rather than a write.

Conflicts are only recorded; the caller decides whether they are an error.
"""
function resolve_rules(vals, rules::AbstractVector;
                       targets=nothing, seedprov=Returns(:provided), maxpasses=1000)
    values = Dict{Symbol,Float64}(vals)
    provenance = Dict{Symbol,Symbol}(s => seedprov(s) for s in keys(values))

    res = ResolutionResult(values, provenance, length(rules))

    g = rule_graph(rules)
    keep = _prune_rules(g, rules, targets)
    res.pruned .= .!keep

    for block in _rule_blocks(g, keep)
        if length(block) == 1
            # there are no self-loops, so a lone rule gets exactly one chance
            _try_fire!(res, rules, only(block))
        else
            _resolve_block!(res, rules, block; maxpasses)
        end
    end

    for i in eachindex(rules)
        (res.fired[i] || res.pruned[i] || rules[i].optional) && continue
        push!(res.unfired, (; rule=i, unknown=[s for s in rules[i].sym if !haskey(res.vals, s)]))
    end

    # `writer` already knows who produced what, so one scan finds every non-finite value; the
    # `any` guard keeps the collect+sort off the hot path, the sort keeps the order stable.
    if any(!isfinite, Base.values(res.vals))
        for s in sort!(collect(keys(res.vals)))
            isfinite(res.vals[s]) ||
                push!(res.nonfinite, (; sym=s, rule=get(res.writer, s, 0)))
        end
    end
    res
end

"""
    rule_graph(rules) -> SimpleDiGraph

Dependency graph over `rules`: an edge `i → j` means rule `i` writes something rule `j` reads.

Vertices are rules, not variables, so a formula with several outputs stays a single node.

Edges mean "run before", not "requires": if two rules could write the same symbol, both are
ordered ahead of everything that reads it.
"""
function rule_graph(rules)
    writers = Dict{Symbol,Vector{Int}}()
    for (i, r) in enumerate(rules), s in r.outsym
        push!(get!(writers, s, Int[]), i)
    end

    g = SimpleDiGraph(length(rules))
    for (j, r) in enumerate(rules), s in r.sym
        haskey(writers, s) || continue
        for i in writers[s]
            add_edge!(g, i, j) # a self-edge is impossible, a rule cannot read its own output
        end
    end
    g
end

# walk backwards from the rules writing a target and keep everything reachable.
# `targets === nothing` keeps all rules.
function _prune_rules(g, rules, targets)
    isnothing(targets) && return trues(length(rules))

    targetset = Set{Symbol}(targets)
    seeds = [i for (i, r) in enumerate(rules) if !isdisjoint(r.outsym, targetset)]
    keep = falses(length(rules))
    for v in Graphs.BFSIterator(g, seeds; neighbors_type=Graphs.inneighbors)
        keep[v::Int] = true # the iterator's eltype is `Any`
    end
    keep
end

# Groups of mutually dependent rules, in execution order. Pruning has to happen before Tarjan
# runs, otherwise the groups come out bigger than they really are.
function _rule_blocks(g, keep)
    sg, vmap = Graphs.induced_subgraph(g, findall(keep))
    sccs = reverse(Graphs.strongly_connected_components_tarjan(sg))
    [sort!(vmap[scc]) for scc in sccs]
end

# Rules in a cycle cannot be put in order, so the group is re-scanned until nothing changes and
# whichever rule happens to have its inputs ready goes first. An overwrite re-arms everyone who
# read the old value.
#
# This settles because a rule fires at most once and an overwrite always raises the target's
# rank. `maxpasses` is only a safety net.
function _resolve_block!(res, rules, block; maxpasses)
    overwritten = Set{Symbol}()
    for _ in 1:maxpasses
        changed = false
        for i in block
            res.fired[i] && continue
            _try_fire!(res, rules, i; overwritten) && (changed = true)
        end
        changed || return nothing
        if !isempty(overwritten)
            for i in block
                res.fired[i] || continue
                isdisjoint(rules[i].sym, overwritten) && continue
                res.fired[i] = false # its inputs moved, its result is stale
            end
            empty!(overwritten)
        end
    end
    error("Resolution of a rule block did not settle within $maxpasses passes: the values of \
           $(unique(Symbol[s for i in block for s in rules[i].outsym])) keep changing. This \
           should not happen — please report it, together with the model that produced it.")
end

# Returns whether the rule fired. Firing depends on readiness (all inputs set?)
function _try_fire!(res, rules, i; overwritten=nothing)
    r = rules[i]
    all(s -> haskey(res.vals, s), r.sym) || return false

    u = Float64[res.vals[s] for s in r.sym]
    out = zeros(Float64, length(r.outsym))
    r.f(out, u)

    res.fired[i] = true
    for (k, s) in enumerate(r.outsym)
        _write_value!(res, s, out[k], r.provenance, i; overwritten)
    end
    true
end

# The write policy, and the only place where precedence is decided.
function _write_value!(res, s, v, provenance, i; overwritten)
    if !haskey(res.vals, s)
        _store!(res, s, v, provenance, i; overwritten=nothing)
        return nothing
    end
    held = _precedence(res.provenance[s])
    offered = _precedence(provenance)
    if held > offered
        # recorded, not silent: for a weak formula this *is* the outcome to report
        push!(res.yields, (; sym=s, offered=v, rule=i))
        return nothing
    elseif held == offered
        # same rank never overwrites, we only check that the two agree
        _agree(res.vals[s], v) ||
            push!(res.conflicts, (; sym=s, held=res.vals[s], offered=v, rule=i))
        return nothing
    end
    # the rule that got here first is the one that lost
    haskey(res.writer, s) && push!(res.yields, (; sym=s, offered=res.vals[s], rule=res.writer[s]))
    _store!(res, s, v, provenance, i; overwritten)
    nothing
end

# `overwritten` collects the symbols whose value actually moved, so a block knows what to
# re-read. Compared with a tolerance, otherwise jitter in the last bits would never settle.
function _store!(res, s, v, provenance, i; overwritten)
    if !isnothing(overwritten) && haskey(res.vals, s) && !_agree(res.vals[s], v)
        push!(overwritten, s)
    end
    res.vals[s] = v
    res.provenance[s] = provenance
    res.writer[s] = i
    nothing
end
