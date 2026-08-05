####
#### Rule resolution: one graph for observed equations, output equations and init formulas
####
# In short: every observed equation, every output
# equation and every user formula is the same thing — an algebraic rule `out = f(in…)`. Put
# them in one bucket, add inverse rules where a relation is invertible, and derive the
# execution order from which symbols are known. Orientation becomes a property of the query
# instead of something baked into the stored model.
#
# This file is the *executor* only: it knows nothing about components, MTK or metadata. What
# rules exist and what to do with the result is the caller's business — `initialize_component`
# errors where `NWState` reconstruction skips, and neither policy belongs here.

"""
    ResolutionRule(f, outsym, sym; provenance, optional, label=nothing)

One algebraic rule `outsym = f(sym…)` in the resolution graph.

The payload is plain vector-in, vector-out:

    f(out::Vector{Float64}, u::Vector{Float64})::Nothing

which is what a `build_function` result already looks like, so an observed equation wraps with
no adaptation at all. A user formula indexes by symbol instead, and the constructor taking an
[`InitFormula`](@ref) puts the [`SymbolicView`](@ref) wrapping *inside* the payload. The
executor therefore never sees symbolic indexing, and one loop drives both kinds.

`outsym`/`sym` are **canonical** names. The `InitFormula` constructor below canonicalizes them
and lets the payload keep speaking the names the user wrote — which is the whole of what the
old `normalize` pass did to a formula's two ends.

- `provenance` — what this rule stamps on the values it writes, see [`_precedence`](@ref).
- `optional` — whether failing to fire is acceptable. Observed and output rules are optional
  (nothing says their inputs must ever be known); user init formulas are not.
- `label` — short human reference for diagnostics, may be `nothing`.

A rule can be droppable from either end, and the two are independent — all four combinations are
meaningful:

- **weak** (a `provenance` below `:provided`) drops it on the *target* side: the value is already
  there, so the rule yields.
- **optional** drops it on the *source* side: the inputs never arrived, so the rule is skipped.

A rule fires *atomically*: it runs once all its inputs are known and produces all its outputs,
but each output is then written independently (see [`resolve_rules`](@ref)). That is what makes
a multi-output formula a plain single node of the graph rather than a special case — a
collision on one output cannot strand its siblings.
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

    :strong_formula > :provided > :derived > :weak_formula

`:derived` is what a rule of the graph produces from other values. It sits *below* `:provided`
because an observed or output equation is definitional — it computes a value, it does not
assert one, so it must never displace something the user wrote down. It sits *above*
`:weak_formula` because a weak formula's purpose is to yield to anything that actually
determines its target, and a derived value determines it.
"""
function _precedence(provenance::Symbol)
    provenance === :strong_formula ? 4 :
    provenance === :provided       ? 3 :
    provenance === :derived        ? 2 :
    provenance === :weak_formula   ? 1 :
    throw(ArgumentError("Unknown value provenance :$provenance"))
end

# `for [:out…]` reads badly mid-sentence, so name a rule by its label where it has one —
# same convention as `_formula_ref` for the formula types.
_rule_ref(r::ResolutionRule) = isnothing(r.label) ? "rule for $(r.outsym)" : "`$(r.label)`"

"""
    ResolutionRule(f::InitFormula, am::AliasMap)
    ResolutionRule(f::GuessFormula, am::AliasMap)

Wrap a user formula as a rule: canonicalize both symbol lists through `am`, and keep calling
the user's own closure under the names they wrote.

Canonicalizing renames, same order and same length, so the payload needs no arithmetic and no
permutation — it just presents the buffers under the original names. A scaled relation is not
an alias and never reaches here; it is an invertible pair of rules in the graph.

A `GuessFormula` writes at the lowest precedence: a guess is a seed, not an assertion, so it yields
to anything the init pass determined. The rest of the guess pass's collision policy is
`resolve_rules`' `overwrite_equal`.
"""
ResolutionRule(f::InitFormula, am::AliasMap) =
    _wrap_formula(f, am; provenance=f.weak ? :weak_formula : :strong_formula)
ResolutionRule(f::GuessFormula, am::AliasMap) =
    _wrap_formula(f, am; provenance=:weak_formula)

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

# Members of one alias class are one variable, so two of them in a rule's output list would
# mean writing the same variable twice — caught here rather than as a mystery collision.
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

What [`resolve_rules`](@ref) produced, and what it observed on the way. Deliberately free of
policy: it records what happened and leaves erroring, warning or ignoring to the caller.

- `vals` — the resolved values, a fresh `Dict{Symbol,Float64}`. Directly usable as the
  `defaults`/`guesses` dict downstream.
- `provenance` — where every entry of `vals` came from.
- `writer` — which rule wrote each entry, as an index into the rule vector; a seeded entry that
  no rule touched has no key here. Provenance alone does not identify the writer (two formulas
  share one), so this is what a diagnostic must key on to credit a value to a rule.
- `fired` / `pruned` — per rule, indexed like the rule vector handed to `resolve_rules`.
- `unfired` — non-optional rules that never ran, with the inputs still unknown. This is the
  "to compute `:δ`, provide one of {…}" material. A rule that *ran* and had all of its writes
  yield counts as `fired`: it was tested against what was already there, not dropped.
- `nonfinite` — the symbols holding a non-finite value, with the rule that wrote each (`0` for a
  seed nothing touched). Purely a diagnostic: such a value is written and propagated like any
  other, and this only says where it entered.
- `conflicts` — two rules of equal precedence landing on one target with disagreeing values.
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
        @NamedTuple{sym::Symbol, held::Float64, offered::Float64, rule::Int}[]
    )
end

"""
    resolve_rules(vals, rules; targets, seed_provenance, overwrite_equal, maxpasses) -> ResolutionResult

Resolve everything derivable from `vals` by firing rules whose inputs are known, in an order
derived from the rules themselves. `vals` is the root set and is **not** mutated: the result
carries a fresh `Dict{Symbol,Float64}`.

The walk:

1. Seed, with provenance `seed_provenance` (`:provided` for the init pass).
2. Prune to what can still matter, if `targets` is given: reverse reachability from the targets
   over the rules' outputs. Callers should pass targets rather than rely on the walk being cheap.
3. Tarjan over the *rule* graph (a rule is one node however many outputs it has), then walk the
   condensation in topological order. Singleton blocks fire at most once; multi-member blocks
   are re-scanned until nothing changes, with readiness deciding which member goes first.
4. Write policy per target, comparing the rule's provenance `rp` against the target's current
   `tp`: empty slot → write; `tp > rp` → yield silently; `tp == rp` → tolerance check, recorded
   as a conflict on disagreement; `tp < rp` → overwrite.

`overwrite_equal=true` relaxes the equal-precedence case to an overwrite without a check, which
is what the guess pass wants: refining a guess is desirable, not a contradiction.
"""
function resolve_rules(vals, rules::AbstractVector;
                       targets=nothing, seed_provenance=:provided,
                       overwrite_equal=false, maxpasses=1000)
    _precedence(seed_provenance)
    values = Dict{Symbol,Float64}(vals)
    provenance = Dict{Symbol,Symbol}(s => seed_provenance for s in keys(values))

    res = ResolutionResult(values, provenance, length(rules))

    g = rule_graph(rules)
    keep = _prune_rules(g, rules, targets)
    res.pruned .= .!keep

    for block in _rule_blocks(g, keep)
        if length(block) == 1
            # a rule cannot read its own output, so a singleton block has no self-loop and
            # gets exactly one chance
            _try_fire!(res, rules, only(block); overwrite_equal)
        else
            _resolve_block!(res, rules, block; overwrite_equal, maxpasses)
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

Dependency graph over `rules`: vertex `i` is `rules[i]`, and there is an edge `i → j` whenever
an output of rule `i` is an input of rule `j`, i.e. "`i` must run before `j`".

**Nodes are rules, not variables.** A rule with several outputs is one vertex, which is what
spares multi-output formulas any special handling downstream — Tarjan then treats such a
formula atomically without any merge-then-split bookkeeping.

An edge means "run before", **not** "requires": two rules writing one symbol both get an edge to
every reader, though either alone would determine it. That is deliberate — it puts all candidate
writers ahead of any reader, so precedence settles before the value is read. Need is checked at
firing time instead, per symbol. The cost is that a redundant writer sitting *downstream* of a
reader lands both in one SCC, i.e. redundancy can enlarge blocks.

Edges are fully determined by the order of `rules`: no dict iteration reaches the graph, so the
execution plan does not shift with Julia's hash order.
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

# Pruning is plain **reverse reachability in the rule graph**: seed with the rules writing a
# target, then walk backwards, because a predecessor is by construction a rule writing something
# a kept rule reads. `neighbors_type=inneighbors` is what walks the edges the wrong way round,
# so no reversed copy of the graph is materialized.
#
# `targets === nothing` means "keep everything" — for callers that genuinely want every
# observable resolved, and for tests.
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

# The kept rules as blocks of mutually dependent ones, in execution order, as indices into
# `rules`. Restricting to `keep` first matters: dropping a vertex can split an SCC, so running
# Tarjan on the full graph and filtering afterwards would leave blocks coarser than they are.
function _rule_blocks(g, keep)
    sg, vmap = Graphs.induced_subgraph(g, findall(keep))
    sccs = reverse(Graphs.strongly_connected_components_tarjan(sg))
    [sort!(vmap[scc]) for scc in sccs]
end

# Inside a block a writer of `x` may sit downstream of a reader of `x`, so one pass is not
# enough: re-scan until nothing changes. Readiness is what selects — the rule whose inputs
# happen to be known goes first, and the order of `rules` only decides between two rules ready
# at the same moment.
#
# Termination: a rule fires at most once, so a block settles in at most `|block|` firings —
# except that an *overwrite* re-arms the rules reading the overwritten symbol, so their stale
# results get repaired. That cannot run away either, because an overwrite strictly raises the
# target's provenance and the precedence levels are bounded. `maxpasses` guards a hole in that
# argument, it is not part of the design.
function _resolve_block!(res, rules, block; overwrite_equal, maxpasses)
    overwritten = Set{Symbol}()
    for _ in 1:maxpasses
        changed = false
        for i in block
            res.fired[i] && continue
            _try_fire!(res, rules, i; overwrite_equal, overwritten) && (changed = true)
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
function _try_fire!(res, rules, i; overwrite_equal, overwritten=nothing)
    r = rules[i]
    all(s -> haskey(res.vals, s), r.sym) || return false

    u = Float64[res.vals[s] for s in r.sym]
    out = zeros(Float64, length(r.outsym))
    r.f(out, u)

    res.fired[i] = true
    for (k, s) in enumerate(r.outsym)
        _write_value!(res, s, out[k], r.provenance, i; overwrite_equal, overwritten)
    end
    true
end

# The write policy, and the one place precedence is decided. Purely local to the target: the
# topological order has already put every writer of `s` ahead of every reader of it, so an
# overwrite here can never strand a value someone computed from the old one — outside a block,
# where `_resolve_block!` re-arms the readers instead.
function _write_value!(res, s, v, provenance, i; overwrite_equal, overwritten)
    if !haskey(res.vals, s)
        _store!(res, s, v, provenance, i; overwritten=nothing)
        return nothing
    end
    held = _precedence(res.provenance[s])
    offered = _precedence(provenance)
    if held > offered
        return nothing # yield silently: something stronger already determined this
    elseif held == offered && !overwrite_equal
        # equal precedence is a check, never a write — that is what bounds the block fixpoint, and it
        # is the free half of the consistency checking that pruning otherwise costs us
        _agree(res.vals[s], v) ||
            push!(res.conflicts, (; sym=s, held=res.vals[s], offered=v, rule=i))
        return nothing
    end
    _store!(res, s, v, provenance, i; overwritten)
    nothing
end

# `overwritten` collects what a block has to re-read, so the test is whether the *value* moved,
# not whether the precedence did — an `overwrite_equal` write strands every result computed
# from the old number just as a precedence-raising one does. Tolerance rather than `!=`
# deliberately: a value jittering in the last bits would otherwise re-arm its readers forever.
function _store!(res, s, v, provenance, i; overwritten)
    if !isnothing(overwritten) && haskey(res.vals, s) && !_agree(res.vals[s], v)
        push!(overwritten, s)
    end
    res.vals[s] = v
    res.provenance[s] = provenance
    res.writer[s] = i
    nothing
end
