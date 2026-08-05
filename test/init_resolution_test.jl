using NetworkDynamics
using NetworkDynamics: ResolutionRule, resolve_rules, rule_graph, AliasMap, _precedence
using Graphs: Graphs
using Test

# The resolution executor is a pure function of a rule vector and a value dict: no component,
# no Symbolics, no ModelingToolkit. Everything here is hand-built, which is the point — the
# semantics that are hard to see through a compiled bus (readiness inside a block, precedence,
# termination) are one-liners at this level.

# rules take plain vectors, so a test rule is just a closure over indices
rule(f, outsym, sym; provenance=:derived, optional=true, label=nothing) =
    ResolutionRule(f, outsym, sym; provenance, optional, label)
# `out = factor * u[1]`, the shape both alias directions and most test rules take
scale(factor, out, in; kw...) = rule((o, u) -> (o[1] = factor * u[1]; nothing), [out], [in]; kw...)

@testset "rule construction" begin
    @test_throws ArgumentError rule((o, u) -> nothing, [:x], [:x])
    @test_throws ArgumentError ResolutionRule((o, u) -> nothing, [:x], [:y];
                                              provenance=:nonsense, optional=true)

    # tiers are symbols, but they are ordered
    @test _precedence(:strong_formula) > _precedence(:provided) > _precedence(:guess) >
          _precedence(:derived) > _precedence(:weak_formula)

    # concrete, so a Vector{ResolutionRule} is not a Vector{Any} in disguise
    @test isconcretetype(ResolutionRule)
end

@testset "formula wrapping" begin
    am = AliasMap(:alias_out => :y, :alias_in => :x)
    f = @initformula :alias_out = 2 * :alias_in
    r = ResolutionRule(f, am)

    # the rule speaks canonical names, the user's closure keeps speaking theirs
    @test r.outsym == [:y]
    @test r.sym == [:x]
    @test r.provenance == :strong_formula
    @test !r.optional
    out = zeros(1)
    r.f(out, [3.0])
    @test out == [6.0]

    @test ResolutionRule(@initformula(weak = true, :alias_out = :alias_in), am).provenance == :weak_formula
    @test ResolutionRule(@guessformula(:alias_out = :alias_in), am).provenance == :guess

    # two members of one alias class in one list are one variable, not two — on either end
    am2 = AliasMap(:a => :x, :b => :x)
    collapsing = @initformula begin
        :a = 1.0
        :b = 2.0
    end
    @test_throws ArgumentError ResolutionRule(collapsing, am2)
    @test_throws ArgumentError ResolutionRule(@initformula(:out = :a + :b), am2)

    # and the case the canonicalization *creates*: two names that look independent in the
    # formula but land on one symbol across the two ends, which is a rule reading its own output
    @test_throws ArgumentError ResolutionRule(@initformula(:x = 2 * :alias_in), am)
    @test_throws ArgumentError ResolutionRule(@initformula(:alias_in = 2 * :x), am)
end

@testset "rule graph" begin
    rules = [scale(2.0, :y, :x),                                    # 1: x → y
             scale(3.0, :z, :y),                                    # 2: y → z
             rule((o, u) -> (o[1] = u[1]; o[2] = u[1]; nothing), [:a, :b], [:z])]  # 3: z → a,b
    g = rule_graph(rules)
    @test Graphs.nv(g) == 3
    @test Set(Tuple.(Graphs.edges(g))) == Set([(1, 2), (2, 3)])

    # a multi-output rule is *one* vertex — that is what spares it any special handling
    both = [rule((o, u) -> (o[1] = u[1]; o[2] = u[1]; nothing), [:p, :q], [:seed]),
            scale(1.0, :r, :p), scale(1.0, :s, :q)]
    gb = rule_graph(both)
    @test Graphs.nv(gb) == 3
    @test Set(Tuple.(Graphs.edges(gb))) == Set([(1, 2), (1, 3)])

    # two rules writing the same symbol are independent — no edge between them
    @test Graphs.ne(rule_graph([scale(2.0, :q, :p), scale(3.0, :q, :p)])) == 0
end

@testset "non-finite seeds are ordinary values" begin
    # "unknown" is the absence of a key, so nothing is filtered: `NaN` and `Inf` are seeded and
    # carry `:provided` like any number. `nothing`/`missing` are turned away one level up, in
    # `_numeric_seed`, so `resolve_rules` never has to test for them.
    res = resolve_rules(Dict(:a => 1.0, :b => NaN, :e => Inf), ResolutionRule[])
    @test keys(res.vals) == Set([:a, :b, :e])
    @test res.provenance[:b] == :provided
    @test isnan(res.vals[:b])
    @test res.vals[:e] == Inf
    # the seeds themselves show up in `nonfinite`, credited to no rule
    @test res.nonfinite == [(; sym=:b, rule=0), (; sym=:e, rule=0)]

    # a rule mapping the infinite input to a finite limit contributes its value; one that stays
    # infinite contributes that, and is named together with the rule that produced it
    limit = [rule((o, u) -> (o[1] = exp(-u[1]); nothing), [:decayed], [:e]),
             scale(2.0, :ramp, :e)]
    res2 = resolve_rules(Dict(:e => Inf), limit; targets=[:decayed, :ramp])
    @test res2.vals[:decayed] == 0.0
    @test res2.vals[:ramp] == Inf
    @test res2.nonfinite == [(; sym=:e, rule=0), (; sym=:ramp, rule=2)]
end

@testset "orientation is a property of the query" begin
    # The motivating case: `i_r ~ -term_i_r` is the output equation (invertible), and
    # `term_i_r ~ I*cos(δ)` the observed equation defining the same observable. A formula
    # computing δ *from* term_i_r closes a cycle with the observed equation — the block is
    # mixed, and readiness is what resolves it.
    rules = [scale(-1.0, :term_i_r, :i_r; label="inverse output eq"),
             scale(-1.0, :i_r, :term_i_r; label="output eq"),
             rule((o, u) -> (o[1] = u[1] * cos(u[2]); nothing), [:term_i_r], [:I, :δ]; label="obs"),
             rule((o, u) -> (o[1] = acos(u[1] / u[2]); nothing), [:δ], [:term_i_r, :I];
                  provenance=:strong_formula, optional=false, label="δ formula")]

    res = resolve_rules(Dict(:i_r => -0.5, :I => 1.0), rules; targets=[:δ, :i_r, :I])
    @test res.vals[:term_i_r] ≈ 0.5      # from the inverted output equation
    @test res.vals[:δ] ≈ acos(0.5)       # from the formula, which could not have expanded
    @test all(res.fired)
    # the observed equation became ready only after δ was known, and agreed — a free check
    @test isempty(res.conflicts)
    @test isempty(res.unfired)

    # the mirror query, which the old compile-time orientation could not serve at the same
    # time: hand it term_i_r instead and the output equation runs forward
    mirror = resolve_rules(Dict(:term_i_r => 0.5, :I => 1.0), rules; targets=[:δ, :i_r, :I])
    @test mirror.vals[:i_r] ≈ -0.5
    @test mirror.vals[:δ] ≈ acos(0.5)
    @test isempty(mirror.conflicts)
end

@testset "pruning" begin
    obs = [scale(2.0, :y, :x), scale(3.0, :z, :y)]

    # nothing reaches a settable target, so the whole graph goes and nothing is computed —
    # this is what keeps a component without formulas free of cost
    res = resolve_rules(Dict(:x => 1.0), obs; targets=[:x])
    @test all(res.pruned)
    @test !haskey(res.vals, :y)

    # a target downstream keeps exactly the chain that reaches it
    res2 = resolve_rules(Dict(:x => 1.0), obs; targets=[:z])
    @test !any(res2.pruned)
    @test res2.vals[:z] ≈ 6.0

    # no targets at all means "resolve everything"
    res3 = resolve_rules(Dict(:x => 1.0), obs)
    @test res3.vals[:y] ≈ 2.0

    # a pruned rule never fires, so it is not reported as unfired either
    required = [scale(2.0, :dead, :x; provenance=:strong_formula, optional=false)]
    res4 = resolve_rules(Dict(:x => 1.0), required; targets=[:x])
    @test all(res4.pruned)
    @test isempty(res4.unfired)
end

@testset "precedence" begin
    rules = [scale(99.0, :a, :seed),
             scale(3.5, :b, :seed; provenance=:strong_formula, optional=false),
             scale(5.0, :c, :seed; provenance=:weak_formula, optional=false)]
    seeded = Dict(:seed => 2.0, :a => 1.0, :b => 1.0, :c => 1.0)
    res = resolve_rules(seeded, rules; targets=[:a, :b, :c])

    @test res.vals[:a] == 1.0 && res.provenance[:a] == :provided        # derived yields
    @test res.vals[:b] == 7.0 && res.provenance[:b] == :strong_formula  # strong overwrites
    @test res.vals[:c] == 1.0 && res.provenance[:c] == :provided        # weak yields
    @test isempty(res.conflicts)

    # with the targets free instead, every rule lands on an empty slot
    free = resolve_rules(Dict(:seed => 2.0), rules; targets=[:a, :b, :c])
    @test free.vals[:a] == 198.0
    @test free.vals[:c] == 10.0
    @test free.provenance[:c] == :weak_formula

    # a derived value still beats a weak formula: the weak one yields to anything that
    # actually determines its target
    over = [scale(2.0, :w, :seed), scale(9.0, :w, :seed; provenance=:weak_formula, optional=false)]
    @test resolve_rules(Dict(:seed => 2.0), over; targets=[:w]).vals[:w] == 4.0
end

@testset "a superseded rule is recorded, whichever way round it fired" begin
    # The same fact reaches the executor two ways depending on the walk: the weak rule lands on a
    # value that already outranks it (a yield), or it got there first and is overwritten. Both are
    # recorded under the *losing* rule, so a caller's log does not depend on the firing order.
    weak_last = [scale(2.0, :w, :seed; provenance=:strong_formula, optional=false),
                 scale(9.0, :w, :seed; provenance=:weak_formula, optional=false)]
    a = resolve_rules(Dict(:seed => 2.0), weak_last; targets=[:w])
    @test a.vals[:w] == 4.0
    @test only(a.yields) == (; sym=:w, offered=18.0, rule=2)

    weak_first = reverse(weak_last)
    b = resolve_rules(Dict(:seed => 2.0), weak_first; targets=[:w])
    @test b.vals[:w] == 4.0
    @test only(b.yields) == (; sym=:w, offered=18.0, rule=1)   # the weak rule, now index 1

    # yielding to a seed is the same thing with no rule on the other side
    seeded = resolve_rules(Dict(:seed => 2.0, :w => 1.0), [weak_last[2]]; targets=[:w])
    @test seeded.vals[:w] == 1.0
    @test only(seeded.yields) == (; sym=:w, offered=18.0, rule=1)

    # nothing is recorded when a rule simply wins
    @test isempty(resolve_rules(Dict(:seed => 2.0), [weak_last[1]]; targets=[:w]).yields)
end

@testset "collisions at equal precedence" begin
    dup = [scale(2.0, :q, :p), scale(3.0, :q, :p)]
    res = resolve_rules(Dict(:p => 1.0), dup; targets=[:q])

    # one of them assigns, the other lands on it and is taken as a check. *Which* is not
    # promised — the graph does not order two independent writers, and at equal provenance the
    # check has to agree within tolerance anyway. Only the disagreement is the news.
    @test res.vals[:q] == 3.0
    @test only(res.conflicts) == (; sym=:q, held=3.0, offered=2.0, rule=1)

    # agreeing writers are not a conflict — that is the redundancy the old structural
    # duplicate-writer check rejected outright
    agree = [scale(2.0, :q, :p), scale(2.0, :q, :p)]
    @test isempty(resolve_rules(Dict(:p => 1.0), agree; targets=[:q]).conflicts)
end

@testset "a non-finite output fires and propagates" begin
    rules = [rule((o, u) -> (o[1] = u[1] / 0 - u[1] / 0; nothing), [:t], [:x];
                  provenance=:strong_formula, optional=false, label="nan rule"),
             scale(2.0, :downstream, :t)]
    res = resolve_rules(Dict(:x => 1.0), rules; targets=[:t, :downstream])

    # the rule ran, so it fired and wrote — suppressing the write would only hide the origin
    @test res.fired == [true, true]
    @test isnan(res.vals[:t])
    @test isnan(res.vals[:downstream]) # the poison spreads, by design
    @test isempty(res.unfired)
    # both are reported, each credited to the rule that wrote it
    @test res.nonfinite == [(; sym=:downstream, rule=2), (; sym=:t, rule=1)]
end

@testset "a NaN inside a block does not stall the fixpoint" begin
    # `_agree(NaN, NaN)` has to hold, or `_store!` would see the value move on every pass and
    # `_resolve_block!` would burn `maxpasses` and error out
    cyc = [rule((o, u) -> (o[1] = u[1] * NaN; nothing), [:a], [:b];
                provenance=:strong_formula, optional=false),
           rule((o, u) -> (o[1] = u[1] + 1; nothing), [:b], [:a];
                provenance=:strong_formula, optional=false),
           scale(1.0, :b, :seed)]
    res = resolve_rules(Dict(:seed => 1.0), cyc; targets=[:a, :b])
    # `:b` is seeded 1.0 by the `:derived` rule, then overwritten NaN by the strong one that
    # reads `:a` — the poison closes the cycle and the block still settles
    @test isnan(res.vals[:a]) && isnan(res.vals[:b])
    @test [nf.sym for nf in res.nonfinite] == [:a, :b]
end

@testset "unfired rules name what they are missing" begin
    rules = [rule((o, u) -> (o[1] = u[1] + u[2]; nothing), [:out], [:known, :absent];
                  provenance=:strong_formula, optional=false)]
    res = resolve_rules(Dict(:known => 1.0), rules; targets=[:out])
    @test only(res.unfired) == (; rule=1, unknown=[:absent])

    # an optional rule that cannot fire is not an error, it is the normal case for observables
    optional = [rule((o, u) -> (o[1] = u[1]; nothing), [:o], [:absent])]
    @test isempty(resolve_rules(Dict{Symbol,Float64}(), optional; targets=[:o]).unfired)
end

@testset "multi-output formulas are one node" begin
    rules = [rule((o, u) -> (o[1] = u[1]; o[2] = 2 * u[1]; nothing), [:g, :h], [:k];
                  provenance=:strong_formula, optional=false),
             rule((o, u) -> (o[1] = u[1] + u[2]; nothing), [:s], [:g, :h])]
    res = resolve_rules(Dict(:k => 3.0), rules; targets=[:s])
    @test res.vals[:s] ≈ 9.0

    # a collision on one output must not strand its siblings: :g is already provided, :h is
    # not, and the formula still delivers :h
    partial = resolve_rules(Dict(:k => 3.0, :g => 100.0), rules; targets=[:s, :g, :h])
    @test partial.vals[:g] == 3.0     # strong displaces the provided 100.0
    @test partial.provenance[:g] == :strong_formula
    @test partial.vals[:h] == 6.0     # ... and the sibling output lands either way
end

@testset "blocks terminate" begin
    # same tier both ways, seeded above them: the closing rule yields, nothing loops
    yielding = [scale(2.0, :m, :n), scale(0.6, :n, :m)]
    res = resolve_rules(Dict(:n => 1.0), yielding; targets=[:m, :n])
    @test res.vals == Dict(:m => 2.0, :n => 1.0)
    @test isempty(res.conflicts)

    # both strong, so the second one *does* raise :n above the seed, which re-arms the first
    # rule and repairs its stale result. The repaired value disagrees, and that inconsistency
    # is reported rather than iterated forever.
    raising = [scale(2.0, :m, :n; provenance=:strong_formula, optional=false),
               scale(0.6, :n, :m; provenance=:strong_formula, optional=false)]
    res2 = resolve_rules(Dict(:n => 1.0), raising; targets=[:m, :n])
    @test res2.vals[:n] ≈ 1.2
    @test only(res2.conflicts).sym == :m

    # a consistent cycle closes cleanly: 0.5 * 2 == 1
    consistent = [scale(2.0, :m, :n; provenance=:strong_formula, optional=false),
                  scale(0.5, :n, :m; provenance=:strong_formula, optional=false)]
    res3 = resolve_rules(Dict(:n => 1.0), consistent; targets=[:m, :n])
    @test res3.vals[:m] ≈ 2.0
    @test isempty(res3.conflicts)
end

@testset "the plan does not depend on hash order" begin
    # Same rules, same values, eight independent blocks. Which of them assigns is decided by
    # the graph walk and not promised here — but it must be *stable*, i.e. a function of the
    # rule vector alone. Dict iteration must never reach the plan.
    mk() = [scale(Float64(i), :dup, :seed) for i in 1:8]
    res = resolve_rules(Dict(:seed => 1.0), mk(); targets=[:dup])
    @test length(res.conflicts) == 7 # one assignment, the rest land on it as checks

    for _ in 1:20
        again = resolve_rules(Dict(:seed => 1.0), mk(); targets=[:dup])
        @test again.vals[:dup] == res.vals[:dup]
        @test [c.rule for c in again.conflicts] == [c.rule for c in res.conflicts]
    end
end
