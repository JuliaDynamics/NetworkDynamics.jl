using NetworkDynamics
using NetworkDynamics: ResolutionRule, resolve_rules, has_obsrules, get_obsrules,
                       get_aliasmap, settable_symbols, obssym, sym, psym,
                       _get_component_observed
using ModelingToolkitBase
using ModelingToolkitBase: t_nounits as t, D_nounits as Dt
using SciCompDSL
using Test

# The `:derived` half of the resolution graph, built from a component's symbolic equations at
# compile time. What is tested here is the *translation* — which equations become which rules —
# plus the one invariant that anchors it: firing every rule from a full root set must reproduce
# `obsf` exactly. The executor itself is tested in `init_resolution_test.jl`, without MTK.

# One bus carrying every shape the translation has to tell apart: an identity alias, a
# sign-flipped relation, a chain of them, a genuinely algebraic observable, a multi-root one,
# one that depends explicitly on time, and a constant one.
@mtkmodel ObsRuleBus begin
    @variables begin
        u_r(t) = 1.0
        u_i(t) = 0.0
        term_u_r(t); θ(t); scaled(t); nl(t); Vmeas(t); summed(t); timed(t); konst(t)
        i_r(t), [input=true]
        P(t), [output=true]
    end
    @parameters begin
        Vset = 1.0
    end
    @equations begin
        Dt(u_r) ~ -u_r + i_r
        Dt(u_i) ~ -u_i
        term_u_r ~ u_r              # identity alias: one variable, no rule
        θ ~ -u_r                    # sign flipped: two rules, one per direction
        scaled ~ 2*θ                # likewise
        nl ~ u_r^2                  # not invertible: forward only
        Vmeas ~ Vset                # a parameter defines it — still an invertible relation
        summed ~ u_r + u_i + Vset   # multi-root
        timed ~ u_r * t             # explicitly time dependent
        konst ~ 42.0                # no inputs
        P ~ u_r + nl                # the output equation
    end
end
@named _orb = ObsRuleBus()
VM = VertexModel(_orb, [:i_r], [:P])
RULES = get_obsrules(VM)

# rules as `out => in` pairs, which is all most assertions below care about
pairs_of(rules) = [(only(r.outsym) => r.sym) for r in rules]

@testset "which equations become which rules" begin
    @test has_obsrules(VM)
    @test RULES isa Vector{ResolutionRule}
    @test all(r -> r.provenance == :derived && r.optional, RULES)

    P = pairs_of(RULES)

    # an invertible relation is a *pair* of rules, and that pair is the whole point: which of
    # the two ends determines the other is decided at init time, not here
    @test (:θ => [:u_r]) in P && (:u_r => [:θ]) in P
    @test (:scaled => [:θ]) in P && (:θ => [:scaled]) in P
    @test (:Vmeas => [:Vset]) in P && (:Vset => [:Vmeas]) in P

    # ... while a non-invertible one is forward only
    @test (:nl => [:u_r]) in P
    @test !((:u_r => [:nl]) in P)

    # an identity alias is not a relation between two variables, it is one variable under two
    # names. Both directions would be `u_r = u_r`, so neither is a rule — and nothing else
    # mentions the alias, `pick_best_alias_names` substituted it away everywhere.
    @test get_aliasmap(VM) == NetworkDynamics.AliasMap(:term_u_r => :u_r)
    @test !any(p -> first(p) == :term_u_r || :term_u_r in last(p), P)

    # the independent variable is an ordinary input, but never an output: it is an argument
    # of the model, not an unknown of it
    @test (:timed => [:u_r, :t]) in P
    @test !any(p -> first(p) == :t, P)

    # a multi-root rule reads its variables in whatever order `get_variables` hands them back
    summedrule = only(filter(r -> r.outsym == [:summed], RULES))
    @test Set(summedrule.sym) == Set([:u_r, :u_i, :Vset])
    @test (:konst => Symbol[]) in P
    @test (:P => [:u_r, :nl]) in P

    # `:observed` stays exactly what `obsf` computes; the rules are a second view of the same
    # equations, never a rewrite of them
    @test length(NetworkDynamics.get_metadata(VM, :observed)) == length(obssym(VM))
end

@testset "firing everything reproduces obsf" begin
    vals = Dict(:u_r => 0.7, :u_i => -0.3, :Vset => 1.3, :i_r => 0.2, :t => 2.5)
    res = resolve_rules(vals, RULES; targets=nothing)

    @test isempty(res.conflicts)
    @test isempty(res.nonfinite)
    @test all(res.fired)

    # an observable sharing a slot with a settable symbol carries no entry of its own — the
    # invariant is "reproduces `obsf` *up to slot identity*"
    am = get_aliasmap(VM)
    reference = _get_component_observed(VM, vals; t=2.5)
    for (s, v) in zip(obssym(VM), reference)
        @test res.vals[NetworkDynamics.canonicalize(am, s)] ≈ v
    end
    @test res.vals[:P] ≈ 0.7 + 0.7^2

    # nothing derived may displace what the caller provided
    @test res.provenance[:u_r] == :provided
    @test res.provenance[:θ] == :derived
end

@testset "time is an input like any other" begin
    # without a `t` the time dependent observable simply does not resolve — it is an optional
    # rule, so that is a skip and not a failure
    res = resolve_rules(Dict(:u_r => 2.0), RULES; targets=nothing)
    @test !haskey(res.vals, :timed)
    @test isempty(res.unfired) # everything here is optional
    @test res.vals[:nl] ≈ 4.0
end

# The case the whole design exists for: on a component whose *output* is a scaled alias of an
# observable, which end determines the other depends on what the query writes. As a rule pair
# both queries are answerable from one rule set.
@mtkmodel InvBus begin
    @variables begin
        x(t) = 1.0
        y(t)
        u(t), [input=true]
        o(t), [output=true]
    end
    @equations begin
        Dt(x) ~ -x + u
        y ~ x^2
        o ~ -y
    end
end
@named _ib = InvBus()
IVM = VertexModel(_ib, [:u], [:o])

@testset "output equations are rules, in both directions" begin
    P = pairs_of(get_obsrules(IVM))
    @test (:y => [:x]) in P
    @test (:o => [:y]) in P
    @test (:y => [:o]) in P

    # forward: from the state, through the observable, onto the output
    fwd = resolve_rules(Dict(:x => 3.0), get_obsrules(IVM); targets=nothing)
    @test fwd.vals[:y] ≈ 9.0 && fwd.vals[:o] ≈ -9.0

    # mirror: from the output back onto the observable, which no forward-only rule set could do
    back = resolve_rules(Dict(:o => -4.0), get_obsrules(IVM); targets=nothing)
    @test back.vals[:y] ≈ 4.0
    @test !haskey(back.vals, :x) # `y ~ x^2` is not invertible and is not inverted
end

@testset "pruning to the targets" begin
    # a component with no formulas has nothing reading its observables, so with the `NWState`
    # target set — states and parameters only — the whole graph prunes away
    res = resolve_rules(Dict(:x => 3.0), get_obsrules(IVM);
                        targets=union(Set(sym(IVM)), Set(psym(IVM))))
    @test all(res.pruned)
    @test !haskey(res.vals, :y)

    # with the outputs among the targets, the output rule seeds the walk and drags its whole
    # upstream cone back in
    res = resolve_rules(Dict(:x => 3.0), get_obsrules(IVM); targets=settable_symbols(IVM))
    @test !any(res.pruned)
    @test res.vals[:o] ≈ -9.0
end

@testset "components not compiled from MTK carry no rules" begin
    cf = VertexModel(f=(du, u, in, p, t) -> du .= u, g=1:1, sym=[:x], outsym=[:o])
    @test !has_obsrules(cf)
    @test isempty(get_obsrules(cf))
end
