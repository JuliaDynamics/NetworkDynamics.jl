using NetworkDynamics
using NetworkDynamics: psym, obssym, resolve_default_from,
                       set_default_from!, get_default_from, has_default_from,
                       SymbolicView, set_default!
using Graphs
using ModelingToolkitBase
using ModelingToolkitBase: t_nounits as t, D_nounits as D, System, @variables, @parameters, @named
using Test

@__MODULE__() == Main ? includet(joinpath(pkgdir(NetworkDynamics), "test", "ComponentLibrary.jl")) : (const Lib = Main.Lib)

# ---------------------------------------------------------------------------------------------
# Plain (non-MTK) routing models. `S_b` is a per-unit base parameter — inert in these trivial
# dynamics, it only exists to be copied around by `default_from`. CHECK_COMPONENT is toggled off
# while building them because an intentionally-unused base param would otherwise be flagged.
# ---------------------------------------------------------------------------------------------
hubf(du, u, agg, p, t) = (du[1] = -u[1] + agg[1]; nothing)
linef(e, vs, vd, p, t) = (e[1] = p[1] * (vs[1] - vd[1]); nothing)

function _models()
    old = NetworkDynamics.CHECK_COMPONENT[]
    NetworkDynamics.CHECK_COMPONENT[] = false
    try
        mkhub(sb, nm) = VertexModel(f=hubf, sym=[:u=>0.0], insym=[:i=>0.0], g=1:1,
                                    psym=[:S_b=>sb], name=nm)
        # a line whose end base copies from `dir` (:src / :dst); `ge=true` pins the graph
        # placement (needed only when the network is built from a vector of named vertices).
        mkline(dir; sb_default=nothing, ge=false) = begin
            gekw = ge ? (; src=:hub1, dst=:hub2) : (;)
            em = EdgeModel(; g=Directed(linef), outsym=[:P=>0.0],
                           insym=(src=[:us=>0.0], dst=[:ud=>0.0]),
                           psym=[:Y=>0.5, :S_b], name=:line, gekw...)
            set_default_from!(em, :S_b, (dir, :S_b))
            isnothing(sb_default) || set_default!(em, :S_b, sb_default)
            em
        end
        (; mkhub, mkline)
    finally
        NetworkDynamics.CHECK_COMPONENT[] = old
    end
end
M = _models()

@testset "default_from metadata survives Network construction" begin
    em = M.mkline(:src)
    @test has_default_from(em, :S_b)
    @test get_default_from(em, :S_b) == (:src, :S_b)
    nw = Network(path_graph(2), [M.mkhub(100.0, :hub1), M.mkhub(200.0, :hub2)], em)
    @test get_default_from(nw[EIndex(1)], :S_b) == (:src, :S_b)
end

@testset "src/dst copy, reinit tracking, and weak yield" begin
    # :src copies from vertex 1, :dst from vertex 2
    for (dir, exp) in ((:src, 100.0), (:dst, 200.0))
        nw = Network(path_graph(2), [M.mkhub(100.0, :hub1), M.mkhub(200.0, :hub2)], M.mkline(dir))
        s = initialize_componentwise(nw)
        @test s[EIndex(1, :S_b)] == exp
    end

    # reinit tracking: change the src base, the copied value must follow (weak output was not
    # frozen into a persistent default by the first init)
    nw = Network(path_graph(2), [M.mkhub(100.0, :hub1), M.mkhub(200.0, :hub2)], M.mkline(:src))
    @test initialize_componentwise(nw)[EIndex(1, :S_b)] == 100.0
    set_default!(nw[VIndex(1)], :S_b, 250.0)
    @test initialize_componentwise(nw)[EIndex(1, :S_b)] == 250.0

    # weak yield: a user-set default on the target wins on *every* init
    nw2 = Network(path_graph(2), [M.mkhub(100.0, :hub1), M.mkhub(200.0, :hub2)],
                  M.mkline(:src; sb_default=42.0))
    @test initialize_componentwise(nw2)[EIndex(1, :S_b)] == 42.0
    set_default!(nw2[VIndex(1)], :S_b, 250.0)
    @test initialize_componentwise(nw2)[EIndex(1, :S_b)] == 42.0
end

@testset "composition-layer override wins over source default" begin
    nw = Network(path_graph(2), [M.mkhub(100.0, :hub1), M.mkhub(200.0, :hub2)], M.mkline(:src))
    s = initialize_componentwise(nw; default_overrides=Dict(VIndex(1, :S_b) => 999.0))
    @test s[EIndex(1, :S_b)] == 999.0
end

@testset "valueless source is skipped (no formula built)" begin
    # source param carries no default and no override → nothing to copy, so no weak formula is
    # injected for the target (asserted at the resolution level; the target keeps its own default).
    old = NetworkDynamics.CHECK_COMPONENT[]; NetworkDynamics.CHECK_COMPONENT[] = false
    try
        srchub = VertexModel(f=hubf, sym=[:u=>0.0], insym=[:i=>0.0], g=1:1, psym=[:S_b], name=:hub1)
        nw = Network(path_graph(2), [srchub, M.mkhub(200.0, :hub2)], M.mkline(:src; sb_default=7.0))
        merged = resolve_default_from(nw, nothing, nothing)
        @test isempty(merged) || !haskey(merged, EIndex(1))
    finally
        NetworkDynamics.CHECK_COMPONENT[] = old
    end
end

@testset "full per-unit flow (Lib.pu_injector_network): hub base → injector + line, reinit" begin
    # Realistic hub/injector/line slice from the ComponentLibrary. Two hubs carry distinct system
    # bases (100 / 200); the base must flow into each injector (`default_from (:hub,…)`), each
    # injector's rating `Sn` (`initf_weak`, weak), and the line (`default_from (:src,…)`). The
    # hubs and the line additionally carry `bound_to` aliases (eliminated at compile) — so this
    # one network exercises all three per-unit features together.
    nw, s0 = Lib.pu_injector_network()
    ov = interface_values(s0)   # base is observable-only, so s0's interface == the powerflow's
    report(s) = (inj1_Sb=s[VIndex(:inj1, :S_b)], inj1_Sn=s[VIndex(:inj1, :Sn)],
                 inj2_Sb=s[VIndex(:inj2, :S_b)], inj2_Sn=s[VIndex(:inj2, :Sn)],
                 line_Sb=s[EIndex(:line, :S_b)])

    @test report(s0) == (inj1_Sb=100.0, inj1_Sn=100.0, inj2_Sb=200.0, inj2_Sn=200.0, line_Sb=100.0)

    # change hub1's base + reinit → inj1 (via :hub), its rating Sn (via initf_weak) and the line
    # (via :src) all follow; hub2's dependents (inj2) are untouched
    set_default!(nw[VIndex(1)], :busbar₊S_b, 350.0)
    s1 = initialize_componentwise!(nw; verbose=false, subverbose=false, default_overrides=ov)
    @test report(s1) == (inj1_Sb=350.0, inj1_Sn=350.0, inj2_Sb=200.0, inj2_Sn=200.0, line_Sb=350.0)

    # pin inj2's rating explicitly → the weak initf_weak yields, Sn stays put while S_b still tracks
    set_default!(nw[VIndex(4)], :Sn, 555.0)
    s2 = initialize_componentwise!(nw; verbose=false, subverbose=false, default_overrides=ov)
    @test report(s2) == (inj1_Sb=350.0, inj1_Sn=350.0, inj2_Sb=200.0, inj2_Sn=555.0, line_Sb=350.0)
end

@testset "NWState(nw) re-application drops weak formulas (yields to a user default)" begin
    # Regression guard for `_get_appropriate_dict` (`src/symbolicindexing.jl`): building an
    # `NWState(nw)` re-applies a component's init formulas, and must run `drop_weak_formulas`
    # exactly like `initialize_component` — otherwise a weak formula fires unconditionally and
    # clobbers a user default it is supposed to defer to. inj1 carries a weak `initf_weak` on `Sn`
    # (rating weakly follows the base `S_b`); base is 100. (Same shape a future `default_from`
    # applied at `NWState` construction would exercise — reuse this testset for it.)
    nw, _ = Lib.pu_injector_network()
    inj1 = nw[VIndex(3)]
    @test any(f.weak for f in NetworkDynamics.get_initformulas(inj1))   # the weak formula is present

    # no user default → the weak formula fires: Sn tracks the base
    @test NWState(nw)[VIndex(3, :Sn)] == NWState(nw)[VIndex(3, :S_b)] == 100.0

    # a user default on the weak target → the weak formula yields, the default survives
    set_default!(inj1, :Sn, 777.0)
    @test NWState(nw)[VIndex(3, :Sn)] == 777.0
end

@testset "NWState(nw) materializes default_from copies" begin
    # `default_from` lowers to a weak formula that is resolved *network-wide* (`resolve_default_from`);
    # the NWState constructor must run that pre-pass too, so the copied neighbor value shows up on
    # reconstruction — not only after a full init. The line's `S_b` copies from its src hub (100).
    nw = Network(path_graph(2), [M.mkhub(100.0, :hub1), M.mkhub(200.0, :hub2)], M.mkline(:src))
    @test NWState(nw)[EIndex(1, :S_b)] == 100.0

    # a user default on the target → the weak default_from yields, the user value survives
    nw2 = Network(path_graph(2), [M.mkhub(100.0, :hub1), M.mkhub(200.0, :hub2)],
                  M.mkline(:src; sb_default=42.0))
    @test NWState(nw2)[EIndex(1, :S_b)] == 42.0
end

@testset "structural / resolution errors" begin
    hubs() = [M.mkhub(100.0, :hub1), M.mkhub(200.0, :hub2)]

    # unresolvable source symbol → error naming available params
    let em = M.mkline(:src)
        set_default_from!(em, :S_b, (:src, :NOPE))
        nw = Network(path_graph(2), hubs(), em)
        @test_throws ArgumentError initialize_componentwise(nw)
    end

    # :src / :dst on a *vertex* parameter → error
    let
        old = NetworkDynamics.CHECK_COMPONENT[]; NetworkDynamics.CHECK_COMPONENT[] = false
        try
            bad = M.mkhub(100.0, :hub1); set_default_from!(bad, :S_b, (:src, :S_b))
            nw = Network(path_graph(2), [bad, M.mkhub(200.0, :hub2)], M.mkline(:dst))
            @test_throws ArgumentError initialize_componentwise(nw)
        finally
            NetworkDynamics.CHECK_COMPONENT[] = old
        end
    end

    # :hub on an *edge* parameter → error
    let em = M.mkline(:src)
        set_default_from!(em, :S_b, (:hub, :S_b))
        nw = Network(path_graph(2), hubs(), em)
        @test_throws ArgumentError initialize_componentwise(nw)
    end

    # :hub on a non-injector vertex → error
    let
        old = NetworkDynamics.CHECK_COMPONENT[]; NetworkDynamics.CHECK_COMPONENT[] = false
        try
            bad = M.mkhub(100.0, :hub1); set_default_from!(bad, :S_b, (:hub, :S_b))
            nw = Network(path_graph(2), [bad, M.mkhub(200.0, :hub2)], M.mkline(:dst))
            @test_throws ArgumentError initialize_componentwise(nw)
        finally
            NetworkDynamics.CHECK_COMPONENT[] = old
        end
    end

    # default_from on a non-parameter (a state) → error
    let
        old = NetworkDynamics.CHECK_COMPONENT[]; NetworkDynamics.CHECK_COMPONENT[] = false
        try
            bad = M.mkhub(100.0, :hub1); set_default_from!(bad, :u, (:hub, :S_b))
            nw = Network(path_graph(2), [bad, M.mkhub(200.0, :hub2)], M.mkline(:dst))
            @test_throws ArgumentError initialize_componentwise(nw)
        finally
            NetworkDynamics.CHECK_COMPONENT[] = old
        end
    end
end

@testset "component show surfaces default_from" begin
    str = sprint(show, MIME"text/plain"(), M.mkline(:src))
    @test occursin("default_from", str)
    @test occursin(":S_b ← src S_b", str)
end

# ---------------------------------------------------------------------------------------------
# MTK spelling: `@parameters S_b [default_from = (:src, :busbar₊S_b)]` must bridge into ND
# metadata, and compose with a `bound_to` internal quantity (eliminated within the edge).
# ---------------------------------------------------------------------------------------------
busbar_sys(; name) = (@parameters S_b=100.0; @variables ub(t) [guess=1.0]; System([D(ub) ~ S_b - ub], t; name))
function pu_bus(; name)
    @named busbar = busbar_sys()
    @variables u(t)=1.0 i(t) [guess=0.0] o(t) [guess=1.0]
    @parameters P=1.0
    @parameters S_b [bound_to = :busbar₊S_b]
    System([D(u) ~ -u + i + P/S_b, o ~ u], t; name, systems=[busbar])
end
function pu_line(; name)
    @variables src_u(t) [guess=1.0] dst_u(t) [guess=1.0] flow(t) [guess=0.0]
    @parameters Y=0.5
    @parameters S_b [default_from = (:src, :busbar₊S_b)]
    @parameters Y_pu [bound_to = :S_b]
    System([flow ~ Y * Y_pu * (src_u - dst_u) / S_b], t; name)
end

@testset "MTK spelling bridges to ND metadata + composes with bound_to" begin
    @named line = pu_line()
    em = EdgeModel(line, [:src_u], [:dst_u], AntiSymmetric([:flow]))
    @test :S_b ∈ psym(em)                       # default_from target stays a live parameter
    @test :Y_pu ∉ psym(em)                      # bound_to internal is eliminated …
    @test :Y_pu ∈ obssym(em)                    # … and reappears as an observable
    @test get_default_from(em, :S_b) == (:src, :busbar₊S_b)

    # resolution: the MTK-declared default_from resolves against the src bus's *namespaced*
    # busbar base and tracks it (asserted at the pre-pass level — the surrounding component
    # physics are underdetermined in isolation, and end-to-end init is covered by the plain
    # models above).
    busv() = VertexModel(pu_bus(name=:bus), [:i], [:o])
    nw = Network(path_graph(2), [busv(), busv()], em)
    bake(nw) = let m = resolve_default_from(nw, nothing, nothing)
        f = only(m[EIndex(1)])
        out = SymbolicView(zeros(1), f.outsym); f(out, SymbolicView(Float64[], Symbol[]))
        out[:S_b]
    end
    @test bake(nw) == 100.0
    set_default!(nw[VIndex(1)], :busbar₊S_b, 275.0)
    @test bake(nw) == 275.0
end
