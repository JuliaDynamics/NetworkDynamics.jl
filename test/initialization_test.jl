using NetworkDynamics, Graphs
using SteadyStateDiffEq, OrdinaryDiffEq


@testset "test find_fixpoint" begin

    function swing_equation!(dv, v, esum, (M, P, D), t)
        dv[1] = v[2]
        dv[2] = 1/M *(P - D * v[2] + esum[1])
        nothing
    end
    swing_vertex = ODEVertex(swing_equation!; sym=[:θ, :ω], psym=[:M=>1, :P, :D=>0.1])

    function powerflow!(e, v_s, v_d, (K,), t)
        e[1] = K * sin(v_s[1] - v_d[1])
    end
    powerflow_edge = StaticEdge(powerflow!; dim=1, psym=[:K=>6], coupling=AntiSymmetric())
    g = watts_strogatz(4, 2, 0.0)

    nd = Network(g, swing_vertex, powerflow_edge)
    p = NWParameter(nd)
    p.v[1:4, :P] = [1, -1, 1, -1]

    s0 = find_fixpoint(nd, p)
    @test pflat(s0) == pflat(p)

    s1 = find_fixpoint(nd, p; alg=DynamicSS(Rodas5P()))
    @test diff(s0.v[1:nv(g),:θ]) ≈ diff(s1.v[1:nv(g),:θ])
end
