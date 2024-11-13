using NetworkDynamics, Graphs
using SteadyStateDiffEq, OrdinaryDiffEqRosenbrock
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as Dt


@testset "test find_fixpoint" begin
    function swing_equation!(dv, v, esum, (M, P, D), t)
        dv[1] = v[2]
        dv[2] = 1/M *(P - D * v[2] + esum[1])
        nothing
    end
    swing_vertex = VertexFunction(f=swing_equation!, g=1, sym=[:θ, :ω], psym=[:M=>1, :P, :D=>0.1])

    function powerflow!(e, v_s, v_d, (K,), t)
        e[1] = K * sin(v_s[1] - v_d[1])
    end
    powerflow_edge = EdgeFunction(g=AntiSymmetric(powerflow!); outdim=1, psym=[:K=>6])
    g = watts_strogatz(4, 2, 0.0)

    nd = Network(g, swing_vertex, powerflow_edge)
    p = NWParameter(nd)
    p.v[1:4, :P] = [1, -1, 1, -1]

    s0 = find_fixpoint(nd, p)
    @test pflat(s0) == pflat(p)

    s1 = find_fixpoint(nd, p; alg=DynamicSS(Rodas5P()))
    @test isapprox(diff(s0.v[1:nv(g),:θ]), diff(s1.v[1:nv(g),:θ]); atol=1e-8)
end

@testset "test component initialization" begin
    @mtkmodel InitSwing begin
        @variables begin
            u_r(t)=1, [description="bus d-voltage", output=true]
            u_i(t)=0.1, [description="bus q-voltage", output=true]
            i_r(t)=1, [description="bus d-current (flowing into bus)", input=true]
            i_i(t)=0.1, [description="bus d-current (flowing into bus)", input=true]
            ω(t), [guess=0.0, description="Rotor frequency"]
            θ(t), [guess=0.0, bounds=[-π, π], description="Rotor angle"]
            Pel(t), [guess=1, description="Electrical Power injected into the grid"]
        end
        @parameters begin
            M=0.005, [description="Inertia"]
            D=0.1, [description="Damping"]
            V=sqrt(u_r^2 + u_i^2), [description="Voltage magnitude"]
            ω_ref=0, [description="Reference frequency"]
            Pm, [guess=0.1,description="Mechanical Power"]
        end
        @equations begin
            Dt(θ) ~ ω - ω_ref
            Dt(ω) ~ 1/M * (Pm - D*ω - Pel)
            Pel ~ u_r*i_r + u_i*i_i
            u_r ~ V*cos(θ)
            u_i ~ V*sin(θ)
        end
    end
    sys = InitSwing(name=:swing)
    vf = VertexFunction(sys, [:i_r, :i_i], [:u_r, :u_i])

    @test vf.symmetadata[:u_r][:default] == 1
    @test vf.symmetadata[:u_i][:default] == 0.1
    @test vf.symmetadata[:Pm][:guess] == 0.1
    @test vf.symmetadata[:θ][:bounds] == [-π, π]
    @test vf.symmetadata[:i_r][:default] == 1
    @test vf.symmetadata[:i_i][:default] == 0.1

    NetworkDynamics.initialize_component!(vf)
    @test NetworkDynamics.init_residual(vf) < 1e-8
    @test init_residual(vf) ≈ init_residual(vf; recalc=true)
end
