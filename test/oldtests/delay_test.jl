using Test
using NetworkDynamics
using Graphs
using OrdinaryDiffEq
using DelayDiffEq

N = 10
g = complete_graph(N) # more edges than vertices!

@inline function diffusionedge!(e, v_s, v_d, p, t)
    e[1] = v_s[1] - v_d[1]
    nothing
end

@inline function diffusionvertex!(dv, v, esum, p, t)
    dv[1] = esum[1]
    nothing
end

nd_diffusion_vertex = ODEVertex(; f=diffusionvertex!, dim=1, pdim=0)
nd_diffusion_edge = StaticEdge(; f=diffusionedge!, dim=1, pdim=0, coupling=AntiSymmetric())
nd! = Network(g, nd_diffusion_vertex, nd_diffusion_edge)


@testset "ODE diffusion" begin
    nd! = Network(g, nd_diffusion_vertex, nd_diffusion_edge)

    x0 = [ones(N ÷ 2); -ones(N ÷ 2)]
    tspan = (0.0, 100.0)

    dx0 = similar(x0)
    prob = ODEProblem(nd!, x0, tspan, nothing)
    sol = solve(prob, Tsit5())

    #@btime solve($prob, Tsit5(), abstol=1e-6)
    @test isapprox(sol[end], zeros(N); atol=1e-6)
end

using NetworkDynamics.KernelAbstractions

@kernel function ekernel!(@Const(batch), du, @Const(u), @Const(srcrange), @Const(dstrange), @Const(p), @Const(t))
    I = @index(Global)
    @inline apply_edge_unbuffered!(type, batch, I, du, u, srcrange, dstrange, p, t)
end
@kernel function ekernel!(@Const(batch), (du, @Const(u), @Const(srcrange), @Const(dstrange), @Const(p)), @Const(t))
    I = @index(Global)
    @inline apply_edge_unbuffered!(type, batch, I, du, u, srcrange, dstrange, p, t)
end

NetworkDynamics.executionstyle(nd!)



@testset "Homogeneous DDE diffusion" begin
    @inline function delayedge!(e, v_s, v_d, h_v_s, h_v_d, p, t)
        hist1 = h_v_s(t - 1.0; idxs=1)
        e[1] = hist1 - v_d[1]
        nothing
    end

    nd_delay_edge = StaticDelayEdge(; f=delayedge!, dim=1, coupling=:undirected)

    dnd! = network_dynamics(nd_diffusion_vertex, nd_delay_edge, g)
    x0 = [ones(N ÷ 2); -ones(N ÷ 2)]
    tspan = (0.0, 100.0)

    # constant initial conditions on interval
    h(p, t; idxs=nothing) = x0[idxs]

    dx0 = similar(x0)
    #@btime $dnd!.f($dx0,$x0,$h,nothing,0.)
    prob = DDEProblem(dnd!, x0, h, tspan; constant_lags=[1.0])
    sol = solve(prob, MethodOfSteps(Tsit5()))
    #@btime solve($prob, MethodOfSteps(Tsit5()))

    @test isapprox(sol[end], zeros(N); atol=1e-5)
end


@testset "Vertex DDE diffusion" begin

    @inline function delayvertex!(dv, v, edges, h_v, p, t)
        dv[1] = h_v(t - 1.0; idxs=1)
        for e in edges
            dv[1] += e[1]
        end
        nothing
    end
    nd_delay_vertex = DDEVertex(; f=delayvertex!, dim=1)

    dnd! = network_dynamics(nd_delay_vertex, nd_diffusion_edge, g)
    dnd! = ODEFunction(dnd!.f)
    x0 = [ones(N ÷ 2); -ones(N ÷ 2)]
    tspan = (0.0, 100.0)

    # constant initial conditions on interval
    h(p, t; idxs=nothing) = x0[idxs]
    dx0 = similar(x0)
    #@btime $dnd!.f($dx0,$x0,$h,nothing,0.)
    prob = DDEProblem(dnd!, x0, h, tspan; constant_lags=[1.0])
    sol = solve(prob, MethodOfSteps(Tsit5()))
    #@btime solve($prob, MethodOfSteps(Tsit5()))

    @test isapprox(sol[end], zeros(N); atol=1e-5)
end

@testset "Heterogeneous DDE diffusion" begin
    @inline function delayedge!(e, v_s, v_d, h_v_s, h_v_d, p, t)
        hist1 = h_v_s(t - p; idxs=1)
        e[1] = hist1 - v_d[1]
        nothing
    end

    nd_delay_edge = StaticDelayEdge(; f=delayedge!, dim=1, coupling=:undirected)

    dnd! = network_dynamics(nd_diffusion_vertex, nd_delay_edge, g)

    x0 = [ones(N ÷ 2); -ones(N ÷ 2)]
    tspan = (0.0, 100.0)

    # constant initial conditions on interval
    h(p, t; idxs=nothing) = 0.0
    # (vertexp, edgep)
    p = (nothing, collect(1:ne(g)) ./ ne(g))

    dx0 = similar(x0)
    #@btime $dnd!.f($dx0,$x0,$h,1.,0.)

    prob = DDEProblem(dnd!, x0, h, tspan, p; constant_lags=p[2])
    sol = solve(prob, MethodOfSteps(Tsit5()))

    # for too many lags, declaring lags is very slow
    #prob = DDEProblem(dnd!, x0, h, tspan, p)
    #@btime solve($prob, MethodOfSteps(Tsit5()),save_everystep=false)

    # Not sure if the failure to reach 0 is due to loss of precision or due to the problem formulation
    @test isapprox(sol[end], zeros(N); atol=1e-2)
end
