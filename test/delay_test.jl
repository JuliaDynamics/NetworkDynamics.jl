using Test
using NetworkDynamics
using LightGraphs
using OrdinaryDiffEq
using DelayDiffEq
#using BenchmarkTools

N = 10
g = complete_graph(N) # more edges than vertices!

@inline function diffusionedge!(e, v_s, v_d, p, t)
    e[1] = v_s[1] - v_d[1]
    nothing
end



@inline function diffusionvertex!(dv, v, edges, p, t)
    dv[1] = 0.
    for e in edges
        dv[1] += e[1]
    end
    nothing
end

nd_diffusion_vertex = ODEVertex(f! = diffusionvertex!, dim = 1)
nd_diffusion_edge = StaticEdge(f! = diffusionedge!, dim = 1, coupling = :undirected)
nd! = network_dynamics(nd_diffusion_vertex, nd_diffusion_edge, g)


@testset "ODE diffusion" begin
    nd! = network_dynamics(nd_diffusion_vertex, nd_diffusion_edge, g)

    x0 = [ones(N÷2); -ones(N÷2)]
    tspan = (0.,100.)

    dx0 = similar(x0)
    #@btime $nd!.f($dx0,$x0,nothing,0.)
    prob = ODEProblem(nd!, x0, tspan, nothing)
    sol = solve(prob, Tsit5())

    #@btime solve($prob, Tsit5(), abstol=1e-6)
    @test isapprox(sol[end], zeros(N), atol=1e-6)
end

# @testset "Homogeneous DDE diffusion old" begin
#     @inline function delayedge!(e, v_s, v_d, h_v_s, h_v_d, p, t)
#         e[1] = h_v_s[1] - v_d[1]
#         nothing
#     end
#     nd_delay_edge = StaticDelayEdge(f! = delayedge!, dim = 1)
#
#     dnd! = network_dynamics(nd_diffusion_vertex, nd_delay_edge, g)
#
#     x0 = [ones(N÷2); -ones(N÷2)]
#     tspan = (0.,100.)

#     # constant initial conditions on interval
#     h!(out, p, t) = out .= x0
#     # (vertexp, edgep, τ)
#     p = (nothing, nothing, 1.)
#
#     dx0 = similar(x0)
#     # @btime $dnd!.f($dx0,$x0,$h!,$p,0.)
#
#     prob = DDEProblem(dnd!, x0, h!, tspan, p)#, constant_lags=[p[end]])
#     sol = solve(prob, MethodOfSteps(Tsit5()), abstol=1e-6)
#     # @btime solve($prob, MethodOfSteps(Tsit5()), abstol=1e-6)
#     @test isapprox(sol[end], zeros(N), atol=1e-4)
# end



@testset "Homogeneous DDE diffusion" begin
    @inline function delayedge!(e, v_s, v_d, h, v_s_idx, v_d_idx, p, t)
        hist1 = h(p, t - 1., idxs=v_s_idx[1])
        e[1] = hist1 - v_d[1]
        nothing
    end

    nd_delay_edge = StaticDelayEdge(f! = delayedge!, dim = 1, coupling = :undirected)

    dnd! = network_dynamics(nd_diffusion_vertex, nd_delay_edge, g)
    dnd! = ODEFunction(dnd!.f)
    x0 = [ones(N÷2); -ones(N÷2)]
    tspan = (0.,100.)

    # constant initial conditions on interval
    h(p, t; idxs=nothing) = x0[idxs]
    # (vertexp, edgep, τ)

    dx0 = similar(x0)
    #@btime $dnd!.f($dx0,$x0,$h,nothing,0.)
    prob = DDEProblem(dnd!, x0, h, tspan, constant_lags=[1.])
    sol = solve(prob, MethodOfSteps(Tsit5()))
    #@btime solve($prob, MethodOfSteps(Tsit5()))

    @test isapprox(sol[end], zeros(N), atol=1e-6)
end

@testset "Heterogeneous DDE diffusion" begin
    @inline function delayedge!(e, v_s, v_d, h, v_s_idx, v_d_idx, p, t)
        hist1 = h(p, t - p, idxs=v_s_idx[1])
        e[1] = hist1 - v_d[1]
        nothing
    end

    nd_delay_edge = StaticDelayEdge(f! = delayedge!, dim = 1, coupling = :undirected)

    dnd! = network_dynamics(nd_diffusion_vertex, nd_delay_edge, g)

    x0 = [ones(N÷2); -ones(N÷2)]
    tspan = (0.,100.)

    # constant initial conditions on interval
    h(p, t; idxs=nothing) = 0.
    # (vertexp, edgep)
    p = (nothing, collect(1:ne(g)) ./ ne(g))

    dx0 = similar(x0)
    #@btime $dnd!.f($dx0,$x0,$h,1.,0.)

    prob = DDEProblem(dnd!, x0, h, tspan, p, constant_lags = p[2])
    sol = solve(prob, MethodOfSteps(Tsit5()))

    # for too many lags, declaring lags is very slow
    #prob = DDEProblem(dnd!, x0, h, tspan, p)
    #@btime solve($prob, MethodOfSteps(Tsit5()),save_everystep=false)

    # With inhomogenities we lose a lot of precision
    # increase reltol?
    @test isapprox(sol[end], zeros(N), atol=1e-2)
end
