using Test
using NetworkDynamics
using LightGraphs
using OrdinaryDiffEq
using DelayDiffEq
using BenchmarkTools

N = 10
g = watts_strogatz(N,4,0.) # more edges than vertices!

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
nd_diffusion_edge = StaticEdge(f! = diffusionedge!, dim = 1, coupling = :antisymmetric)
nd! = network_dynamics(nd_diffusion_vertex, nd_diffusion_edge, g)


@testset "ODE diffusion" begin
    nd! = network_dynamics(nd_diffusion_vertex, nd_diffusion_edge, g)

    x0 = [ones(N÷2); -ones(N÷2)]
    tspan = (0.,5000.)
    dx0 = similar(x0)
    @btime $nd!.f($dx0,$x0,nothing,0.)
    prob = ODEProblem(nd!, x0, tspan, nothing)
    sol = solve(prob, Tsit5(), abstol=1e-6)

    @btime solve($prob, Tsit5(), abstol=1e-6)
    @test isapprox(sol[end], zeros(N), atol=1e-4)
end

@testset "Homogeneous DDE diffusion" begin
    @inline function delayedge!(e, v_s, v_d, h_v_s, h_v_d, p, t)
        e[1] = h_v_s[1] - v_d[1]
        nothing
    end
    nd_delay_edge = StaticDelayEdge(f! = delayedge!, dim = 1)

    dnd! = network_dynamics(nd_diffusion_vertex, nd_delay_edge, g)

    x0 = [ones(N÷2); -ones(N÷2)]
    tspan = (0.,5000.)
    # constant initial conditions on interval
    h!(out, p, t) = out .= x0
    # (vertexp, edgep, τ)
    p = (nothing, nothing, 1.)

    dx0 = similar(x0)
    @btime $dnd!.f($dx0,$x0,$h!,$p,0.)

    prob = DDEProblem(dnd!, x0, h!, tspan, p)#, constant_lags=[p[end]])
    sol = solve(prob, MethodOfSteps(Tsit5()), abstol=1e-6)
    @btime solve($prob, MethodOfSteps(Tsit5()), abstol=1e-6)
    @test isapprox(sol[end], zeros(N), atol=1e-4)
end



@testset "New homogeneous DDE diffusion" begin
    @inline function delayedge!(e, v_s, v_d, h, v_s_idx, v_d_idx, p, t)
        hist1 = h(nothing, t-1., idxs=v_s_idx[1])::Float64
        e[1] = hist1 - v_d[1]
        hist2 = h(nothing, t-1., idxs=v_d_idx[1])::Float64
        e[2] = hist2 - v_s[1]
        nothing
    end

    nd_delay_edge = StaticDelayEdge(f! = delayedge!, dim = 2, coupling = :fiducial)

    dnd! = network_dynamics(nd_diffusion_vertex, nd_delay_edge, g)

    x0 = [ones(N÷2); -ones(N÷2)]
    tspan = (0.,5000.)
    # constant initial conditions on interval
    h(p, t; idxs=nothing) = x0[idxs]
    # (vertexp, edgep, τ)
    p = (nothing, nothing, 1.)

    dx0 = similar(x0)
    @btime $dnd!.f($dx0,$x0,$h,$p,0.)

    prob = DDEProblem(dnd!, x0, h, tspan, p)#, constant_lags=[p[end]])
    sol = solve(prob, MethodOfSteps(Tsit5()), abstol=1e-6)
    @btime solve($prob, MethodOfSteps(Tsit5()), abstol=1e-6)

    @test isapprox(sol[end], zeros(N), atol=1e-3)
end



## Allocations

@inline function delayedge!(e, v_s, v_d, f::F, v_s_idx, v_d_idx, p, t) where F
    hist1 = f(nothing, t-1., idxs=v_s_idx[1])::Float64
    e[1] = hist1 - v_d[1]
    hist2 = f(nothing, t-1., idxs=v_d_idx[1])::Float64
    e[2] = hist2 - v_s[1]
    nothing
end

nd_delay_edge = StaticDelayEdge(f! = delayedge!, dim = 2, coupling = :fiducial)


begin
    # (vertexp, edgep, τ)
    p = (nothing, nothing, 1.)
    x0 = [ones(N÷2); -ones(N÷2)]
    const x1 = copy(x0)
    dx0 = similar(x0)
    const dx1 = similar(dx0)
end
function hh(p::Nothing, t::Float64; idxs::Int64=nothing)::Float64
    return x1[idxs]
end
@btime delayedge!(dx1,dx1,dx1,hh,1:1,1:1,p,0)


dnd! = network_dynamics(nd_diffusion_vertex, nd_delay_edge, g)

@btime $dnd!.f($dx1,$x1,$hh,$p,0.)
#@btime $nd!.f($dx0,$x0,$hh,$p,0.)

tspan = (0.,5000.)
prob = DDEProblem(dnd!, x1, hh, tspan, p)#, constant_lags=[p[end]])
sol = solve(prob, MethodOfSteps(Tsit5()), abstol=1e-6)
#@btime solve($prob, MethodOfSteps(Tsit5()), abstol=1e-6)

@test isapprox(sol[end], zeros(N), atol=1e-3)
