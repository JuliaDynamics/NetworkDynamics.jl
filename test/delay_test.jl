using Test
using NetworkDynamics
using LightGraphs
using OrdinaryDiffEq
using DelayDiffEq

N = 10
g = watts_strogatz(N,2,0.)

@inline function diffusionedge!(e, v_s, v_d, p, t)
    e[1] = v_s[1] - v_d[1]
    nothing
end

@inline function delayedge!(e, v_s, v_d, h_v_s, h_v_d, p, t)
    e[1] = h_v_s[1] - v_d[1]
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
nd_delay_edge = StaticDelayEdge(f! = delayedge!, dim = 1)


nd! = network_dynamics(nd_diffusion_vertex, nd_diffusion_edge, g)

x0 = [ones(N÷2); -ones(N÷2)]
tspan = (0.,200.)
prob = ODEProblem(nd!, x0, tspan, nothing)
sol = solve(prob, Tsit5(), abstol=1e-6)


@test isapprox(sol[end], zeros(N), atol=1e-6)



dnd! = network_dynamics(nd_diffusion_vertex, nd_delay_edge, g)

x0 = [ones(N÷2); -ones(N÷2)]
tspan = (0.,200.)
# constant initial conditions on interval
h!(out, p, t) = out .= x0
# (vertexp, edgep, τ)
p = (nothing, nothing, 1.)

prob = DDEProblem(dnd!, x0, h!, tspan, p)#, constant_lags=[p[end]])
sol = solve(prob, MethodOfSteps(Tsit5()), abstol=1e-6)

isapprox(sol[end], zeros(N), atol=1e-6)

dx0 = similar(x0)

nd!.f

@code_warntype dnd!.f(dx0, x0,nothing,0.)

@allocated nd!(dx0, x0,nothing,0.)
@allocated dnd!(dx0, x0, p, 0.)
