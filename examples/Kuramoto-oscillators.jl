using Pkg
Pkg.activate(@__DIR__)
using Revise

using NetworkDynamics
using LightGraphs
using LinearAlgebra
using OrdinaryDiffEq

g = barabasi_albert(50,5)

#= vertex! is basically dv = sum(e_d) - sum(e_s), so basically simple diffusion with the addition
of staticedge! and odeedge! below. =#

@inline function kuramoto_edge!(e,v_s,v_d,p,t)
    e[1] = sin(v_s[1] - v_d[1]) * (1. + t * p.T_inv)
    # e[2] = sin(v_d[1] - v_s[1])
    nothing
end

struct kuramoto_parameters
    ω
    T_inv
end

@inline function kuramoto_vertex!(dv, v, e_s, e_d, p, t)
    # Note that e_s and e_d might be empty, the code needs to be able to deal
    # with this situation.
    dv .= p.ω
    for e in e_s
        dv .-= e[1]
    end
    for e in e_d
        dv .+= e[1] # .-= e[2]
    end
    nothing
end

odevertex = ODEVertex(f! = kuramoto_vertex!, dim = 1)
staticedge = StaticEdge(f! = kuramoto_edge!, dim = 1)

parameters = [kuramoto_parameters(10. * randn(), 0.) for v in vertices(g)]
append!(parameters, [kuramoto_parameters(0., 1. /30.) for v in edges(g)])

kuramoto_network! = network_dynamics(odevertex,staticedge,g)

x0 = randn(nv(g))
dx = similar(x0)

kuramoto_network!(dx, x0, parameters, 0.)

prob = ODEProblem(kuramoto_network!, x0, (0.,150.), parameters)

sol = solve(prob, Rodas4()) # ForwardDiff error
sol2 = solve(prob, Tsit5())

using Plots

plot(sol, tspan=(147.,150.))
plot(sol.t[:],mod2pi.(sol'[:,:]),legend=false)
plot(sol2)

using Profile
using ProfileView

@profile sol = solve(prob, Tsit5())

ProfileView.view()
