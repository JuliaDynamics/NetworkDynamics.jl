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
    e[1] = sin(v_s[1] - v_d[1]) * (1. + t * p)
    # e[2] = sin(v_d[1] - v_s[1])
    nothing
end

@inline function kuramoto_vertex!(dv, v, e_s, e_d, p, t)
    # Note that e_s and e_d might be empty, the code needs to be able to deal
    # with this situation.
    dv .= p
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

v_pars = [1. * randn() for v in vertices(g)]
e_pars = [1. /3. for e in edges(g)]

parameters = (v_pars, e_pars)

kuramoto_network! = network_dynamics(odevertex,staticedge,g)

x0 = randn(nv(g))
dx = similar(x0)

kuramoto_network!(dx, x0, parameters, 0.)

prob = ODEProblem(kuramoto_network!, x0, (0.,15.), parameters)

sol = solve(prob, Rodas4P()) # ForwardDiff error
sol2 = solve(prob, Tsit5())
sol3 = solve(prob, TRBDF2())

using Plots

plot(sol)
plot(sol2)
plot(sol3)

using Profile
using ProfileView

@profile sol = solve(prob, Tsit5())

ProfileView.view()
