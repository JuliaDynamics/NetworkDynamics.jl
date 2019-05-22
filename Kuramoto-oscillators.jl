include("src/NetworkDynamics.jl")
using .NetworkDynamics
using LightGraphs
using LinearAlgebra
using DifferentialEquations

g = barabasi_albert(500,5)

#= vertex! is basically dv = sum(e_d) - sum(e_s), so basically simple diffusion with the addition
of staticedge! and odeedge! below. =#

function kuramoto_edge!(e,v_s,v_d,p,t)
    e[1] = sin(v_s[1] - v_d[1])
    # e[2] = sin(v_d[1] - v_s[1])
    nothing
end

struct kuramoto_parameters
    ω
end

function kuramoto_vertex!(dv, v, e_s, e_d, p, t)
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

vertexes = [odevertex for v in vertices(g)]
edgices = [staticedge for e in edges(g)]

parameters = [kuramoto_parameters(10. * randn()) for v in vertices(g)]
append!(parameters, [kuramoto_parameters(0.) for v in edges(g)])

kuramoto_network! = network_dynamics(vertexes,edgices,g)

using ForwardDiff
T_dual = ForwardDiff.Dual{nothing,Float64, ForwardDiff.pickchunksize(length(kuramoto_network!.f.e_int))}
kuramoto_network! = network_dynamics(vertexes,edgices,g; T=T_dual)

x0 = randn(nv(g))
dx = similar(x0)

kuramoto_network!(dx, x0, parameters, 0.)

prob = ODEProblem(kuramoto_network!, x0, (0.,5.), parameters)

sol = solve(prob, Rodas4())

using Profile
using ProfileView

@profile sol = solve(prob)

ProfileView.view()

using Plots

plot(sol, vars=1:10:500)
