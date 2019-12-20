# The most elaborate example is Kuramoto-inertial.jl check that out for an
# overiew of the package.

using Pkg
Pkg.activate(@__DIR__)
using Revise

using NetworkDynamics
using LightGraphs
using LinearAlgebra
using OrdinaryDiffEq

g = barabasi_albert(10,5)

#= vertex! is basically dv = sum(e_d) - sum(e_s), so basically simple diffusion with the addition
of staticedge! and odeedge! below. =#

function vertex!(dv, v, e_s, e_d, p, t)
    # Note that e_s and e_d might be empty, the code needs to be able to deal
    # with this situation.
    dv .= 0.
    for e in e_s
        dv .-= e
    end
    for e in e_d
        dv .+= e
    end
    nothing
end

odeedge! = (de,e,v_s,v_d,p,t) -> de .= 1000. * (v_s - v_d - e)
staticedge! = (e,v_s,v_d,p,t) -> e .= v_s - v_d

# We construct the Vertices and Edges with dimension 2. This means we will have
# two parallel diffusions on the network.

odevertex = ODEVertex(f! = vertex!,dim = 2, sym = [:v,:w])
odeedge = ODEEdge(f! = odeedge!, dim = 2, sym = [:v,:w])

staticedge = StaticEdge(f! = staticedge!, dim = 2, sym = [:v,:w])

test = network_dynamics(odevertex, staticedge, g)

x0 = rand(20)

test_prob = ODEProblem(test,x0,(0.,5.))

test_sol = solve(test_prob, Tsit5())

using Plots

plot(test_sol, vars = syms_containing(test, "v"))
plot(test_sol, vars = syms_containing(test, "w"))
