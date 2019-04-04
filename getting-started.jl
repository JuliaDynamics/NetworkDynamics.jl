include("src/NetworkDynamics.jl")
using .NetworkDynamics
using LightGraphs
using LinearAlgebra
using DifferentialEquations

g = barabasi_albert(10,5)

#= vertex! is basically dv = sum(e_d) - sum(e_s), so basically simple diffusion with the addition
of staticedge! and odeedge! below. =#
function vertex!(dv, v, e_s, e_d, p, t)
    # Note that e_ss and e_ds might be empty, the code needs to be able to deal
    # with this situation.
    dv .= 0
    for e in e_s
        dv .-= e
    end
    for e in e_d
        dv .+= e
    end
    nothing
end

#Pls ignore everything that has a dde as a prefix for now :).
function ddevertex!(dv,v,h_v,e_s,e_d,h_s,h_d,p,t)
    dv .= -v
    for h in e_s
        dv .-= h
    end
    for h in e_d
        dv .+= h
    end
end

odeedge! = (dl,l,v_s,v_d,p,t) -> dl .= 1000*(v_s - v_d - l)
staticedge! = (l,v_s,v_d,p,t) -> l .= v_s - v_d
ddeedge! = (de,e,h_e,v_s,v_d,h_s,h_d,p,t) -> de .= v_s - v_d

# We construct the Vertices and Edges with dimension 2.

odevertex = ODEVertex(vertex!,2)
odeedge = ODEEdge(odeedge! ,2)
staticedge = StaticEdge(staticedge!, 2)
ddevertex = DDEVertex(ddevertex!, w, 0, 0)
ddeedge = DDEEdge(ddeedge!, 1)


vertexes = [odevertex for v in vertices(g)]
edgices = [odeedge for e in edges(g)]


test = network_dynamics(vertexes,edgices,g)

x0 = rand(70)
h0(p,t; idxs = 1:70) = x0

test_prob = ODEProblem(test,x0,(0.,20.))

test(x0,x0,nothing,0.)

test_sol = solve(test_prob)

using Plots

plot(test_sol, legend = false, vars = 1:20)
