using NetworkDynamics
using LightGraphs
using LinearAlgebra
using DifferentialEquations

g = barabasi_albert(10,5)

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

odeedge! = (dl,l,v_s,v_d,p,t) -> dl .= 1000*(v_s - v_d - l)
staticedge! = (l,v_s,v_d,p,t) -> l .= v_s - v_d

odevertex = ODEVertex(vertex!,2, [2 1; -1 1], nothing)
odeedge = ODEEdge(odeedge! ,2, [1 -1; 2 0.3], nothing)
staticedge = StaticEdge(staticedge!, 2)

vertexes = [odevertex for v in vertices(g)]
edgices = [staticedge for e in edges(g)]


test = network_dynamics(vertexes,edgices,g)

test_prob = ODEProblem(test,rand(20),(0.,2.))

test(rand(70),rand(70),nothing,0.)

test_sol = solve(test_prob, Rodas4(autodiff = false))

@assert 1 == 1

using NetworkDynamics
