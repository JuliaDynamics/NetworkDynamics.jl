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

function ddevertex!(dv,v,h_v,e_s,e_d,h_s,h_d,p,t)
    dv .= 0
    for h in h_s
        dv .-= h
    end
    for h in h_d
        dv .+= h
    end
end

odeedge! = (dl,l,v_s,v_d,p,t) -> dl .= 1000*(v_s - v_d - l)
staticedge! = (l,v_s,v_d,p,t) -> l .= v_s - v_d
ddeedge! = (de,e,h_e,v_s,v_d,h_s,h_d,p,t) -> de .= v_s - v_d

odevertex = ODEVertex(vertex!,2, [2 1; -1 1], nothing)
odeedge = ODEEdge(odeedge! ,2, [1 -1; 2 0.3], nothing)
staticedge = StaticEdge(staticedge!, 2)
ddevertex = DDEVertex(ddevertex!, 2, [0.01,0.02],[0.02,0.03])
ddeedge = DDEEdge(ddeedge!, 2)


vertexes = [ddevertex for v in vertices(g)]
edgices = [ddeedge for e in edges(g)]


test = network_dynamics(vertexes,edgices,g)

test_prob = DDEProblem(test,rand(70),(0.,2.))

test(rand(70),rand(70),nothing,0.)

test_sol = solve(test_prob)
