using NetworkDynamics
using LightGraphs
using LinearAlgebra
using SparseArrays
using DifferentialEquations
using Plots

struct SolAnalytic
    L
end
function (sa::SolAnalytic)(x0, p, t)
    exp(- t * sa.L) * x0
end

g = barabasi_albert(10,5)

L = laplacian_matrix(g)

sol_analytic = SolAnalytic(Symmetric(Array(L)))

diff_network_L = ODEFunction((dx, x, p, t) -> dx .= - L * x, analytic=sol_analytic)

x0 = rand(10)

prob_L = ODEProblem(diff_network_L,x0,(0.,5.))
sol_L = solve(prob_L)

println(sol_L.errors)

plot(sol_L)

sol_ana = zeros(length(sol_L.t), length(x0))
for i in 1:length(sol_L.t)
    sol_ana[i, :] .= sol_analytic(x0, nothing, sol_L.t[i])
end
plot(sol_L.t, sol_ana)

#Now for NetworkDynamics

diffusion_edge! = (e,v_s,v_d,p,t) -> @. e = v_s - v_d

function diffusion_vertex!(dv, v, e_s, e_d, p, t)
    dv .= edge_sum(e_s, e_d) # Oriented sum of the incoming and outgoing edges
end

odevertex = ODEVertex(f! = diffusion_vertex!,dim = 1)
staticedge = StaticEdge(f! = diffusion_edge!, dim = 1)
odeedge = ODEEdge(staticedge) # We promote the static edge to an ODEEdge artifically

vertex_list = [odevertex for v in vertices(g)]
edge_list = [staticedge for e in edges(g)]

diff_network_st = network_dynamics(vertex_list,edge_list,g)

x0 = rand(10)

prob_st = ODEProblem(diff_network_st,x0,(0.,5.))

sol_st = solve(prob_st)

plot(sol_st)
