import NetworkDynamics
using LightGraphs
using LinearAlgebra
using SparseArrays
using DifferentialEquations
using Plots
using BenchmarkTools

const ND = NetworkDynamics

g = barabasi_albert(10^4,5)
L = laplacian_matrix(g)

function diffusion_dyn(dx, x, p, t)
    dx .= - L * x
    nothing
end

diff_network_L = ODEFunction(diffusion_dyn)

x0 = rand(nv(g))

# prob_L = ODEProblem(diff_network_L,x0,(0.,5.))
# sol_L = solve(prob_L)

#Now for NetworkDynamics

@inline function diffusion_edge!(e,v_s,v_d,p,t)
    e .= v_s .- v_d
    nothing
end

@inline function diffusion_vertex!(dv, v, e_s, e_d, p, t)
    dv .= 0.
    ND.edge_sum!(dv, e_s, e_d) # Oriented sum of the incoming and outgoing edges
    nothing
end

odevertex = ND.ODEVertex(f! = diffusion_vertex!,dim = 1)
staticedge = ND.StaticEdge(f! = diffusion_edge!, dim = 1)

vertex_list = [odevertex for v in vertices(g)]
edge_list = [staticedge for e in edges(g)]

diff_network_st = ND.network_dynamics(vertex_list,edge_list,g)

x0 = rand(nv(g))
dx_L = similar(x0)
dx_st = similar(x0)

@code_warntype diff_network_st(dx_st, x0, nothing, 0.)

diff_network_st(dx_st, x0, nothing, 0.)
diff_network_L(dx_L, x0, nothing, 0.)

println("Benachmarking")

@btime diff_network_st(dx_st, x0, nothing, 0.)
@btime diff_network_L(dx_L, x0, nothing, 0.)
# @btime diffusion_dyn(dx_L, x0, nothing, 0.)
