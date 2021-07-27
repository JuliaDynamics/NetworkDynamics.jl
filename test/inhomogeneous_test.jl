using Test
using NetworkDynamics
using LightGraphs
using LinearAlgebra
using SparseArrays
using OrdinaryDiffEq

# We will work on a random graph:
g = barabasi_albert(10,5)
L = laplacian_matrix(g)

N = nv(g)

printstyled("--- Network dynamics with a static vertex ---\n", bold=true, color=:white)


@inline function diffusion_edge!(e,v_s,v_d,p,t)
    e .= v_s .- v_d
    nothing
end
@inline function diffusion_edge2!(e,v_s,v_d,p,t)
    e .= v_s .- v_d
    nothing
end
@inline function diffusion_vertex!(dv, v, edges, p, t)
    dv .= 0.
    sum_coupling!(dv, edges) # sum incoming edges
    nothing
end

statvertex = StaticVertex(f! = (v, edges, p, t) -> v .= pi, dim = 1)

odevertex = ODEVertex(f! = diffusion_vertex!, dim = 1)
staticedge = StaticEdge(f! = diffusion_edge!, dim = 1, coupling = :antisymmetric)
staticedge2 = StaticEdge(f! = diffusion_edge2!, dim = 1, coupling = :antisymmetric)

vertex_list = [statvertex, odevertex]
append!(vertex_list, [odevertex for i in 1:N-2])
edge_list = StaticEdge[staticedge for e in edges(g)]
edge_list[2] = staticedge2

diff_network_st_ver = network_dynamics(vertex_list, edge_list, g)

x0 = rand(nv(g))
x0[1] = pi


prob_st_ver = ODEProblem(diff_network_st_ver, x0, (0.,500.))
@time sol_st_ver = solve(prob_st_ver, Rodas4())


function test(vertex_list, edge_list, g)
    nd = network_dynamics(vertex_list, edge_list, g)

    x0 = rand(nv(g))
    x0[1] = pi


    prob_st_ver = ODEProblem(nd, x0, (0.,500.))
    @time sol_st_ver = solve(prob_st_ver, Rodas4())
    @time sol_st_ver = solve(prob_st_ver, Rodas4(), saveat=[250,500])
end
test(vertex_list, edge_list, barabasi_albert(10,5))

println("These dynamics should flow to π, at t=500. they are there up to $(maximum(abs.(sol_st_ver(500.) .- pi)))")
@test maximum(abs.(sol_st_ver(500.) .- pi)) < 10^-7


dx = similar(x0)
e = [[ones(N)],[ones(N)]]
foo! = diff_network_st_ver.f
foo!(dx, x0, nothing, 0.)
@test (@allocated foo!(dx, x0, nothing, 0.)) == 0




### Inspect unrolled loop
# using Unrolled
# using .nd_ODE_Static_mod
# nd_ODE_Static_mod.vertex_loop!
# d = prob_st_ver.f.f
# @code_unrolled(vertex_loop!(x0, nothing, 0,
#     d.graph_data, d.graph_structure,
#     d.unique_vertices!, d.unique_v_indices, d.parallel))
