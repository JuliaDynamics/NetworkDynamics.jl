using Test
using NetworkDynamics
using Graphs
using LinearAlgebra
using SparseArrays
using OrdinaryDiffEq
using Random
Random.seed!(42)

# We will work on a random graph:
g = barabasi_albert(10, 5)
L = laplacian_matrix(g)

N = nv(g)

printstyled("--- Network dynamics with a static vertex ---\n"; bold=true, color=:white)


@inline function diffusion_edge!(e, v_s, v_d, p, t)
    e .= v_s .- v_d
    nothing
end
@inline function diffusion_edge2!(e, v_s, v_d, p, t)
    e .= v_s .- v_d
    nothing
end
@inline function diffusion_vertex!(dv, v, esum, p, t)
    dv .= esum
    nothing
end

statvertex = StaticVertex((v, edges, p, t) -> v .= pi; dim=1) |> ODEVertex
odevertex = ODEVertex(; f=diffusion_vertex!, dim=1)

staticedge = StaticEdge(; f=diffusion_edge!, dim=1, coupling=AntiSymmetric())
staticedge2 = StaticEdge(; f=diffusion_edge2!, dim=1, coupling=AntiSymmetric())

vertex_list = [statvertex, odevertex]
append!(vertex_list, [odevertex for i in 1:N-2])
edge_list = StaticEdge[staticedge for e in edges(g)]
edge_list[2] = staticedge2

diff_network_st_ver = Network(g, vertex_list, edge_list)

x0 = rand(nv(g))
NWState(diff_network_st_ver, x0).v[1,1] = pi

prob_st_ver = ODEProblem(diff_network_st_ver, x0, (0.,500.))
Main.test_execution_styles(prob_st_ver) # testing all ex styles #src
sol_st_ver = solve(prob_st_ver, Rodas4());

println("These dynamics should flow to π, at t=500. they are there up to $(maximum(abs.(sol_st_ver(500.) .- pi)))")
@test maximum(abs.(sol_st_ver(500.0) .- pi)) < 10^-7

dx = similar(x0)
f = diff_network_st_ver
f(dx, x0, nothing, 0.0)
@test (@allocations f(dx, x0, nothing, 0.0)) == 0


estates1 = sol_st_ver(sol_st_ver.t[end], idxs=EIndex(1:ne(g),1))
estates2 = NWState(diff_network_st_ver, sol_st_ver[end]).e[1:ne(g),1]
estates3 = NWState(sol_st_ver, sol_st_ver[end]).e[1:ne(g),1]
@test estates1 == estates2 == estates3
