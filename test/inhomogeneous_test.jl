using Test
using NetworkDynamics
using Graphs
using LinearAlgebra
using SparseArrays
using OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqNonlinearSolve
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

statvertex = VertexModel(;g=(v, edges, p, t) -> v .= pi, outdim=1) |> ff_to_constraint
odevertex = VertexModel(; f=diffusion_vertex!, g=1, dim=1)

staticedge = EdgeModel(; g=AntiSymmetric(diffusion_edge!), outdim=1)
staticedge2 = EdgeModel(; g=AntiSymmetric(diffusion_edge2!), outdim=1)

vertex_list = [statvertex, odevertex]
append!(vertex_list, [odevertex for i in 1:N-2])
edge_list = EdgeModel[staticedge for e in edges(g)]
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


estates1 = sol_st_ver(sol_st_ver.t[end], idxs=EIndex(1:ne(g),:o))
estates2 = NWState(diff_network_st_ver, sol_st_ver.u[end]).e[1:ne(g),:o]
estates3 = NWState(sol_st_ver, sol_st_ver.u[end]).e[1:ne(g),:o]
@test estates1 == estates2 == estates3


####
#### Many distinct component models
####
# Above `MAX_UNROLLED_BATCHES` the batches live in a Vector instead of a Tuple and each
# batch call becomes a dynamic dispatch. That must not cost an allocation per batch:
# the rhs allocation count has to stay flat as the network grows.
let
    manyvertex(c) = VertexModel(; f=(dv, v, esum, p, t) -> (dv[1] = esum[1] - c*v[1]; nothing),
                                g=StateMask(1:1), dim=1, pdim=0)
    function manynetwork(N)
        g = watts_strogatz(N, 4, 0.0)  # ring lattice, β=0 => deterministic
        # a distinct closure per vertex => one batch per vertex
        Network(g, [manyvertex(i/N) for i in 1:N], staticedge;
                execution=SequentialExecution{false}())
    end
    function rhs_allocs(nw)
        x0 = rand(dim(nw)); dx = similar(x0)
        nw(dx, x0, nothing, 0.0)
        @allocations nw(dx, x0, nothing, 0.0)
    end

    small = manynetwork(2 * NetworkDynamics.MAX_UNROLLED_BATCHES)
    large = manynetwork(8 * NetworkDynamics.MAX_UNROLLED_BATCHES)
    @test length(small.vertexbatches) > NetworkDynamics.MAX_UNROLLED_BATCHES
    @test small.vertexbatches isa AbstractVector
    @test rhs_allocs(small) == rhs_allocs(large) ≤ 4
end
