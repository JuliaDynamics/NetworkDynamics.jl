using NetworkDynamics
using Graphs
using Test
using LinearAlgebra

N = 5
g = star_digraph(N)
add_edge!(g, 2, 3)

static_vertex = VertexFunction(; g=(x, edges, p, t) -> nothing, outdim=1) |> ff_to_constraint
no_mm_vertex  = VertexFunction(; f=(dx, x, edges, p, t) -> nothing, dim=1, g=1)
num_mm_vertex = VertexFunction(; f=(dx, x, edges, p, t) -> nothing, dim=1, g=1, mass_matrix=0)
vec_mm_vertex = VertexFunction(; f=(dx, x, edges, p, t) -> nothing, dim=2, g=1, mass_matrix=Diagonal([1, 0]))
mat_mm_vertex = VertexFunction(; f=(dx, x, edges, p, t) -> nothing, dim=2, g=1, mass_matrix=[[1, 1] [0, 0]])

vertex_list = [static_vertex, no_mm_vertex, num_mm_vertex, vec_mm_vertex, mat_mm_vertex]

static_edge = EdgeFunction(; g=Directed((e, v_s, v_d, p, t) -> nothing), outdim=1)
no_mm_edge = EdgeFunction(; f=(de, e, v_s, v_d, p, t) -> nothing, g=Directed(1), dim=1)
num_mm_edge = EdgeFunction(; f=(de, e, v_s, v_d, p, t) -> nothing, g=Directed(1), dim=1, mass_matrix=0)
vec_mm_edge = EdgeFunction(; f=(de, e, v_s, v_d, p, t) -> nothing, g=Directed(1), dim=2, mass_matrix=Diagonal([1, 0]))
mat_mm_edge = EdgeFunction(; f=(de, e, v_s, v_d, p, t) -> nothing, g=Directed(1), dim=2, mass_matrix=[[1, 1] [0, 0]])
zero_mm_edge = EdgeFunction(; f=(de, e, v_s, v_d, p, t) -> nothing, g=Directed(1), dim=1, mass_matrix=0)

edge_list = [static_edge, no_mm_edge, num_mm_edge, vec_mm_edge, mat_mm_edge]
ode_edge_list = [zero_mm_edge, no_mm_edge, num_mm_edge, vec_mm_edge, mat_mm_edge]

MM = zeros(7, 7)
MM[2, 2] = 1
MM[4, 4] = 1
MM[6, 6] = 1
MM[7, 6] = 1

@testset "Construction of mass matrix" begin

    @test all(Network(g, static_vertex, static_edge).mass_matrix .== zeros(N))
    @test Network(g, no_mm_vertex, static_edge).mass_matrix == I
    @test Network(g, no_mm_vertex, no_mm_edge).mass_matrix == I
    @test all(Network(g, vertex_list, static_edge).mass_matrix .== MM)

    @test all(Network(g, vertex_list, ode_edge_list).mass_matrix .==
               [MM zeros(size(MM)); zeros(size(MM)) MM])
    @test all(Network(g, static_vertex, ode_edge_list).mass_matrix .==
               [zeros(N, N) zeros(N, 7); zeros(7, N) MM])
end
