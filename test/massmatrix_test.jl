using NetworkDynamics
using Graphs
using Test
using LinearAlgebra

N = 5
g = star_digraph(N)
add_edge!(g, 2, 3)

static_vertex = StaticVertex(; f=(x, edges, p, t) -> nothing, dim=1, pdim=0) |> ODEVertex
no_mm_vertex  = ODEVertex(; f=(dx, x, edges, p, t) -> nothing, dim=1, pdim=0)
num_mm_vertex = ODEVertex(; f=(dx, x, edges, p, t) -> nothing, dim=1, pdim=0,
    mass_matrix=0)
vec_mm_vertex = ODEVertex(; f=(dx, x, edges, p, t) -> nothing, dim=2, pdim=0,
    mass_matrix=Diagonal([1, 0]))
mat_mm_vertex = ODEVertex(; f=(dx, x, edges, p, t) -> nothing, dim=2, pdim=0,
    mass_matrix=[[1, 1] [0, 0]])

vertex_list = [static_vertex, no_mm_vertex, num_mm_vertex, vec_mm_vertex, mat_mm_vertex]

static_edge = StaticEdge(; f=(e, v_s, v_d, p, t) -> nothing, dim=1, pdim=0, coupling=Directed())
no_mm_edge = ODEEdge(; f=(de, e, v_s, v_d, p, t) -> nothing, dim=1, pdim=0, coupling=Directed())
num_mm_edge = ODEEdge(; f=(de, e, v_s, v_d, p, t) -> nothing, dim=1, pdim=0,
    mass_matrix=0, coupling=Directed())
vec_mm_edge = ODEEdge(; f=(de, e, v_s, v_d, p, t) -> nothing, dim=2, pdim=0,
    mass_matrix=Diagonal([1, 0]), coupling=Directed())
mat_mm_edge = ODEEdge(; f=(de, e, v_s, v_d, p, t) -> nothing, dim=2, pdim=0,
    mass_matrix=[[1, 1] [0, 0]], coupling=Directed())
zero_mm_edge = ODEEdge(; f=(de, e, v_s, v_d, p, t) -> nothing, dim=1, pdim=0, coupling=Directed(), mass_matrix=0)

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
