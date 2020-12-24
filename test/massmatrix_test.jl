using NetworkDynamics
using LightGraphs
using Test
using LinearAlgebra

N = 5
g = star_digraph(N)
add_edge!(g, 2, 3)

static_vertex = StaticVertex(f! = (x, e_s, e_d, p, t) -> nothing, dim = 1)
no_mm_vertex  = ODEVertex(f! = (dx, x, e_s, e_d, p, t) -> nothing, dim = 1)
num_mm_vertex = ODEVertex(f! = (dx, x, e_s, e_d, p, t) -> nothing, dim = 1,
                          mass_matrix = 0)
vec_mm_vertex = ODEVertex(f! = (dx, x, e_s, e_d, p, t) -> nothing, dim = 2,
                          mass_matrix = [1, 0])
mat_mm_vertex = ODEVertex(f! = (dx, x, e_s, e_d, p, t) -> nothing, dim = 2,
                          mass_matrix = [[1,1] [0,0]])

vertex_list = [static_vertex, no_mm_vertex, num_mm_vertex, vec_mm_vertex, mat_mm_vertex]

static_edge = StaticEdge(f! = (e, v_s, v_d, p, t) -> nothing, dim = 1)
no_mm_edge  = ODEEdge(f! = (de, e, v_s, v_d, p, t) -> nothing, dim = 1, coupling= :directed)
num_mm_edge = ODEEdge(f! = (de, e, v_s, v_d, p, t) -> nothing, dim = 1,
                      mass_matrix = 0, coupling= :directed)
vec_mm_edge = ODEEdge(f! = (de, e, v_s, v_d, p, t) -> nothing, dim = 2,
                      mass_matrix = [1, 0], coupling= :directed)
mat_mm_edge = ODEEdge(f! = (de, e, v_s, v_d, p, t) -> nothing, dim = 2,
                      mass_matrix = [[1,1] [0,0]], coupling= :directed)
zero_mm_edge  = ODEEdge(f! = (de, e, v_s, v_d, p, t) -> nothing, dim = 1, coupling= :directed, mass_matrix = 0)

edge_list = [static_edge, no_mm_edge, num_mm_edge, vec_mm_edge, mat_mm_edge]
ode_edge_list = [zero_mm_edge, no_mm_edge, num_mm_edge, vec_mm_edge, mat_mm_edge]

MM = zeros(7,7)
MM[2,2] = 1
MM[4,4] = 1
MM[6,6] = 1
MM[7,6] = 1

@testset "Construction of mass matrix" begin

        @test prod(network_dynamics(static_vertex, static_edge, g).mass_matrix .== zeros(N))
        @test network_dynamics(no_mm_vertex, static_edge, g).mass_matrix == I
        @test network_dynamics(no_mm_vertex, no_mm_edge, g).mass_matrix == I
        @test prod(network_dynamics(vertex_list, static_edge, g).mass_matrix .== MM)

        @test prod(network_dynamics(vertex_list, ode_edge_list, g).mass_matrix .==
                   [MM zeros(size(MM));zeros(size(MM)) MM])
        @test prod(network_dynamics(static_vertex, ode_edge_list, g).mass_matrix .==
                   [zeros(N,N) zeros(N,7); zeros(7,N) MM])
end


# This test is broken since static_edges with undefined coupling cant be promoted to ODEdges
@test_broken prod(network_dynamics(vertex_list, edge_list, g).mass_matrix .==
           [MM zeros(size(MM));zeros(size(MM)) MM])

# Testing for coupling type errors, this should go in another file
@test_throws ArgumentError network_dynamics(no_mm_vertex, no_mm_edge, SimpleGraph(g))
