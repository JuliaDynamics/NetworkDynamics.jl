using NetworkDynamics
using LightGraphs
using Test
using LinearAlgebra

N = 5
g = star_graph(N)
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
no_mm_edge  = ODEEdge(f! = (de, e, v_s, v_d, p, t) -> nothing, dim = 1)
num_mm_edge = ODEEdge(f! = (de, e, v_s, v_d, p, t) -> nothing, dim = 1,
                      mass_matrix = 0)
vec_mm_edge = ODEEdge(f! = (de, e, v_s, v_d, p, t) -> nothing, dim = 2,
                      mass_matrix = [1, 0])
mat_mm_edge = ODEEdge(f! = (de, e, v_s, v_d, p, t) -> nothing, dim = 2,
                      mass_matrix = [[1,1] [0,0]])

edge_list = [static_edge, no_mm_edge, num_mm_edge, vec_mm_edge, mat_mm_edge]

@test prod(network_dynamics(static_vertex, static_edge, g).mass_matrix .== zeros(N))

@test network_dynamics(no_mm_vertex, static_edge, g).mass_matrix == I

@test network_dynamics(no_mm_vertex, no_mm_edge, g).mass_matrix == I

MM = zeros(7,7)
MM[2,2] = 1
MM[4,4] = 1
MM[6,6] = 1
MM[7,6] = 1


@test prod(network_dynamics(vertex_list, static_edge, g).mass_matrix .== MM)

@test prod(network_dynamics(vertex_list, edge_list, g).mass_matrix .==
           [MM 0*similar(MM);0*similar(MM) MM])

@test prod(network_dynamics(static_vertex, edge_list, g).mass_matrix .==
           [zeros(N,N) zeros(N,7); zeros(7,N) MM])
