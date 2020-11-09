using Test
using LightGraphs
using NetworkDynamics

@testset "Test GraphData Accessors" begin
    g = SimpleGraph(5)
    add_edge!(g, (1,2))
    add_edge!(g, (1,4))
    add_edge!(g, (1,5))
    add_edge!(g, (2,3))
    add_edge!(g, (2,4))
    add_edge!(g, (2,5))
    add_edge!(g, (3,4))
    add_edge!(g, (3,5))
    v_dims = [2 for i in vertices(g)]
    e_dims = [2 for i in edges(g)]
    gs = GraphStruct(g, v_dims, e_dims, [:t], [:t])

    v_array = rand(sum(v_dims))
    e_array = rand(sum(e_dims))
    gd = GraphData(v_array, e_array, gs)

    @test get_vertex(gd, 1) == v_array[1:2]
    @test get_vertex(gd, 2) == v_array[3:4]
    @test get_vertex(gd, 3) == v_array[5:6]
    @test get_vertex(gd, 4) == v_array[7:8]
    @test get_vertex(gd, 5) == v_array[9:10]
    @test get_edge(gd, 1) == e_array[1:2]
    @test get_edge(gd, 2) == e_array[3:4]
    @test get_edge(gd, 3) == e_array[5:6]
    @test get_edge(gd, 4) == e_array[7:8]
    @test get_edge(gd, 5) == e_array[9:10]
    @test get_edge(gd, 6) == e_array[11:12]
    @test get_edge(gd, 7) == e_array[13:14]
    @test get_edge(gd, 8) == e_array[15:16]

    @test get_src_vertex(gd, 1) == v_array[1:2]
    @test get_dst_vertex(gd, 1) == v_array[3:4]
    @test get_src_vertex(gd, 8) == v_array[5:6]
    @test get_dst_vertex(gd, 8) == v_array[9:10]

    @test get_out_edges(gd, 1) == [e_array[1:2],e_array[3:4],e_array[5:6]]
    @test get_out_edges(gd, 3) == [e_array[13:14],e_array[15:16]]

    @test get_in_edges(gd, 1) == []
    @test get_in_edges(gd, 3) == [e_array[7:8]]
    @test get_in_edges(gd, 5) == [e_array[5:6],e_array[11:12],e_array[15:16]]
end
