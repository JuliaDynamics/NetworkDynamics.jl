using NDPrototype
using Test
using LightGraphs

using NDPrototype: VertexBatch, parameter_range

@testset "NDPrototype.jl" begin
    @test "constructor" begin
        g = complete_graph(10)
        vertexf = ODEVertex(; f=x->x^2, dim=1, pdim=2)
        edgef = StaticEdge(; f=x->x^2, dim=2, pdim=3)

        nd = networkdynamic(g, vertexf, edgef; verbose=false);

        g = SimpleGraph(4)
        add_edge!(g, 1, 2)
        add_edge!(g, 1, 2)
    end

    @test "Vertex batch" begin
        vb = VertexBatch([1,2,3,4], # vertices
                         sum, # function
                         3, # dimension
                         1, # first index in state vector
                         2, # p dim
                         4 # first index of p
                         )
        @test parameter_range(vb, 1) == 4:5
        @test parameter_range(vb, 2) == 6:7
        @test parameter_range(vb, 3) == 8:9
        @test parameter_range(vb, 4) == 10:11
    end
end

@testset "batch_identical" begin
    using NDPrototype: batch_identical
    v = :foo
    idx = [1,7,2,5]
    @test batch_identical(v, idx) == (:foo, idx)
    v = [:foo, :foo, :bar, :baz, :foo]
    idx = [5,4,3,2,1]
    @test batch_identical(v, idx) == ([:foo, :bar, :baz], [[5,4,1], [3], [2]])
end

@testset "test Index manager" begin
    using NDPrototype: _lastinrange, IndexManager, _lastindex_data, _lastindex_para, isdense
    im = IndexManager()
    @test isdense(im)

    @test _lastindex_data(im) == 0
    @test _lastindex_para(im) == 0
    im.v_data[1] = 1:6
    im.e_para[1] = 1:3
    @test isdense(im)
    @test _lastindex_data(im) == 6
    @test _lastindex_para(im) == 3
    im.e_data[2] = 8:8
    im.v_para[2] = 4:5
    @test !isdense(im)
    @test _lastindex_data(im) == 8
    @test _lastindex_para(im) == 5

    using NDPrototype: register_edges!, register_vertices!
    im = IndexManager()
    @test (1, 1) == register_edges!(im, [1,2,3,4], 2, 3)
    @test (9, 13) == register_edges!(im, [5,6], 1, 1)
    @test im.e_data.keys == [1,2,3,4,5,6]
    @test im.e_para.keys == [1,2,3,4,5,6]
    @test all(values(im.e_data) .== [1:2, 3:4, 5:6, 7:8, 9:9, 10:10])
    @test all(values(im.e_para) .== [1:3, 4:6, 7:9, 10:12, 13:13, 14:14])
    @test (11, 15) == register_vertices!(im, [1,2], 2, 3)
    @test (15, 21) == register_vertices!(im, [5,6], 1, 1)
    @test im.v_data.keys == [1,2,5,6]
    @test im.v_para.keys == [1,2,5,6]
    @test all(values(im.v_data) .== [11:12,13:14,15:15,16:16])
    @test all(values(im.v_para) .== [15:17,18:20,21:21,22:22])
    @test_throws ErrorException register_vertices!(im, [1,2], 1, 1)
    @test isdense(im)
end

@testset "greedy edge coloring" begin
    using NDPrototype: color_edges_greedy, isvalid
    for i in 1:20
        g = complete_graph(i)
        gc = color_edges_greedy(g)
        @test isvalid(gc)
    end
end

@testset "invert cdpath!" begin
    using NDPrototype: ColoredGraph, setcolor!, isvalid, invert_cdpath!, color
    g = ColoredGraph(complete_graph(5))
    setcolor!(g, 1=>2, 1)
    @test isvalid(g)
    invert_cdpath!(g, 1, 2, 1)
    @test isvalid(g)
    @test color(g, 1=>2) == 2
    g = ColoredGraph(complete_graph(5))
    setcolor!(g, 1=>2, 1)
    setcolor!(g, 2=>3, 2)
    setcolor!(g, 3=>4, 1)
    setcolor!(g, 4=>5, 2)
    @test isvalid(g)
    invert_cdpath!(g, 1, 2, 1)
    @test isvalid(g)
    @test color(g, 1=>2) == 2
    @test color(g, 2=>3) == 1
    @test color(g, 3=>4) == 2
    @test color(g, 4=>5) == 1
end
