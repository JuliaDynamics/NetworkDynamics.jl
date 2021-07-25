using NDPrototype
using Test
using LightGraphs

using NDPrototype: VertexBatch, parameter_range

@testset "NDPrototype.jl" begin
    @test "network layer" begin

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
