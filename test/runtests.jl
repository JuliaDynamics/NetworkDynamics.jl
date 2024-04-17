using NDPrototype
using Test
using Graphs

using NDPrototype: VertexBatch, parameter_range

@testset "NDPrototype.jl" begin
    @testset "constructor" begin
        using NDPrototype: StateType, statetype
        g = complete_graph(10)
        vertexf = ODEVertex(; f=x->x^2, dim=1, pdim=2)
        @test statetype(vertexf) == StateType.dynamic

        edgef = StaticEdge(; f=x->x^2, dim=2, pdim=3, coupling=AntiSymmetric())
        statetype(edgef) == StateType.static

        nd = Network(g, vertexf, edgef; verbose=true);

        @test statetype(only(nd.vertexbatches)) == StateType.dynamic
        @test statetype(only(nd.nl.edgebatches)) == StateType.static
        @test isdense(nd.im)
        @test nd.im.size_dynamic == nv(g)
        @test nd.im.size_static == ne(g)*2
        nd.cachepool
        lbc = LazyBufferCache()
        buff = lbc[zeros(3)]
        buff[1] = 4
        buff = lbc[ones(3), 10]
        buff = lbc[ones(3)]
    end

    @testset "Vertex batch" begin
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
    using NDPrototype: _batch_identical
    v = :foo
    idx = [1,7,2,5]
    @test _batch_identical(v, idx) == ([:foo], [idx])
    v = [:foo, :foo, :bar, :baz, :foo]
    idx = [5,4,3,2,1]
    @test _batch_identical(v, idx) == ([:foo, :bar, :baz], [[5,4,1], [3], [2]])
end

@testset "greedy edge coloring" begin
    using NDPrototype: color_edges_greedy, isvalid
    for i in 1:20
        g = complete_graph(i)
        colors = color_edges_greedy(g)
        @test isvalid(g, colors)
    end
end
