using NDPrototype
using Test
using Graphs
using SafeTestsets

using NDPrototype: VertexBatch, parameter_range

(isinteractive() ? includet : include)("ComponentLibrary.jl")

@testset "Test Component Library" begin
    using NDPrototype: compf
    a = Lib.diffusion_edge()
    b = Lib.diffusion_edge()
    @test compf(a) == compf(b)
    a = Lib.diffusion_edge_closure()
    b = Lib.diffusion_edge_closure()
    @test compf(a) != compf(b)
    Lib.diffusion_edge_fid()
end

@testset "subscript" begin
    using NDPrototype: subscript
    @test subscript(10) == "₁₀"
    @test subscript(5) == "₅"
end

using NDPrototype: stylesymbolarray
syms = [:a,:b,:c]
defaults = [1, 2, nothing]
stylesymbolarray(syms,defaults, Dict(1=>:red,2=>:orange))
stylesymbolarray(syms,defaults, Dict(1=>:red,2=>:red))

@testset "NDPrototype.jl" begin
    @testset "constructor" begin
        using NDPrototype: StateType, statetype, isdense
        g = complete_graph(10)
        vertexf = ODEVertex(; f=x -> x^2, dim=1, pdim=2)
        @test statetype(vertexf) == NDPrototype.Dynamic()

        edgef = StaticEdge(; f=x -> x^2,
                           dim=2, pdim=3,
                           coupling=AntiSymmetric())
        @test statetype(edgef) == NDPrototype.Static()

        nd = Network(g, vertexf, edgef; verbose=true)

        @test statetype(only(nd.vertexbatches)) == NDPrototype.Dynamic()
        @test statetype(only(nd.layer.edgebatches)) == NDPrototype.Static()
        @test isdense(nd.im)
        @test nd.im.lastidx_dynamic == nv(g)
        @test nd.im.lastidx_static == nd.im.lastidx_dynamic + ne(g) * 2
        @test nd.vertexbatches isa Tuple
        @test nd.layer.edgebatches isa Tuple
    end

    @testset "constructor" begin
        using NDPrototype: statetype
        g = complete_graph(10)
        vertexf = ODEVertex(; f=x -> x^2, dim=1, pdim=2)
        edgef = StaticEdge(; f=x -> x^2, dim=2, pdim=3, coupling=AntiSymmetric())

        using NDPrototype: SequentialExecution
        nd = Network(g, vertexf, edgef; verbose=true, execution=SequentialExecution{true}())

        nd = Network(g, vertexf, edgef; verbose=true,
                     execution=SequentialExecution{false}())
    end

    @testset "Vertex batch" begin
        using NDPrototype: BatchStride, VertexBatch, parameter_range
        vb = VertexBatch{nothing, typeof(sum)}([1, 2, 3, 4], # vertices
            sum, # function
            BatchStride(1, 3),
            BatchStride(4, 2),
            BatchStride(0, 0))
        @test parameter_range(vb, 1) == 4:5
        @test parameter_range(vb, 2) == 6:7
        @test parameter_range(vb, 3) == 8:9
        @test parameter_range(vb, 4) == 10:11
    end
end

@testset "batch_identical" begin
    using NDPrototype: _batch_identical
    v = :foo
    idx = [1, 7, 2, 5]
    @test _batch_identical(v, idx) == [idx]
    v = [:foo, :foo, :bar, :baz, :foo]
    idx = [5, 4, 3, 2, 1]
    @test _batch_identical(v, idx) == [[5, 4, 1], [3], [2]]
end

@testset "greedy edge coloring" begin
    using NDPrototype: color_edges_greedy, isvalid
    for i in 1:20
        g = complete_graph(i)
        colors = color_edges_greedy(g)
        @test isvalid(g, colors)
    end
end

@safetestset "Aggregation Tests" begin include("aggregators_test.jl") end
