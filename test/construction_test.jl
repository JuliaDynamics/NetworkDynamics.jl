@testset "nd construction tests" begin
    using NDPrototype: StateType, statetype, isdense
    g = complete_graph(10)
    vertexf = ODEVertex(; f=x -> x^2, dim=1, pdim=2)
    @test statetype(vertexf) == NDPrototype.Dynamic()

    edgef = StaticEdge(; f=x -> x^2,
        dim=2, pdim=3,
        coupling=AntiSymmetric())
    @test statetype(edgef) == NDPrototype.Static()

    nd = Network(g, vertexf, edgef; verbose=false)

    @test statetype(only(nd.vertexbatches)) == NDPrototype.Dynamic()
    @test statetype(only(nd.layer.edgebatches)) == NDPrototype.Static()
    @test isdense(nd.im)
    @test nd.im.lastidx_dynamic == nv(g)
    @test nd.im.lastidx_static == nd.im.lastidx_dynamic + ne(g) * 2
    @test nd.vertexbatches isa Tuple
    @test nd.layer.edgebatches isa Tuple

    using NDPrototype: statetype
    g = complete_graph(10)
    vertexf = ODEVertex(; f=x -> x^2, dim=1, pdim=2)
    edgef = StaticEdge(; f=x -> x^2, dim=2, pdim=3, coupling=AntiSymmetric())

    using NDPrototype: SequentialExecution
    nd = Network(g, vertexf, edgef; verbose=false, execution=SequentialExecution{true}())

    nd = Network(g, vertexf, edgef; verbose=false,
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
    @test false
end
