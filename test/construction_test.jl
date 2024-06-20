using NetworkDynamics
using Graphs

@testset "nd construction tests" begin
    using NetworkDynamics: StateType, statetype, isdense
    g = complete_graph(10)
    vertexf = ODEVertex(; f=x -> x^2, dim=1, pdim=2)
    @test statetype(vertexf) == NetworkDynamics.Dynamic()

    edgef = StaticEdge(; f=x -> x^2,
        dim=2, pdim=3,
        coupling=AntiSymmetric())
    @test statetype(edgef) == NetworkDynamics.Static()

    nd = Network(g, vertexf, edgef; verbose=false)

    @test statetype(only(nd.vertexbatches)) == NetworkDynamics.Dynamic()
    @test statetype(only(nd.layer.edgebatches)) == NetworkDynamics.Static()
    @test isdense(nd.im)
    @test nd.im.lastidx_dynamic == nv(g)
    @test nd.im.lastidx_static == nd.im.lastidx_dynamic + ne(g) * 2
    @test nd.vertexbatches isa Tuple
    @test nd.layer.edgebatches isa Tuple

    using NetworkDynamics: statetype
    g = complete_graph(10)
    vertexf = ODEVertex(; f=x -> x^2, dim=1, pdim=2)
    edgef = StaticEdge(; f=x -> x^2, dim=2, pdim=3, coupling=AntiSymmetric())

    using NetworkDynamics: SequentialExecution
    nd = Network(g, vertexf, edgef; verbose=false, execution=SequentialExecution{true}())

    nd = Network(g, vertexf, edgef; verbose=false,
        execution=SequentialExecution{false}())

    _du = rand(dim(nd))
    _u = rand(dim(nd))
    _p = rand(pdim(nd))
    @test_throws ArgumentError nd(rand(dim(nd)+1), _u, _p, 0.0)
    @test_throws ArgumentError nd(_du, rand(dim(nd)+1), _p, 0.0)
    @test_throws ArgumentError nd(_du, _u, rand(pdim(nd)+1), 0.0)
end

@testset "Vertex batch" begin
    using NetworkDynamics: BatchStride, VertexBatch, parameter_range
    vb = VertexBatch{ODEVertex, typeof(sum)}([1, 2, 3, 4], # vertices
        sum, # function
        BatchStride(1, 3),
        BatchStride(4, 2),
        BatchStride(0, 0))
    @test parameter_range(vb, 1) == 4:5
    @test parameter_range(vb, 2) == 6:7
    @test parameter_range(vb, 3) == 8:9
    @test parameter_range(vb, 4) == 10:11
end

@testset "massmatrix construction test" begin
    using LinearAlgebra: I, UniformScaling, Diagonal
    v1 = ODEVertex(x->x^1, 2, 0; mass_matrix=I)
    v2 = ODEVertex(x->x^2, 2, 0; mass_matrix=Diagonal([2,0]))
    v3 = ODEVertex(x->x^3, 2, 0; mass_matrix=[1 2;3 4])
    v4 = ODEVertex(x->x^4, 2, 0; mass_matrix=UniformScaling(0))
    v5 = ODEVertex(x->x^5, 2, 0; mass_matrix=I)
    e1 = ODEEdge(x->x^1, 2, 0, Fiducial(); mass_matrix=I)
    e2 = ODEEdge(x->x^2, 2, 0, AntiSymmetric(); mass_matrix=Diagonal([2,0]))
    e3 = ODEEdge(x->x^3, 2, 0, NetworkDynamics.Symmetric(); mass_matrix=[1 2;3 4])
    e4 = ODEEdge(x->x^3, 2, 0, Directed(); mass_matrix=UniformScaling(0))
    nd = Network(path_graph(5), [v1,v2,v3,v4,v5], [e1,e2,e3,e4])

    mm = Matrix(Diagonal([1,1,2,0,1,4,0,0,1,1,1,1,2,0,1,4,0,0]))
    mm[5,6] = 2
    mm[6,5] = 3
    mm[15,16] = 2
    mm[16,15] = 3

    @test nd.mass_matrix == mm

    nd = Network(path_graph(4), [v1,v2,v4,v5], [e1,e2,e4])
    @test nd.mass_matrix isa Diagonal

    nd = Network(path_graph(4), [v1,v1,v1,v1], [e1,e1,e1])
    @test nd.mass_matrix == I && nd.mass_matrix isa UniformScaling
end

@testset "gbuf map construction" begin
    using NetworkDynamics: gbuf_range
    e1 = StaticEdge(x->x^1, 1, 0, AntiSymmetric())
    e2 = StaticEdge(x->x^1, 1, 0, AntiSymmetric())
    v = ODEVertex(x->x^1, 1, 0)
    g = path_graph(4)
    nd = Network(g, v, [e1,e2,e1])
    map = nd.layer.gather_map
    for batch in nd.layer.edgebatches
        for batch_subi in 1:length(batch)
            eidx = batch.indices[batch_subi]
            src, dst = map[gbuf_range(batch, batch_subi), :]
            edge = collect(edges(g))[eidx]
            @test src == edge.src
            @test dst == edge.dst
        end
    end
end

@testset "test componen function constructors" begin
    using LinearAlgebra
    v = ODEVertex(identity; dim=2, pdim=3)
    @test v.name == :ODEVertex
    @test v.obsf == nothing
    @test length(v.sym) == 2
    @test v.mass_matrix == LinearAlgebra.I

    v = ODEVertex(identity; sym=[:foo,:bar], pdim=3)
    @test v.dim == 2
    @test length(v.psym) == 3

    v = ODEVertex(identity; sym=[:foo,:bar], psym=[:a])
    @test v.pdim==1
    @test v.dim==2

    @test_throws ArgumentError ODEVertex(identity; dim=1, pdim=1, mass_matrix=[1 2;3 4])

    @test_throws ArgumentError ODEVertex(identity, 1, 0; obsf=identity)
    v = ODEVertex(identity, 1, 0; obsf=identity, obssym=[:foo])
    @test_throws ArgumentError ODEVertex(identity, 1, 0; obsf=nothing, obssym=[:foo])

    @test_throws ArgumentError ODEVertex(identity, 1, 0; depth=2)

    StaticEdge(identity, 5, 0, Fiducial(); depth=2)
    StaticEdge(identity, 5, 0, Fiducial(); depth=1)
    @test_throws ArgumentError StaticEdge(identity, 5, 0, Fiducial(); depth=3)

    e = StaticEdge(identity, 5, 0, Directed())
    @test e.name == :StaticEdge
    e = ODEEdge(identity, 5, 0, Directed())
    StaticVertex(identity, 5, 0)

    v = ODEVertex(identity, 2, 3)
    @test v.dim == 2
    @test v.pdim == 3
    v = ODEVertex(identity, [:foo, :bar], 3)
    @test v.dim == 2
    @test v.pdim == 3
    v = ODEVertex(identity, 2, [:a, :b, :c])
    @test v.dim == 2
    @test v.pdim == 3
    v = ODEVertex(identity, [:foo, :bar], [:a, :b, :c])
    @test v.dim == 2
    @test v.pdim == 3

    using NetworkDynamics: _has_defaults
    @test _has_defaults([:a,:b,:c]) == false
    @test _has_defaults([:a=>1,:b,:c]) == true
    @test _has_defaults([:a=>1,:b=>2,:c=>3]) == true

    v = ODEVertex(identity, [:foo=>1, :bar], [:a=>2, :b, :c=>7])
    @test v.def == [1,nothing]
    @test v.pdef == [2,nothing,7]

    @test_throws ArgumentError ODEVertex(identity, [:foo=>1]; def=[1])
    @test_throws ArgumentError ODEVertex(identity, 1, [:foo=>1]; pdef=[1])

    v = ODEVertex(identity, :foo, :bar)
    @test v.sym == [:foo]
    @test v.psym == [:bar]

    v = ODEVertex(identity, :foo=>1, :bar)
    @test v.sym == [:foo]
    @test v.def == [1]
    @test v.psym == [:bar]

    v = ODEVertex(identity, :foo, :bar=>1)
    @test v.sym == [:foo]
    @test v.psym == [:bar]
    @test v.pdef == [1]
end
