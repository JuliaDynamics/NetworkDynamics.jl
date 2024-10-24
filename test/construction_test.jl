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

@testset "graphless constructor" begin
    @test_throws ArgumentError ODEVertex(x->x^1, 2, 0; metadata="foba")
    v1 = ODEVertex(x->x^1, 2, 0; metadata=Dict(:graphelement=>1), name=:v1)
    @test has_graphelement(v1) && get_graphelement(v1) == 1
    v2 = ODEVertex(x->x^2, 2, 0; name=:v2)
    set_graphelement!(v2, 3)
    v3 = ODEVertex(x->x^3, 2, 0; name=:v3)
    set_graphelement!(v3, 2)

    e1 = StaticEdge(nothing, 0, Symmetric(); graphelement=(;src=1,dst=2))
    @test get_graphelement(e1) == (;src=1,dst=2)
    e2 = StaticEdge(nothing, 0, Symmetric())
    set_graphelement!(e2, (;src=:v3,dst=:v2))
    e3 = StaticEdge(nothing, 0, Symmetric())

    @test_throws ArgumentError Network([v1,v2,v3], [e1,e2,e3])
    set_graphelement!(e3, (;src=3,dst=1))

    nw = Network([v1,v2,v3], [e1,e2,e3])
    @test nw.im.vertexf == [v1, v3, v2]
    g = SimpleDiGraph(3)
    add_edge!(g, 1, 2)
    add_edge!(g, 2, 3)
    add_edge!(g, 3, 1)
    @test nw.im.g == g

    set_graphelement!(e3, (;src=1,dst=2))
    @test_throws ArgumentError Network([v1,v2,v3], [e1,e2,e3])

    set_graphelement!(e3, (;src=2,dst=1))
    Network([v1,v2,v3], [e1,e2,e3]) # throws waring about 1->2 and 2->1 beeing present

    v1 = ODEVertex(x->x^1, 2, 0; metadata=Dict(:graphelement=>1), name=:v1)
    v2 = ODEVertex(x->x^2, 2, 0; name=:v2, vidx=2)
    v3 = ODEVertex(x->x^3, 2, 0; name=:v3, vidx=3)
    nw = Network([v1,v2,v3], [e1,e2,e3])
    @test nw.im.unique_vnames == Dict(:v1=>1, :v2=>2, :v3=>3)

    v1 = ODEVertex(x->x^1, 2, 0; metadata=Dict(:graphelement=>1), name=:v2)
    v2 = ODEVertex(x->x^2, 2, 0; name=:v2, vidx=2)
    v3 = ODEVertex(x->x^3, 2, 0; name=:v3, vidx=3)
    set_graphelement!(e2, 3=>2)
    nw = Network([v1,v2,v3], [e1,e2,e3])
    @test nw.im.unique_vnames == Dict(:v3=>3)
end
@testset "Vertex batch" begin
    using NetworkDynamics: BatchStride, VertexBatch, parameter_range
    vb = VertexBatch{ODEVertex, typeof(sum), Vector{Int}}([1, 2, 3, 4], # vertices
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

@testset "eager gbuf map construction" begin
    using NetworkDynamics: gbuf_range
    e1 = StaticEdge(x->x^1, 1, 0, AntiSymmetric())
    e2 = StaticEdge(x->x^1, 1, 0, AntiSymmetric())
    v = ODEVertex(x->x^1, 1, 0)
    g = path_graph(4)
    nd = Network(g, v, [e1,e2,e1])
    map = nd.gbufprovider.map
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
    using NetworkDynamics: pdef, def
    v = ODEVertex(identity; dim=2, pdim=3)
    @test v.name == :ODEVertex
    @test v.obsf == nothing
    @test length(v.sym) == 2
    @test v.mass_matrix == LinearAlgebra.I

    v = ODEVertex(identity; sym=[:foo,:bar], pdim=3)
    @test dim(v) == 2
    @test length(v.psym) == 3

    v = ODEVertex(identity; sym=[:foo,:bar], psym=[:a])
    @test pdim(v)==1
    @test dim(v)==2

    @test_throws ArgumentError ODEVertex(identity; dim=1, pdim=1, mass_matrix=[1 2;3 4])

    @test_throws ArgumentError ODEVertex(identity, 1, 0; obsf=identity)
    v = ODEVertex(identity, 1, 0; obsf=identity, obssym=[:foo])

    @test_throws ArgumentError ODEVertex(identity, 1, 0; depth=2)

    StaticEdge(identity, 5, 0, Fiducial(); depth=2)
    StaticEdge(identity, 5, 0, Fiducial(); depth=1)
    @test_throws ArgumentError StaticEdge(identity, 5, 0, Fiducial(); depth=3)

    e = StaticEdge(identity, 5, 0, Directed())
    @test e.name == :StaticEdge
    e = ODEEdge(identity, 5, 0, Directed())
    StaticVertex(identity, 5, 0)

    v = ODEVertex(identity, 2, 3)
    @test dim(v) == 2
    @test pdim(v) == 3
    v = ODEVertex(identity, [:foo, :bar], 3)
    @test dim(v) == 2
    @test pdim(v) == 3
    v = ODEVertex(identity, 2, [:a, :b, :c])
    @test dim(v) == 2
    @test pdim(v) == 3
    v = ODEVertex(identity, [:foo, :bar], [:a, :b, :c])
    @test dim(v) == 2
    @test pdim(v) == 3

    using NetworkDynamics: _has_metadata
    @test _has_metadata([:a,:b,:c]) == false
    @test _has_metadata([:a=>1,:b,:c]) == true
    @test _has_metadata([:a=>1,:b=>2,:c=>3]) == true

    v = ODEVertex(identity, [:foo=>1, :bar], [:a=>2, :b, :c=>7])
    @test def(v) == [1,nothing]
    @test pdef(v) == [2,nothing,7]

    @test_throws ArgumentError ODEVertex(identity, [:foo=>1]; def=[1])
    @test_throws ArgumentError ODEVertex(identity, 1, [:foo=>1]; pdef=[1])

    v = ODEVertex(identity, :foo, :bar)
    @test v.sym == [:foo]
    @test v.psym == [:bar]

    v = ODEVertex(identity, :foo=>1, :bar)
    @test v.sym == [:foo]
    @test def(v) == [1]
    @test v.psym == [:bar]

    v = ODEVertex(identity, :foo, :bar=>1)
    @test v.sym == [:foo]
    @test v.psym == [:bar]
    @test pdef(v) == [1]
end

@testset "test dispatchT isbitstype" begin
    using InteractiveUtils
    using NetworkDynamics: EdgeFunction, VertexFunction, dispatchT, Symmetric

    for st in subtypes(EdgeFunction)
        st = st{Symmetric}
        T = Type{dispatchT(st)}
        @test Core.Compiler.isconstType(T)
    end
    for st in subtypes(VertexFunction)
        T = Type{dispatchT(st)}
        @test Core.Compiler.isconstType(T)
    end
end

@testset "test metadata constructors" begin
    using NetworkDynamics
    using NetworkDynamics: _fill_defaults, _split_metadata

    d1 = _split_metadata([:a=>1, :b=>2])[2]
    d2 = _split_metadata([:a=>(;default=1), :b=>2])[2]
    @test d1==d2

    d1 = _split_metadata([:a=>1, :b])[2]
    d2 = _split_metadata([:a=>(;default=1), :b])[2]
    @test d1==d2


    kwargs = Dict(:sym=>[:a,:b], :psym=>[:a,:d])
    @test_throws ArgumentError _fill_defaults(ODEVertex, kwargs)

    kwargs = Dict(:sym=>[:a,:b], :dim=>3)
    @test_throws ArgumentError _fill_defaults(ODEVertex, kwargs)

    kwargs = Dict(:sym=>[:a,:b], :dim=>3, :pdim=>0)
    @test_throws ArgumentError _fill_defaults(ODEVertex, kwargs)

    kwargs = Dict(:sym=>[:a,:b], :pdim=>0, :psym=>[:c])
    @test_throws ArgumentError _fill_defaults(ODEVertex, kwargs)

    # different length of sym/def
    kwargs = Dict(:sym=>[:a,:b],:def=>[1], :pdim=>0 )
    @test_throws ArgumentError _fill_defaults(ODEVertex, kwargs)

    # different provision of defs
    kwargs = Dict(:sym=>[:a,:b=>2],:def=>[1,nothing], :pdim=>0 )
    _fill_defaults(ODEVertex, kwargs)[:symmetadata]
    kwargs = Dict(:sym=>[:a=>2,:b],:def=>[1,nothing], :pdim=>0 )
    @test_throws ArgumentError _fill_defaults(ODEVertex, kwargs)
end

@testset "test dealias and copy of components" begin
    using NetworkDynamics: aliasgroups
    v1 = ODEVertex(x->x^1, 2, 0; metadata=Dict(:graphelement=>1), name=:v1)
    v2 = ODEVertex(x->x^2, 2, 0; name=:v2, vidx=2)
    v3 = ODEVertex(x->x^3, 2, 0; name=:v3, vidx=3)

    e1 = StaticEdge(nothing, 0, Symmetric(); graphelement=(;src=1,dst=2))
    e2 = StaticEdge(nothing, 0, Symmetric(); src=:v2, dst=:v3)
    e3 = StaticEdge(nothing, 0, Symmetric(); src=:v3, dst=:v1)

    @test isempty(aliasgroups([v1,v2,v3]))
    @test aliasgroups([v1, v1, v3]) == IdDict(v1 => [1,2])
    @test aliasgroups([v3, v1, v3]) == IdDict(v3 => [1,3])

    g = complete_graph(3)
    # with dealiasing / copying
    nw = Network(g, [v1,v1,v3],[e1,e2,e3]; dealias=true, check_graphelement=false)
    @test nw.im.vertexf[1].f == nw.im.vertexf[2].f
    @test nw.im.vertexf[1].metadata == nw.im.vertexf[2].metadata
    @test nw.im.vertexf[1].metadata !== nw.im.vertexf[2].metadata
    @test nw.im.vertexf[1].symmetadata == nw.im.vertexf[2].symmetadata
    @test nw.im.vertexf[1].symmetadata !== nw.im.vertexf[2].symmetadata
    s0 = NWState(nw)
    @test isnan(s0.v[1,1])
    set_default!(nw.im.vertexf[1], :v₁, 3)
    s1 = NWState(nw)
    @test s1.v[1,1] == 3

    # witout dealisasing
    nw = Network(g, [v1,v1,v3],[e1,e2,e1]; check_graphelement=false)
    @test nw.im.vertexf[1] === nw.im.vertexf[2]
    @test nw.im.edgef[1] === nw.im.edgef[3]
    @test keys(nw.im.aliased_vertexfs) == Set([v1])
    @test keys(nw.im.aliased_edgefs) == Set([e1])
    @test only(unique(values(nw.im.aliased_vertexfs))) == (;idxs=[1,2], hash=hash(v1))
    @test only(unique(values(nw.im.aliased_edgefs))) == (;idxs=[1,3], hash=hash(e1))
    s0 = NWState(nw)
    @test isnan(s0.v[1,1])
    prehash = hash(v1)
    set_default!(v1, :v₁, 3)
    posthash = hash(v1)
    @test prehash !== posthash
    @test NetworkDynamics.aliased_changed(nw; warn=false)
    s1 = NWState(nw)

    # test copy
    v = ODEVertex(x->x^1, 2, 0; metadata=Dict(:graphelement=>1), symmetadata=Dict(:x=>Dict(:default=>1)))
    v2 = copy(v)
    @test v !== v2
    @test v == v2
    @test isequal(v, v2)
    @test get_default(v, :x) == get_default(v2, :x)
    set_default!(v, :x, 99)
    @test get_default(v, :x) == 99
    @test get_default(v2, :x) == 1
    @test v != v2
    @test v1 != v2
end

@testset "test network-remake constructor" begin
    v1 = ODEVertex(x->x^1, 2, 0; metadata=Dict(:graphelement=>1), name=:v1)
    v2 = ODEVertex(x->x^2, 2, 0; name=:v2, vidx=2)
    v3 = ODEVertex(x->x^3, 2, 0; name=:v3, vidx=3)

    e1 = StaticEdge(nothing, 0, Symmetric(); graphelement=(;src=1,dst=2))
    e2 = StaticEdge(nothing, 0, Symmetric(); src=:v1, dst=:v3)
    e3 = StaticEdge(nothing, 0, Symmetric(); src=:v2, dst=:v3)

    g = complete_graph(3)
    nw = Network(g, [v1,v2,v3],[e1,e2,e3])
    nw2 = Network(nw)
    nw2 = Network(nw; g=path_graph(3), edgef=[e1, e3])

    for aggT in subtypes(NetworkDynamics.Aggregator)
        @show aggT
        @test hasmethod(NetworkDynamics.get_aggr_constructor, (aggT,))
    end
end
