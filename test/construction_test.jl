using NetworkDynamics
using NetworkDynamics: Symmetric
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

@testset "test component function constructors" begin
    using NetworkDynamics: PureFeedForward, FeedForward, NoFeedForward, PureStateMap
    using NetworkDynamics: fftype
    using NetworkDynamics: _has_metadata
    using NetworkDynamics: get_defaults, outsym
    @testset "VertexFunction construction" begin
        f = (du, u, in, p, t) -> nothing
        g_pff = (o, in, p, t) -> nothing
        g_ff  = (o, u, in, p, t) -> nothing
        g_nff = (o, u, p, t) -> nothing
        g_psm = (o, u) -> nothing

        cf = VertexFunction(g=g_pff, outdim=1, insym=[:a])
        @test fftype(cf) == PureFeedForward()

        cf = VertexFunction(f=f, dim=1, g=g_nff, outdim=1, insym=[:a])
        @test fftype(cf) == NoFeedForward()

        @test_throws ArgumentError VertexFunction(f=f, sym=[:x,:y], g=StateMask([2,1]), outdim=1, insym=[:a])
        cf = VertexFunction(f=f, sym=[:x,:y], g=StateMask([2,1]), insym=[:a])
        @test fftype(cf) == PureStateMap()
        @test NetworkDynamics.outsym(cf) == [:y, :x]

        cf = VertexFunction(f=f, sym=[:x,:y], g=g_ff, outdim=1, insym=[:a], name=:foo)
        @test fftype(cf) == FeedForward()
        @test NetworkDynamics.outdim(cf) == 1
        @test cf.name == :foo
        @test cf.obsf == nothing
        @test cf.mass_matrix == LinearAlgebra.I
        @test pdim(cf) == 0

        cf = VertexFunction(f=f, g=g_ff, dim=1, outdim=1, pdim=3, indim=2, name=:foo)
        @test length(cf.insym) == 2
        @test length(cf.psym) == 3

        # mismatch of massmatrix size
        @test_throws ArgumentError VertexFunction(f=identity, g=g_psm, dim=1, pdim=1, outdim=1, mass_matrix=[1 2;3 4])

        # passing of obs function
        @test_throws ArgumentError VertexFunction(g=g_pff, outdim=1, obsf=identity)
        cf = VertexFunction(g=g_pff, outdim=1, obsf=identity, obssym=[:foo, :bar])
        obssym(cf) == [:foo, :bar]

        # pass graph element
        cf = VertexFunction(g=g_pff, outdim=1, obsf=identity, obssym=[:foo, :bar], vidx=7)
        @test cf.metadata[:graphelement] == 7
        cf = VertexFunction(g=g_pff, outdim=1, obsf=identity, obssym=[:foo, :bar], graphelement=2)
        @test cf.metadata[:graphelement] == 2

        # passing of metadata
        @test _has_metadata([:a,:b,:c]) == false
        @test _has_metadata([:a=>1,:b,:c]) == true
        @test _has_metadata([:a=>1,:b=>2,:c=>3]) == true

        cf = VertexFunction(g=g_pff, outsym=[:foo=>1, :bar], psym=[:a=>2, :b, :c=>7])
        @test get_default(cf, :foo) == 1
        @test get_defaults(cf, outsym(cf)) == [1, nothing]
        @test get_defaults(cf, psym(cf)) == [2, nothing, 7]

        @test_throws ArgumentError VertexFunction(f=f, g=StateMask(1:1), sym=[:foo=>1]; def=[1])
        @test_throws ArgumentError VertexFunction(f=f, g=StateMask(2:2), dim=2, psym=[:foo=>1]; pdef=[1])
    end


    @testset "EdgeFunction Constructor" begin
        using NetworkDynamics
        using NetworkDynamics: Symmetric
        f = (du, u, in1, in2, p, t) -> nothing
        g_pff = (o1, o2, in1, in2, p, t) -> nothing
        g_ff  = (o1, o2, u, in1, in2, p, t) -> nothing
        g_nff = (o1, o2, u, p, t) -> nothing
        g_psm = (o1, o2, u) -> nothing

        g_single_pff = (o2, in1, in2, p, t) -> nothing
        g_single_ff = (o2, u, in1, in2, p, t) -> nothing
        g_single_nff = (o2, u, p, t) -> nothing
        g_single_psm = (o2, u) -> nothing
        for c in (AntiSymmetric, Symmetric, Directed)
            @test fftype(c(g_single_pff)) == PureFeedForward()
            @test fftype(c(g_single_ff)) == FeedForward()
            @test fftype(c(g_single_nff)) == NoFeedForward()
            @test fftype(c(g_single_psm)) == PureStateMap()
            @test fftype(c(1:2)) == PureStateMap()
        end

        cf = EdgeFunction(g=g_pff, outdim=1, insym=[:a])
        @test fftype(cf) == PureFeedForward()

        cf = EdgeFunction(f=f, dim=1, g=g_nff, outdim=1, insym=[:a])
        @test fftype(cf) == NoFeedForward()

        @test_throws ArgumentError EdgeFunction(f=f, sym=[:x,:y], g=StateMask([2,1]), outdim=1, insym=[:a])
        @test_throws ArgumentError EdgeFunction(f=f, sym=[:x,:y], g=StateMask([2,1]), insym=[:a])

        cf = EdgeFunction(f=f, sym=[:x,:y], g=AntiSymmetric(StateMask(1:2)))
        @test fftype(cf) == PureStateMap()
        @test NetworkDynamics.outsym(cf) == (; src=[:src₊x, :src₊y], dst=[:dst₊x, :dst₊y])

        cf = EdgeFunction(f=f, sym=[:x,:y], g=g_ff, outdim=1, insym=[:a], name=:foo)
        @test fftype(cf) == FeedForward()
        @test NetworkDynamics.outdim(cf) == (;src=1, dst=1)
        @test cf.name == :foo
        @test cf.obsf == nothing
        @test cf.mass_matrix == LinearAlgebra.I
        @test pdim(cf) == 0


        # output sym generation
        cf = EdgeFunction(f=f, sym=[:x,:y], g=AntiSymmetric(StateMask(1:2)))
        @test cf.outsym == (; src=[:₋x,:₋y], dst=[:x,:y])
        cf = EdgeFunction(f=f, sym=[:x,:y], g=Fiducial(1:2, 2:-1:1))
        @test cf.outsym == (; src=[:x,:y], dst=[:y,:x])
        cf = EdgeFunction(f=f, g=g_ff, dim=1, outsym=(;src=[:a=>2], dst=[:b=>4]))
        @test cf.outsym == (; src=[:a], dst=[:b])
        cf = EdgeFunction(f=f, g=g_ff, dim=1, outdim=2)
        @test cf.outsym == (; src=[:src₊o₁, :src₊o₂], dst=[:dst₊o₁, :dst₊o₂])
        cf = EdgeFunction(f=f, g=Directed(g_single_ff), dim=1, outdim=2)
        @test cf.outsym == (; src=[], dst=[:o₁, :o₂])
        cf = EdgeFunction(f=f, g=Directed(1:2), dim=2, outdim=2)
        @test isempty(cf.outsym.src)
        @test cf.outsym.dst == cf.sym
        @test_throws ArgumentError EdgeFunction(f=f, g=Symmetric(g_single_ff), dim=1, outsym=[:x=>1])

        # input sym generation
        cf = EdgeFunction(f=f, g=g_ff, dim=1, outdim=1, pdim=3, indim=2, name=:foo)
        @test length(cf.insym) == 2
        @test length(cf.psym) == 3
        @test cf.insym == (;src = [:src₊i₁, :src₊i₂], dst = [:dst₊i₁, :dst₊i₂])

        cf = EdgeFunction(f=f, g=g_ff, dim=1, outdim=1, insym=[:foo])
        @test_throws ArgumentError EdgeFunction(f=f, g=g_ff, dim=1, outdim=1, insym=(src=[:foo], dst=[:foo]))
    end
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
