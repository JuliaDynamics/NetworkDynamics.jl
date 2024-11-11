using NetworkDynamics
using NetworkDynamics: Symmetric
using Graphs
using LinearAlgebra: LinearAlgebra
using Chairmarks: @b

@testset "graphless constructor" begin
    g = (out, in, p, t) -> nothing
    v1 = VertexFunction(; g, outdim=2, metadata=Dict(:graphelement=>1), name=:v1)
    @test has_graphelement(v1) && get_graphelement(v1) == 1
    v2 = VertexFunction(; g, outdim=2, name=:v2)
    set_graphelement!(v2, 3)
    v3 = VertexFunction(; g, outdim=2, name=:v3)
    set_graphelement!(v3, 2)

    ge = (out, src, dst, p, t) -> nothing
    e1 = EdgeFunction(; g=ge, outdim=1, graphelement=(;src=1,dst=2))
    @test get_graphelement(e1) == (;src=1,dst=2)
    e2 = EdgeFunction(; g=ge, outdim=1)
    set_graphelement!(e2, (;src=:v3,dst=:v2))
    e3 = EdgeFunction(; g=ge, outdim=1)

    # missing graphelement on e3
    @test_throws ArgumentError Network([v1,v2,v3], [e1,e2,e3])
    set_graphelement!(e3, (;src=3,dst=1))

    nw = Network([v1,v2,v3], [e1,e2,e3])
    @test nw.im.vertexf == [v1, v3, v2]
    graph = SimpleDiGraph(3)
    add_edge!(graph, 1, 2)
    add_edge!(graph, 2, 3)
    add_edge!(graph, 3, 1)
    @test nw.im.g == graph

    set_graphelement!(e3, (;src=1,dst=2))
    # same graphelement on multiple edges
    @test_throws ArgumentError Network([v1,v2,v3], [e1,e2,e3])

    set_graphelement!(e3, (;src=2,dst=1))
    Network([v1,v2,v3], [e1,e2,e3]) # throws waring about 1->2 and 2->1 beeing present

    v1 = VertexFunction(; g, outdim=2, metadata=Dict(:graphelement=>1), name=:v1)
    v2 = VertexFunction(; g, outdim=2, name=:v2, vidx=2)
    v3 = VertexFunction(; g, outdim=2, name=:v3, vidx=3)
    nw = Network([v1,v2,v3], [e1,e2,e3])
    @test nw.im.unique_vnames == Dict(:v1=>1, :v2=>2, :v3=>3)

    v1 = VertexFunction(; g, outdim=2, metadata=Dict(:graphelement=>1), name=:v2)
    v2 = VertexFunction(; g, outdim=2, name=:v2, vidx=2)
    v3 =VertexFunction(; g, outdim=2, name=:v3, vidx=3)
    set_graphelement!(e2, 3=>2)
    nw = Network([v1,v2,v3], [e1,e2,e3])
    @test nw.im.unique_vnames == Dict(:v3=>3)
end

@testset "massmatrix construction test" begin
    using LinearAlgebra: I, UniformScaling, Diagonal
    fv = (dv, v, in, p, t) -> nothing
    v1 = VertexFunction(; f=fv, g=1:2, dim=2, mass_matrix=I)
    v2 = VertexFunction(; f=fv, g=1:2, dim=2, mass_matrix=Diagonal([2,0]))
    v3 = VertexFunction(; f=fv, g=1:2, dim=2, mass_matrix=[1 2;3 4])
    v4 = VertexFunction(; f=fv, g=1:2, dim=2, mass_matrix=UniformScaling(0))
    v5 = VertexFunction(; f=fv, g=1:2, dim=2, mass_matrix=I)
    fe = (de, e, in1, in2, p, t) -> nothing
    e1 = EdgeFunction(; f=fe, g=Fiducial([1,2],[2,1]), dim=2, mass_matrix=I)
    e2 = EdgeFunction(; f=fe, g=AntiSymmetric(1:2), dim=2, mass_matrix=Diagonal([2,0]))
    e3 = EdgeFunction(; f=fe, g=Symmetric(1:2), dim=2, mass_matrix=[1 2;3 4])
    e4 = EdgeFunction(; f=fe, g=Directed(1:2), dim=2, mass_matrix=UniformScaling(0))
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
        @test NetworkDynamics.outsym(cf) == (; src=[:₋x, :₋y], dst=[:x, :y])

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
        # no metadata without src/dst
        @test_throws ArgumentError EdgeFunction(f=f, g=Symmetric(g_single_ff), dim=1, outsym=[:x=>1])
        cf = EdgeFunction(g=AntiSymmetric(g_single_pff), outdim=2)
        @test cf.outsym == (; src=[:₋o₁, :₋o₂], dst=[:o₁, :o₂])

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
    @test_throws ArgumentError _fill_defaults(VertexFunction, kwargs)

    kwargs = Dict(:sym=>[:a,:b], :dim=>3)
    @test_throws ArgumentError _fill_defaults(VertexFunction, kwargs)

    kwargs = Dict(:sym=>[:a,:b], :dim=>3, :pdim=>0)
    @test_throws ArgumentError _fill_defaults(VertexFunction, kwargs)

    kwargs = Dict(:sym=>[:a,:b], :pdim=>0, :psym=>[:c])
    @test_throws ArgumentError _fill_defaults(VertexFunction, kwargs)

    # different length of sym/def
    kwargs = Dict(:sym=>[:a,:b],:def=>[1], :pdim=>0 )
    @test_throws ArgumentError _fill_defaults(VertexFunction, kwargs)

    # different provision of defs
    kwargs = Dict(:sym=>[:a,:b=>2],:def=>[1,nothing], :pdim=>0, :g=>1:2, :outdim=>2, :f=>identity)
    _fill_defaults(VertexFunction, kwargs)[:symmetadata]
    kwargs = Dict(:sym=>[:a=>2,:b],:def=>[1,nothing], :pdim=>0 )
    @test_throws ArgumentError _fill_defaults(VertexFunction, kwargs)
end

@testset "test dealias and copy of components" begin
    using NetworkDynamics: aliasgroups
    f = (dv, v, ein, p, t) -> nothing
    g = (out, in, p, t) -> nothing
    ge = (out, src, dst, p, t) -> nothing

    v1 = VertexFunction(;f, g=1:1, dim=2, metadata=Dict(:graphelement=>1), name=:v1)
    v2 = VertexFunction(;f, g=1:1, dim=2, name=:v2, vidx=2)
    v3 = VertexFunction(;f, g=1:1, dim=2, name=:v3, vidx=3)
    e1 = EdgeFunction(; g=ge, outdim=1, graphelement=(;src=1,dst=2))
    e2 = EdgeFunction(; g=ge, outdim=1, src=:v2, dst=:v3)
    e3 = EdgeFunction(; g=ge, outdim=1, src=:v3, dst=:v1)

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
    g = (out, in, p, t) -> nothing
    v = VertexFunction(; g, outdim=2, metadata=Dict(:graphelement=>1), name=:v1, symmetadata=Dict(:x=>Dict(:default=>1)))
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
    g = (out, in, p, t) -> nothing
    ge = (out, src, dst, p, t) -> nothing
    v1 = VertexFunction(; g, outdim=2, metadata=Dict(:graphelement=>1), name=:v1)
    v2 = VertexFunction(; g, outdim=2, name=:v2, vidx=2)
    v3 = VertexFunction(; g, outdim=2, name=:v3, vidx=3)
    e1 = EdgeFunction(; g=ge, outdim=1, graphelement=(;src=1,dst=2))
    e2 = EdgeFunction(; g=ge, outdim=1, src=:v1, dst=:v3)
    e3 = EdgeFunction(; g=ge, outdim=1, src=:v2, dst=:v3)

    g = complete_graph(3)
    nw = Network(g, [v1,v2,v3],[e1,e2,e3])
    nw2 = Network(nw)
    nw2 = Network(nw; g=path_graph(3), edgef=[e1, e3])

    for aggT in subtypes(NetworkDynamics.Aggregator)
        @show aggT
        @test hasmethod(NetworkDynamics.get_aggr_constructor, (aggT,))
    end
end

@testset "test conversion of ff vertex to non ff vertex" begin
    using NetworkDynamics: fftype, compfg
    e = EdgeFunction(g=AntiSymmetric((e1, e2, v1, v2, p, t) -> nothing), outsym=[:i_r, :i_i])

    slackg = function(v, esum, (V,δ), t)
        v[1] = V*sin(δ)
        v[2] = V*cos(δ)
        nothing
    end
    v_pff = VertexFunction(g=slackg, outsym=[:u_r, :u_i], psym=[:V, :δ], insym=[:i_r, :i_i])
    @test fftype(v_pff) == PureFeedForward()

    swingf = function(dv, (δ, ω), (i_r, i_i), (V,Pm,M,D), t)
        Pel = V*sin(δ)*i_r + V*cos(δ)*i_i
        dδ = ω
        dω = 1/M*(Pm - D*ω + Pel)
        dv .= dδ, dω
        nothing
    end
    swingg = function(o, (δ, _), (i_r, i_i), (V,), t)
        o[1] = V*sin(δ)
        o[2] = V*cos(δ)
        nothing
    end
    v_ff = VertexFunction(f=swingf, g=swingg, sym=[:δ, :ω],
        outsym=[:u_r, :u_i], insym=[:i_r, :i_i], psym=[:V,:Pm,:M,:D])

    NetworkDynamics.CHECK_COMPONENT[] = true
    @test_throws ArgumentError Network(path_graph(2), v_pff, e)
    @test_throws ArgumentError Network(path_graph(2), v_ff, e)
    NetworkDynamics.CHECK_COMPONENT[] = false

    for v in [v_ff, v_pff]
        out = zeros(outdim(v))
        du = zeros(dim(v))
        u = rand(dim(v))
        p = rand(pdim(v))
        in = rand(indim(v))
        # compfg(v)(out, du, u, in, p, NaN)
        b = @b $(compfg(v))($out, $du, $u, $in, $p, NaN)
        @test b.allocs == 0

        v2 = ff_to_constraint(v)
        out2 = zeros(outdim(v2))
        du2 = zeros(dim(v2))
        u2 = vcat(u, out)
        # compfg(v2)(out2, du2, u2, in, p, NaN)
        b = @b $(compfg(v2))($out2, $du2, $u2, $in, $p, NaN)
        @test b.allocs == 0
        @test out ≈ out2
        @test du2 ≈ vcat(du, zeros(length(out)))
    end
end
