using NDPrototype
using NDPrototype: aggregate!
using Graphs
using Random
using Chairmarks
using InteractiveUtils

(isinteractive() ? includet : include)("ComponentLibrary.jl")

@testset "compare different aggregators" begin
    vtypes = [ODEVertex(sum, 2, 0),
              ODEVertex(sum, 3, 0),
              ODEVertex(sum, 4, 0)]

    etypes = [StaticEdge(; f=sum, dim=2, pdim=0, coupling=Symmetric()),
              StaticEdge(; f=sum, dim=3, pdim=0, coupling=Symmetric()),
              StaticEdge(; f=sum, dim=4, pdim=0, coupling=Symmetric()),
              StaticEdge(; f=sum, dim=2, pdim=0, coupling=AntiSymmetric()),
              StaticEdge(; f=sum, dim=3, pdim=0, coupling=AntiSymmetric()),
              StaticEdge(; f=sum, dim=4, pdim=0, coupling=AntiSymmetric()),
              StaticEdge(; f=sum, dim=4, pdim=0, coupling=Fiducial()),
              StaticEdge(; f=sum, dim=7, pdim=0, coupling=Fiducial()),
              StaticEdge(; f=sum, dim=2, pdim=0, coupling=Directed()),
              StaticEdge(; f=sum, dim=3, pdim=0, coupling=Directed())]

    @test Set(typeof.(NDPrototype.coupling.(etypes))) == Set(subtypes(NDPrototype.Coupling))

    g = watts_strogatz(10_000, 4, 0.8; seed=1, is_directed=true)
    nvec = rand(vtypes, nv(g))
    evec = rand(etypes, ne(g))

    @time nw1 = Network(g, nvec, evec; accumulator=NaiveAggregator(+));
    @time nw2 = Network(g, nvec, evec; accumulator=NNlibScatter(+));
    @time nw3 = Network(g, nvec, evec; accumulator=KAAggregator(+));
    @time nw4 = Network(g, nvec, evec; accumulator=SequentialAggregator(+));
    @time nw5 = Network(g, nvec, evec; accumulator=PolyesterAggregator(+));

    nstates = nw1.im.lastidx_static
    nagg = nw1.im.lastidx_aggr

    states = rand(nstates)

    aggbuf1 = rand(nagg);
    aggbuf2 = rand(nagg);
    aggbuf3 = rand(nagg);
    aggbuf4 = rand(nagg);
    aggbuf5 = rand(nagg);
    b1 = @b aggregate!($nw1.layer.aggregator, $aggbuf1, $states)
    b2 = @b aggregate!($nw2.layer.aggregator, $aggbuf2, $states)
    b3 = @b aggregate!($nw3.layer.aggregator, $aggbuf3, $states)
    b4 = @b aggregate!($nw4.layer.aggregator, $aggbuf4, $states)
    b5 = @b aggregate!($nw5.layer.aggregator, $aggbuf5, $states)

    @test b1.allocs==0
    @test b2.allocs==0
    @test b4.allocs==0
    @test b5.allocs==0
    @test aggbuf1 ≈ aggbuf2 ≈ aggbuf3 ≈ aggbuf4 ≈ aggbuf5
end


@testset "AggregationMap" begin
    g = path_graph(5)
    n = ODEVertex(sum, 1, 0)
    e = [StaticEdge(; f=sum, dim=2, pdim=0, coupling=Symmetric()),
        StaticEdge(; f=sum, dim=2, pdim=0, coupling=AntiSymmetric()),
        StaticEdge(; f=sum, dim=4, pdim=0, coupling=Fiducial()),
        StaticEdge(; f=sum, dim=2, pdim=0, coupling=Directed())]

    nw = Network(g, n, e; accumulator=KAAggregator(+));
    aggrmap = nw.layer.aggregator.m

    @test aggrmap.range == nv(g)+1 : nv(g) + 10
    @test aggrmap.symrange == nv(g)+1 : nv(g) + 10 - 2

    @test aggrmap.map    == [3, 4,  5,  6, 7, 8, 0, 0, 9, 10]
    @test aggrmap.symmap == [1, 2, -3, -4, 0, 0, 5, 6]

    buf = zeros(10)
    NDPrototype.aggregate!(nw.layer.aggregator, buf, collect(1:nw.im.lastidx_static))

    nw2 = Network(g, n, e; accumulator=NaiveAggregator(+));
    buf2 = zeros(10)
    NDPrototype.aggregate!(nw2.layer.aggregator, buf2, collect(1:nw.im.lastidx_static))
    @test maximum(abs.(buf2 .- buf)) == 0
end
