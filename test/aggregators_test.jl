using NDPrototype
using NDPrototype: aggregate!
using Graphs
using Random
using Chairmarks

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


g = watts_strogatz(10_000, 4, 0.8; seed=1, is_directed=true)
nvec = rand(vtypes, nv(g))
evec = rand(etypes, ne(g))

@time nw1 = Network(g, nvec, evec; accumulator=NaiveAggregator(+));
@time nw2 = Network(g, nvec, evec; accumulator=NNlibScatter(+));

nstates = nw1.im.lastidx_static
nagg = nw1.im.lastidx_aggr

states = rand(nstates)
aggbuf1 = zeros(nagg)
aggbuf2 = copy(aggbuf1)

b1 = @b aggregate!($nw1.layer.aggregator, $aggbuf1, $states)
b2 = @b aggregate!($nw2.layer.aggregator, $aggbuf2, $states)
@test b1.allocs==0
@test b2.allocs==0
@test isapprox(aggbuf1, aggbuf2)
