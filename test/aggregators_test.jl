using NetworkDynamics
using NetworkDynamics: aggregate!
using Graphs
using Random
using Chairmarks
using InteractiveUtils
using Test
using StableRNGs

(isinteractive() && @__MODULE__()==Main ? includet : include)("ComponentLibrary.jl")

@testset "compare different aggregators" begin
    vtypes = [ODEVertex(x->1x, 2, 0),
              ODEVertex(x->2x, 3, 0),
              ODEVertex(x->3x, 4, 0)]

    etypes = [StaticEdge(; f=x->1x, dim=2, pdim=0, coupling=Symmetric()),
              StaticEdge(; f=x->2x, dim=3, pdim=0, coupling=Symmetric()),
              StaticEdge(; f=x->3x, dim=4, pdim=0, coupling=Symmetric()),
              StaticEdge(; f=x->4x, dim=2, pdim=0, coupling=AntiSymmetric()),
              StaticEdge(; f=x->5x, dim=3, pdim=0, coupling=AntiSymmetric()),
              StaticEdge(; f=x->6x, dim=4, pdim=0, coupling=AntiSymmetric()),
              StaticEdge(; f=x->7x, dim=4, pdim=0, coupling=Fiducial()),
              StaticEdge(; f=x->8x, dim=7, pdim=0, coupling=Fiducial()),
              StaticEdge(; f=x->9x, dim=2, pdim=0, coupling=Directed()),
              StaticEdge(; f=x->10x, dim=3, pdim=0, coupling=Directed())]

    @test Set(typeof.(NetworkDynamics.coupling.(etypes))) == Set(subtypes(NetworkDynamics.Coupling))

    rng = StableRNG(1)
    g = watts_strogatz(10_000, 4, 0.8; rng, is_directed=true)
    nvec = copy.(rand(rng, vtypes, nv(g)))
    evec = copy.(rand(rng, etypes, ne(g)))

    basenw = Network(g, nvec, evec);
    states = rand(rng, basenw.im.lastidx_static)
    results = Vector{Float64}[]

    # @b AggregationMap($(basenw.im), $(basenw.layer.edgebatches))
    # @b SparseAggregator($(basenw.im), $(basenw.layer.edgebatches))

    for accT in subtypes(NetworkDynamics.Aggregator)
        aggregator = accT(+)
        nw = Network(g, nvec, evec; aggregator);
        aggbuf = rand(rng, nw.im.lastidx_aggr);
        b = @b aggregate!($nw.layer.aggregator, $aggbuf, $states)
        @info "Execute $accT" b
        if accT ∉ [KAAggregator, ThreadedAggregator]
            # @test b.allocs==0
        end
        push!(results, aggbuf)
    end
    # check if all results are equal
    for i in 2:length(results)
        issame = results[1] ≈ results[i]
        @test issame
        if !issame
            println("Error for $(subtypes(NetworkDynamics.Aggregator)[i])")
            println("extrema(Δ) = $(extrema(results[1] .- results[i]))")
            @info "compare" results[1] results[i] results[1]-results[i]
        end
    end
end


@testset "AggregationMap" begin
    g = path_graph(5)
    n = ODEVertex(sum, 1, 0)
    e = [StaticEdge(; f=x->1x, dim=2, pdim=0, coupling=Symmetric()),
        StaticEdge(; f=x->2x, dim=2, pdim=0, coupling=AntiSymmetric()),
        StaticEdge(; f=x->3x, dim=4, pdim=0, coupling=Fiducial()),
        StaticEdge(; f=x->4x, dim=2, pdim=0, coupling=Directed())]

    nw = Network(g, n, e; aggregator=KAAggregator(+))
    aggrmap = nw.layer.aggregator.m

    @test aggrmap.range == nv(g)+1 : nv(g) + 10
    @test aggrmap.symrange == nv(g)+1 : nv(g) + 10 - 2

    @test aggrmap.map    == [3, 4,  5,  6, 7, 8, 0, 0, 9, 10]
    @test aggrmap.symmap == [1, 2, -3, -4, 0, 0, 5, 6]

    buf = zeros(10)
    NetworkDynamics.aggregate!(nw.layer.aggregator, buf, collect(1:nw.im.lastidx_static))

    nw2 = Network(g, n, e; aggregator=NaiveAggregator(+));
    buf2 = zeros(10)
    NetworkDynamics.aggregate!(nw2.layer.aggregator, buf2, collect(1:nw.im.lastidx_static))
    @test maximum(abs.(buf2 .- buf)) == 0
end
