using NetworkDynamics
using NetworkDynamics: aggregate!, NaiveAggregator
using Graphs
using Random
using Chairmarks
using InteractiveUtils
using Test
using StableRNGs

(isinteractive() && @__MODULE__()==Main ? includet : include)("ComponentLibrary.jl")

@testset "compare different aggregators" begin
    fv = (dv, v, ein, p, t) -> nothing
    vtypes = [VertexFunction(; f=fv, dim=2, g=1:2),
              VertexFunction(; f=fv, dim=3, g=1:2),
              VertexFunction(; f=fv, dim=4, g=1:2)]

    gss = (e, vsrc, vdst, p, t) -> nothing
    etypes = [EdgeFunction(; g=Symmetric(gss), outdim=2),
              EdgeFunction(; g=Symmetric(gss), outdim=2),
              EdgeFunction(; g=Symmetric(gss), outdim=2),
              EdgeFunction(; g=AntiSymmetric(gss), outdim=2),
              EdgeFunction(; g=AntiSymmetric(gss), outdim=2),
              EdgeFunction(; g=AntiSymmetric(gss), outdim=2),
              EdgeFunction(; g=Fiducial(gss,gss), outdim=2),
              EdgeFunction(; g=Fiducial(gss,gss), outdim=2),
              EdgeFunction(; g=Directed(gss), outdim=2),
              EdgeFunction(; g=Directed(gss), outdim=2)]

    rng = StableRNG(1)
    g = watts_strogatz(10_000, 4, 0.8; rng, is_directed=true)
    nvec = copy.(rand(rng, vtypes, nv(g)))
    evec = copy.(rand(rng, etypes, ne(g)))

    basenw = Network(g, nvec, evec);
    states = rand(rng, basenw.im.lastidx_out)

    # @b AggregationMap($(basenw.im), $(basenw.layer.edgebatches))
    # @b SparseAggregator($(basenw.im), $(basenw.layer.edgebatches))

    results = Vector{Float64}[]
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
        if !issame
            println("Error for $(subtypes(NetworkDynamics.Aggregator)[i])")
            println("extrema(Δ) = $(extrema(results[1] .- results[i]))")
            @info "compare" results[1] results[i] results[1]-results[i]
        end
        @test issame
    end
end
