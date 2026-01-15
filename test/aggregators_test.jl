using NetworkDynamics
using NetworkDynamics: aggregate!, NaiveAggregator
using Graphs
using Random
using Chairmarks
using InteractiveUtils
using Test
using StableRNGs
using ForwardDiff: Dual
using Symbolics

(isinteractive() && @__MODULE__()==Main ? includet : include)("ComponentLibrary.jl")

@testset "compare different aggregators" begin
    fv = (dv, v, ein, p, t) -> nothing
    vtypes = [VertexModel(; f=fv, dim=2, g=1:2),
              VertexModel(; f=fv, dim=3, g=1:2),
              VertexModel(; f=fv, dim=4, g=1:2)]

    gss = (e, vsrc, vdst, p, t) -> nothing
    etypes = [EdgeModel(; g=Symmetric(gss), outdim=2),
              EdgeModel(; g=Symmetric(gss), outdim=2),
              EdgeModel(; g=Symmetric(gss), outdim=2),
              EdgeModel(; g=AntiSymmetric(gss), outdim=2),
              EdgeModel(; g=AntiSymmetric(gss), outdim=2),
              EdgeModel(; g=AntiSymmetric(gss), outdim=2),
              EdgeModel(; g=Fiducial(gss,gss), outdim=2),
              EdgeModel(; g=Fiducial(gss,gss), outdim=2),
              EdgeModel(; g=Directed(gss), outdim=2),
              EdgeModel(; g=Directed(gss), outdim=2)]

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
        b = @b begin
            fill!($aggbuf, NetworkDynamics._appropriate_zero($aggbuf))
            aggregate!($nw.layer.aggregator, $aggbuf, $states)
        end
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

    # test that aggregate! is additive (doesn't overwrite, but adds to existing values)
    for accT in subtypes(NetworkDynamics.Aggregator)
        println("Testing additivity of $accT")
        aggregator = accT(+)
        nw = Network(g, nvec, evec; aggregator);
        aggbuf_zero = zeros(nw.im.lastidx_aggr)
        aggbuf_offset = collect(Float64, 1:nw.im.lastidx_aggr)
        aggregate!(nw.layer.aggregator, aggbuf_zero, states)
        aggregate!(nw.layer.aggregator, aggbuf_offset, states)
        @test aggbuf_offset - aggbuf_zero ≈ 1:nw.im.lastidx_aggr
    end
end

@testset "Test _appropriate_zero" begin
    @test NetworkDynamics._appropriate_zero([1,2,3]) isa Int
    @test NetworkDynamics._appropriate_zero([1,2,3]) == 0
    @test NetworkDynamics._appropriate_zero([1.0,2,3]) isa Float64
    @test NetworkDynamics._appropriate_zero([1.0,2,3]) == 0.0
    @test NetworkDynamics._appropriate_zero([Dual(1.0), 2, 3]) == Dual(0.0)
    @variables x, y, z
    @test NetworkDynamics._appropriate_zero([x,y,z]) isa Num
    @test NetworkDynamics._appropriate_zero(Any[x,y,z]) isa Float64
end
