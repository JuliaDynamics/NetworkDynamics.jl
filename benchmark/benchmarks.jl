using Pkg

using BenchmarkTools
using Graphs
using NetworkDynamics
using Serialization
using StableRNGs
using SciMLBase
using Test
using ProgressMeter
using Random
(isinteractive() ? includet : include)("benchmark_utils.jl")
(isinteractive() ? includet : include)("benchmark_compat.jl")

if pkgversion(NetworkDynamics) < v"0.9.0"
    (isinteractive() ? includet : include)("benchmark_models_v0.8.jl")
else
    (isinteractive() ? includet : include)("benchmark_models.jl")
end

@info "Benchmark with $(Threads.nthreads()) threads"

bd = BenchmarkDict()

executions = Dict("seq_buf" => SequentialExecution{true}(),
                  "ka_buf" => KAExecution{true}())
aggregations = Dict("nnlib" => NNlibScatter(+),
                    "KA" => KAAggregator(+),
                    "seq" => SequentialAggregator(+),
                    "poly" => PolyesterAggregator(+))

####
#### diffusion benchmarks on a dense watts strogatz graph
####
vertex = diffusion_vertex()
edges = Dict("static_edge" => diffusion_edge(),
             "ode_edge" => diffusion_dedge())

for k in ["static_edge", "ode_edge"]
    edge = edges[k]

    for N in [100, 300, 1000, 3000]  #, 100_000, 1_000_000]
        @info "Benchmark diffusion network with $k on size $N"
        g = watts_strogatz(N, N ÷ 2, 0.0; rng=StableRNG(1))
        b = @be Network($g, $vertex, $edge) evals=1 samples=50 seconds=10
        bd["diffusion", k, "assemble", N] = BenchmarkResult(b)

        _nd = Network(g, vertex, edge)
        _x0 = randx0(_nd)
        _dx = similar(_x0)
        _nd(_dx, _x0, nothing, NaN)

        progress = Progress(length(executions) * length(aggregations))
        for (exname, execution) in executions
            for (aggname, aggregator) in aggregations
                next!(progress)
                nd = try
                    Network(g, vertex, edge; execution, aggregator)
                catch e
                    if e.msg == "execution type not supported"
                        continue
                    else
                        rethrow(e)
                    end
                end
                dx = similar(_x0)
                nd(dx, _x0, nothing, NaN)
                @test dx ≈ _dx

                b = @be $nd($dx, $_x0, nothing, 0.0) seconds=1
                br = BenchmarkResult(b, legacy_order(nd, dx))
                bd["diffusion", k, exname, aggname, N] = br
            end
        end
        finish!(progress)
    end
end

####
#### kuramoto benchmarks on a sparse graph
####

function homogeneous(N)
    g = watts_strogatz(N, 3, 0.8; rng=StableRNG(1))
    edge = static_kuramoto_edge()
    vertex = kuramoto_vertex_2d()
    (vertex, edge, g)
end

function heterogeneous(N)
    rng = StableRNG(1)
    g = watts_strogatz(N, 3, 0.8; rng=StableRNG(1))
    edge = static_kuramoto_edge()
    vertex = [kuramoto_vertex_1d(), kuramoto_vertex_2d()]
    vertices = vertex[shuffle(rng, vcat([1 for _ in 1:N÷2], [2 for _ in 1:N÷2]))]
    (vertices, edge, g)
end

for f in [homogeneous, heterogeneous]
    name = string(f)
    for N in [100, 1_000, 10_000, 100_000]#, 1_000_000]
        @info "Benchmark kuramoto $name on size $N"
        (vert, edg, g) = f(N)
        b = @be Network($g, $vert, $edg) evals=1 samples=50 seconds=10
        bd["kuramoto", name, "assemble", N] = BenchmarkResult(b)

        _nd = Network(g, vert, edg)
        _x0 = randx0(_nd)
        _dx = similar(_x0)
        p = randp(_nd)
        _nd(_dx, _x0, p, NaN)

        progress = Progress(length(executions) * length(aggregations))
        for (exname, execution) in executions
            for (aggname, aggregator) in aggregations
                next!(progress)
                nd = try
                    Network(g, vert, edg; execution, aggregator)
                catch e
                    if e.msg == "execution type not supported"
                        continue
                    else
                        rethrow(e)
                    end
                end
                dx = similar(_x0)
                nd(dx, _x0, p, NaN)
                @test dx ≈ _dx

                b = @be $nd($dx, $_x0, $p, 0.0) seconds=1
                br = BenchmarkResult(b, legacy_order(nd, dx))
                bd["kuramoto", name, exname, aggname, N] = br
            end
        end
        finish!(progress)
    end
end

if !isempty(ARGS)
    name = ARGS[1]
    path = isabspath(name) ? name : joinpath(@__DIR__, name)
    @info "Save results to $path"
    serialize(path, bd)
end
