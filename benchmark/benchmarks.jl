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

using ThreadPinning
pinthreads(:cores)

if pkgversion(NetworkDynamics) < v"0.9.0"
    (isinteractive() ? includet : include)("benchmark_models_v0.8.jl")
else
    (isinteractive() ? includet : include)("benchmark_models.jl")
end

@info "Benchmark with $(Threads.nthreads()) threads"

bd = BenchmarkDict()

executions = Dict()
# try executions["seq_buf"] = SequentialExecution{true}() catch end
# try executions["ka_buf"] = KAExecution{true}() catch end
# try executions["poly_buf"] = PolyesterExecution{true}() catch end
# try executions["threaded_buf"] = ThreadedExecution{true}() catch end
try executions["seq"] = SequentialExecution{false}() catch end
# try executions["ka"] = KAExecution{false}() catch end
# try executions["poly"] = PolyesterExecution{false}() catch end
try executions["threaded"] = ThreadedExecution{false}() catch end

aggregations = Dict()
# try aggregations["nnlib"] = NNlibScatter(+) catch end
# try aggregations["KA"] = KAAggregator(+) catch end
try aggregations["seq"] = SequentialAggregator(+) catch end
# try aggregations["poly"] = PolyesterAggregator(+) catch end
try aggregations["thrd"] = ThreadedAggregator(+) catch end

####
#### diffusion benchmarks on a dense watts strogatz graph
####
vertex = diffusion_vertex()
edges = Dict("static_edge" => diffusion_edge(),
             "ode_edge" => diffusion_dedge())
Ns =  [300, 1000, 3000]  #, 100_000, 1_000_000]

@info "Benchmark diffusion network"
progress = Progress(length(keys(edges)) * length(Ns) * length(executions) * length(aggregations),
                    enabled=!haskey(ENV,"GITHUB_ACTIONS"))
for k in keys(edges)
    edge = edges[k]

    for N in Ns
        g = watts_strogatz(N, N ÷ 2, 0.0; rng=StableRNG(1))
        b = @be Network($g, $vertex, $edge) evals=1 samples=50 seconds=1
        bd["diffusion", k, "assemble", N] = BenchmarkResult(b)

        _nd = Network(g, vertex, edge)
        _x0 = randx0(_nd)
        _dx = similar(_x0)
        _nd(_dx, _x0, nothing, NaN)

        for (exname, execution) in executions
            for (aggname, aggregator) in aggregations
                next!(progress; showvalues = [(:edge, k), (:N, N), (:ex, exname), (:agg, aggname)])
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
                GC.gc()
            end
        end
    end
end
finish!(progress)

####
#### kuramoto benchmarks on a sparse graph
####
Ns = [100, 1_000]#, 10_000, 100_000]#, 1_000_000]

@info "Benchmark kuramoto"
progress = Progress(2 * length(Ns) * length(executions) * length(aggregations),
                    enabled=!haskey(ENV,"GITHUB_ACTIONS"))
for f in [homogeneous, heterogeneous]
    name = string(f)
    for N in Ns
        (vert, edg, g) = f(N)
        b = @be Network($g, $vert, $edg) evals=1 samples=50 seconds=1
        bd["kuramoto", name, "assemble", N] = BenchmarkResult(b)

        _nd = Network(g, vert, edg)
        _x0 = randx0(_nd)
        _dx = similar(_x0)
        p = randp(_nd)
        _nd(_dx, _x0, p, NaN)

        for (exname, execution) in executions
            for (aggname, aggregator) in aggregations
                next!(progress; showvalues = [(:type, name), (:N, N), (:ex, exname), (:agg, aggname)])
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
                GC.gc()
            end
        end
    end
end
finish!(progress)


if !isempty(ARGS)
    name = ARGS[1]
    path = isabspath(name) ? name : joinpath(@__DIR__, name)
    @info "Save results to $path"
    serialize(path, bd)
end
