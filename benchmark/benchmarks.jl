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
try executions["seq_buf"] = SequentialExecution{true}() catch e end
try executions["ka_buf"] = KAExecution{true}() catch e end
try executions["poly_buf"] = PolyesterExecution{true}() catch e end
try executions["thrd_buf"] = ThreadedExecution{true}() catch e end
try executions["seq"] = SequentialExecution{false}() catch e end
try executions["ka"] = KAExecution{false}() catch e end
try executions["poly"] = PolyesterExecution{false}() catch e end
try executions["thrd"] = ThreadedExecution{false}() catch e end

aggregations = Dict()
try aggregations["ka"] = KAAggregator(+) catch e end
try aggregations["seq"] = SequentialAggregator(+) catch e end
try aggregations["poly"] = PolyesterAggregator(+) catch e end
try aggregations["thrd"] = ThreadedAggregator(+) catch e end
try aggregations["sprs"] = SparseAggregator catch e end

configurations = [
    ("seq", "seq"), # default seq
    ("ka", "seq"),
    ("poly", "seq"),
    ("thrd", "seq"),
    ("seq_buf", "seq"),
    ("ka_buf", "seq"),
    ("poly_buf", "seq"),
    ("thrd_buf", "seq"),
    ("seq_buf", "ka"),
    ("seq_buf", "poly"),
    ("seq_buf", "thrd"),
    ("seq_buf", "sprs"),
    ("thrd_buf", "poly"), # default thrd
]

####
#### diffusion benchmarks on a dense watts strogatz graph
####
vertex = diffusion_vertex()
edges = Dict("static_edge" => diffusion_edge(),
             "ode_edge" => diffusion_dedge())
Ns =  [100, 300, 1000, 3000]#, 10000]  #, 100_000, 1_000_000]

@info "Benchmark diffusion network"
progress = Progress(length(keys(edges)) * length(Ns) * length(configurations),
                    enabled=!haskey(ENV,"GITHUB_ACTIONS"))
for k in keys(edges)
    edge = edges[k]

    for N in Ns
        g = watts_strogatz(N, N ÷ 2, 0.0; rng=StableRNG(1))
        b = @be Network($g, $vertex, $edge) evals=1 samples=10 seconds=2
        bd["diffusion", k, "assemble", N] = BenchmarkResult(b)

        _nd = Network(g, vertex, edge)
        _x0 = randx0(_nd)
        _dx = similar(_x0)
        _nd(_dx, _x0, nothing, NaN)

        for (exname, aggname) in configurations
            GC.gc()
            haskey(executions, exname) || continue
            haskey(aggregations, aggname) || continue
            execution = executions[exname]
            aggregator = aggregations[aggname]
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

            b = @be $nd($dx, $_x0, nothing, 0.0) seconds=2
            br = BenchmarkResult(b, legacy_order(nd, dx))
            bd["diffusion", k, exname, aggname, N] = br
        end
    end
end
finish!(progress)

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

Ns = [100, 1_000, 10_000, 100_000]#, 1_000_000]

@info "Benchmark kuramoto"
progress = Progress(2 * length(Ns) * length(configurations),
                    enabled=!haskey(ENV,"GITHUB_ACTIONS"))
for f in [homogeneous, heterogeneous]
    name = string(f)
    for N in Ns
        (vert, edg, g) = f(N)
        b = @be Network($g, $vert, $edg) evals=1 samples=10 seconds=2
        bd["kuramoto", name, "assemble", N] = BenchmarkResult(b)

        _nd = Network(g, vert, edg)
        _x0 = randx0(_nd)
        _dx = similar(_x0)
        p = randp(_nd)
        _nd(_dx, _x0, p, NaN)

        for (exname, aggname) in configurations
            GC.gc()
            haskey(executions, exname) || continue
            haskey(aggregations, aggname) || continue
            execution = executions[exname]
            aggregator = aggregations[aggname]
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

            b = @be $nd($dx, $_x0, $p, 0.0) seconds=2
            br = BenchmarkResult(b, legacy_order(nd, dx))
            bd["kuramoto", name, exname, aggname, N] = br
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
