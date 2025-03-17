using Pkg

using Chairmarks
using Graphs
using NetworkDynamics
using NetworkDynamics: iscudacompatible
using Serialization
using StableRNGs
using SciMLBase
using Test
using ProgressMeter
using Random
using CUDA
using CUDA: adapt

(isinteractive() ? includet : include)("benchmark_utils.jl")
(isinteractive() ? includet : include)("benchmark_compat.jl")
(isinteractive() ? includet : include)("benchmark_models.jl")
(isinteractive() ? includet : include)("benchmark_cases.jl")

# using ThreadPinning
# pinthreads(:cores)

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
    # differet ex over same agg
    ("seq", "seq"),
    ("ka", "seq"),
    # ("ka", "sprs"), # additional for GPU
    # ("ka", "ka"), # additional for GPU
    ("poly", "seq"),
    ("thrd", "seq"),
    ("seq_buf", "seq"), # default seq
    ("ka_buf", "seq"),
    # ("ka_buf", "sprs"), # additional for GPU
    # ("ka_buf", "ka"), # additional for GPU
    ("poly_buf", "seq"),
    ("thrd_buf", "seq"),
    # different agg over same ex
    ("poly_buf", "ka"),
    ("poly_buf", "poly"),
    ("poly_buf", "thrd"),
    ("poly_buf", "sprs"),
    ("poly_buf", "poly"), # default thrd
]

SECONDS = 1

# Main benchmark loop
for case in BENCHMARK_CASES
    @info "Benchmark $(case.name)"
    progress = Progress(length(case.Ns) * length(configurations),
                       enabled=!haskey(ENV,"GITHUB_ACTIONS"))

    for N in case.Ns
        # Create network components
        (g, vertices, edges) = case.constructor(N)
        
        # Benchmark network construction
        b = @be Network($g, $vertices, $edges) evals=1 samples=10 seconds=SECONDS  
        bd[case.name, "assemble", N] = BenchmarkResult(b)

        # Create network instance
        _nd = Network(g, vertices, edges)

        # Setup initial conditions
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
            
            next!(progress; showvalues = [(:case, case.name), (:N, N), (:ex, exname), (:agg, aggname)])
            
            # Create network with specific execution/aggregation
            nd = try
                Network(_nd; execution, aggregator)
            catch e
                if e.msg == "execution type not supported"
                    continue
                else
                    rethrow(e)
                end
            end
            
            # Test correctness
            dx = similar(_x0)
            nd(dx, _x0, p, NaN)
            @test dx ≈ _dx

            # CPU Benchmark
            b = @be $nd($dx, $_x0, $p, 0.0) seconds=SECONDS
            br = BenchmarkResult(b, dx)
            bd[case.name, exname, aggname, N] = br

            # GPU benchmarks if available
            if CUDA.functional() && iscudacompatible(execution) && iscudacompatible(nd.layer.aggregator)
                # Float32 GPU benchmark
                update!(progress; showvalues = [(:case, case.name), (:N, N), 
                                              (:ex, exname*"_gpu"), (:agg, aggname*"_gpu")])
                to = CuArray{Float32}
                nd_d = adapt(to, nd)
                x0_d = adapt(to, _x0)
                dx_d = similar(x0_d)
                p_d = adapt(to, p)

                nd_d(dx_d, x0_d, p_d, NaN)
                if !(maximum(abs.(Vector(dx_d) - _dx)) < 1e-4)
                    @warn "CUDA f32 lead to different results: extrema(Δ) = $(extrema(Vector(dx_d) - _dx))"
                end

                b = @be $nd_d($dx_d, $x0_d, $p_d, 0.0) seconds=SECONDS
                br = BenchmarkResult(b, _dx)  # Store original _dx due to f32 precision
                bd[case.name, exname*"_cuda32", aggname*"_cuda32", N] = br

                # Float64 GPU benchmark
                to = CuArray{Float64}
                nd_d = adapt(to, nd)
                x0_d = adapt(to, _x0)
                dx_d = similar(x0_d)
                p_d = adapt(to, p)

                nd_d(dx_d, x0_d, p_d, NaN)
                @test Vector(dx_d) ≈ _dx

                b = @be $nd_d($dx_d, $x0_d, $p_d, 0.0) seconds=SECONDS
                br = BenchmarkResult(b, Vector(dx_d))
                bd[case.name, exname*"_cuda64", aggname*"_cuda64", N] = br
            end
        end
    end
    finish!(progress)
end

# Save results if filename provided
if !isempty(ARGS)
    name = ARGS[1]
    path = isabspath(name) ? name : joinpath(@__DIR__, name)
    @info "Save results to $path"
    serialize(path, bd)
end
