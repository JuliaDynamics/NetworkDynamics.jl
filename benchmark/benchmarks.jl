using Pkg

using BenchmarkTools
using Graphs
using NetworkDynamics
using Random
using Serialization
using SciMLBase
(isinteractive() ? includet : include)("benchmark_utils.jl")

if pkgversion(NetworkDynamics) < v"0.9.0"
    @info "Define compatibility functions"
    struct SequentialExecution{P} end
    struct KAExecution{P} end
    function Network(g, v, e; execution = SequentialExecution())
        if execution isa KAExecution
            network_dynamics(v, e, g; parallel=true)
        elseif execution isa SequentialExecution
            network_dynamics(v, e, g)
        else
            error("execution type not supported")
        end
    end
    dim(nd::ODEFunction) = length(nd.syms)
    (isinteractive() ? includet : include)("benchmark_models_v0.8.jl")
else
    (isinteractive() ? includet : include)("benchmark_models.jl")
end

@info "Benchmark with $(Threads.nthreads()) threads"

bd = BenchmarkDict()

####
#### diffusion benchmarks on a dense watts strogatz graph
####
vertex = diffusion_vertex()
edges = Dict("static_edge" => diffusion_edge(),
             "ode_edge" => diffusion_dedge())

for k in ["static_edge", "ode_edge"]
    edge = edges[k]

    for N in [100, 1_000]#, 5_000]  #, 100_000, 1_000_000]
        @info "Benchmark diffusion network with $k on size $N"
        g = watts_strogatz(N, N ÷ 2, 0.0; seed=1)
        b = @b Network($g, $vertex, $edge)
        bd["diffusion", k, "assemble", N] = b

        nd1 = Network(g, vertex, edge; execution=SequentialExecution{true}())
        nd2 = Network(g, vertex, edge; execution=KAExecution{true}())
        x0 = randn(dim(nd1))
        dx = similar(x0)

        bd["diffusion", k, "call", N]    = @b $nd1($dx, $x0, nothing, 0.0)
        bd["diffusion", k, "call_mt", N] = @b $nd2($dx, $x0, nothing, 0.0)
    end
end

####
#### kuramoto benchmarks on a sparse graph
####

function homogeneous(N)
    rng = MersenneTwister(1)
    g = watts_strogatz(N, 3, 0.8; seed=1)
    edge = static_kuramoto_edge()
    vertex = kuramoto_vertex_2d()
    if pkgversion(NetworkDynamics) < v"0.9.0"
        p = (randn(rng, nv(g)), randn(rng, ne(g)))
    else
        p = rand(rng, nv(g) + ne(g))
    end
    (p, vertex, edge, g)
end

function heterogeneous(N)
    rng = MersenneTwister(1)
    g = watts_strogatz(N, 3, 0.8; seed=1)
    edge = static_kuramoto_edge()
    vertex = [kuramoto_vertex_1d(), kuramoto_vertex_2d()]
    vertices = vertex[shuffle(rng, vcat([1 for _ in 1:N÷2], [2 for _ in 1:N÷2]))]
    if pkgversion(NetworkDynamics) < v"0.9.0"
        p = (randn(rng, nv(g)), randn(rng, ne(g)))
    else
        p = rand(rng, nv(g) + ne(g))
    end
    (p, vertices, edge, g)
end

for f in [homogeneous, heterogeneous]
    name = string(f)
    for N in [100, 1_000, 10_000]#, 100_000, 1_000_000]
        @info "Benchmark kuramoto $name on size $N"
        (p, v, e, g) = f(N)
        bd["kuramoto", name, "assemble", N] = @b Network($g, $v, $e)

        nd1 = Network(g, v, e; execution=SequentialExecution{true}())
        nd2 = Network(g, v, e; execution=KAExecution{true}())
        x0 = randn(dim(nd1))
        dx = similar(x0)
        bd["kuramoto", name, "call", N]    = @b $nd1($dx, $x0, $p, 0.0)
        bd["kuramoto", name, "call_mt", N] = @b $nd2($dx, $x0, $p, 0.0)
    end
end

if !isempty(ARGS)
    name = ARGS[1]
    path = isabspath(name) ? name : joinpath(@__DIR__, name)
    @info "Save results to $path"
    serialize(path, bd)
end
