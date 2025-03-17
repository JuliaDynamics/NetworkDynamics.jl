"""
    BenchmarkCase

Represents a single benchmark case for NetworkDynamics.

Fields:
- name: String identifier for the benchmark
- constructor: Function that creates the Network for a given N
- Ns: Vector of network sizes to benchmark
- description: Documentation string
"""
struct BenchmarkCase
    name::String
    constructor::Function
    Ns::Vector{Int}
    description::String
end

# Standard benchmark cases
BENCHMARK_CASES = [
    BenchmarkCase(
        "diffusion_static_edge",
        N -> begin
            g = watts_strogatz(N, N ÷ 2, 0.0; rng=StableRNG(1))
            Network(g, diffusion_vertex(), diffusion_edge())
        end,
        # [10,20],
        [100, 300, 1000, 3000],
        "Diffusion network with static edges"
    ),

    BenchmarkCase(
        "diffusion_ode_edge",
        N -> begin
            g = watts_strogatz(N, N ÷ 2, 0.0; rng=StableRNG(1))
            Network(g, diffusion_vertex(), diffusion_dedge())
        end,
        # [10,20],
        [100, 300, 1000, 3000],
        "Diffusion network with ODE edges"
    ),

    BenchmarkCase(
        "kuramoto_homogeneous",
        N -> begin
            g = watts_strogatz(N, 3, 0.8; rng=StableRNG(1))
            Network(g, kuramoto_vertex_2d(), static_kuramoto_edge())
        end,
        # [10,20],
        [100, 1_000, 10_000, 100_000],
        "Homogeneous Kuramoto oscillators"
    ),

    BenchmarkCase(
        "kuramoto_heterogeneous",
        N -> begin
            g = watts_strogatz(N, 3, 0.8; rng=StableRNG(1))
            rng = StableRNG(1)
            vtypes = [kuramoto_vertex_1d(), kuramoto_vertex_2d()]
            vertices = vtypes[shuffle(rng, vcat([1 for _ in 1:N÷2], [2 for _ in 1:N÷2]))]
            Network(g, vertices, static_kuramoto_edge())
        end,
        # [10,20],
        [100, 1_000, 10_000, 100_000],
        "Heterogeneous Kuramoto oscillators"
    ),

    BenchmarkCase(
        "powergrid",
        N -> begin
            g = watts_strogatz(N, 3, 0.8; rng=StableRNG(1))
            rng = StableRNG(1)
            vtypes = [pqnode(), generator()]
            vertices = vtypes[shuffle(rng, vcat([1 for _ in 1:N÷2], [2 for _ in 1:N÷2]))]
            Network(g, vertices, piline())
        end,
        # [10,20],
        [100, 1_000, 10_000, 100_000],
        "Power grid network with PQ nodes and generators"
    )
]
