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
            (g, diffusion_vertex(), diffusion_edge())
        end,
        # [10,20],
        [100, 300, 1000, 3000],
        "Diffusion network with static edges"
    ),

    BenchmarkCase(
        "diffusion_ode_edge",
        N -> begin
            g = watts_strogatz(N, N ÷ 2, 0.0; rng=StableRNG(1))
            (g, diffusion_vertex(), diffusion_dedge())
        end,
        # [10,20],
        [100, 300, 1000, 3000],
        "Diffusion network with ODE edges"
    ),

    BenchmarkCase(
        "kuramoto_homogeneous",
        N -> begin
            g = watts_strogatz(N, 3, 0.8; rng=StableRNG(1))
            (g, kuramoto_vertex_2d(), static_kuramoto_edge())
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
            (g, vertices, static_kuramoto_edge())
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
    ),

    BenchmarkCase(
        "powergrid_inhomogeneous_pq",
        N -> begin
            g = watts_strogatz(N, 3, 0.8; rng=StableRNG(1))
            rng = StableRNG(1)

            # Create fixed P/Q values for half the nodes
            pq_count = N ÷ 2
            P_vals = -1*rand(rng, pq_count)
            Q_vals = -1*rand(rng, pq_count)

            # Create vertex models with baked-in P/Q values
            pq_vertices = [pqnode_inhomogeneous(P_vals[i], Q_vals[i]) for i in 1:pq_count]
            gen_vertices = [generator() for _ in 1:(N-pq_count)]

            # Shuffle the vertex types
            vertices = shuffle(rng, vcat(pq_vertices, gen_vertices))

            (g, vertices, piline())
        end,
        [100, 1_000, 10_000, 100_000],
        "Power grid with heterogeneous PQ nodes having fixed P/Q values"
    )
]
