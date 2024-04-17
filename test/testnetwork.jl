using NDPrototype
using Graphs
using Random
using OrdinaryDiffEq

#############
include("../benchmark/benchmarks.jl")
include("../benchmark/benchmark_utils.jl")

# diffusion test
begin
    N = 100000
    g = watts_strogatz(N, 3, 0.8, seed=1)
    edge = diffusion_dedge()
    vertex = diffusion_vertex()
    nd1 = Network(g, vertex, edge, execution=:seq);
    nd2 = Network(g, vertex, edge, execution=:threaded);
    x0 = randn(dim(nd1));
    dx1 = zero(x0);
    dx2 = zero(x0);
    nd1(dx1, x0, nothing, 0.0)
    nd2(dx2, x0, nothing, 0.0)
    @test dx1 ≈ dx2

    @btime $nd1($dx1, $x0, nothing, 0.0)
    @btime $nd2($dx2, $x0, nothing, 0.0)

    prob1 = ODEProblem(nd1, x0, (0.0, 1.0), nothing)
    prob2 = ODEProblem(nd2, x0, (0.0, 1.0), nothing)
    @time sol1a = solve(prob1, Tsit5());
    @time sol1b = solve(prob1, Rodas4());
    @time sol1c = solve(prob1, Rosenbrock23());
    @time sol2a = solve(prob2, Tsit5());
    @time sol2b = solve(prob2, Rodas4());
end

begin
    N = 100_000
    (p, v, e, g) = heterogeneous(N)
    nd1 = Network(g, v, e, execution=:seq)
    (p, v, e, g) = heterogeneous(N)
    nd2 = Network(g, v, e, execution=:threaded)
    x0 = randn(dim(nd1))
    dx1 = zeros(dim(nd1))
    dx2 = zeros(dim(nd1))
    @time nd1(dx1, x0, p, 0.0) # call to init caches, we don't want to benchmark this
    @time nd1(dx2, x0, p, 0.0) # call to init caches, we don't want to benchmark this
    @test dx1 ≈ dx2

    @time nd2(dx1, x0, p, 0.0) # call to init caches, we don't want to benchmark this
    @time nd2(dx2, x0, p, 0.0) # call to init caches, we don't want to benchmark this
    @test dx1 ≈ dx2

    @btime $nd1($dx1, $x0, $p, 0.0)
  # 509.375 μs (6 allocations: 320 bytes)
    @btime $nd2($dx2, $x0, $p, 0.0)
  # 211.208 μs (222 allocations: 17.84 KiB)

    N = 1000
    g = watts_strogatz(N, Int(N/2), 0.0, seed=1)
    gc = NDPrototype.ColoredGraph(g)
    for e in Graphs.edges(gc)
        c = NDPrototype.pickfree(gc, e.src, e.dst)
        NDPrototype.setcolor!(gc, e, c)
        println(e)
    end
    gc
end
