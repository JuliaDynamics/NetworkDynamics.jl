using NDPrototype
using Graphs
using Random
using OrdinaryDiffEq
using TimerOutputs

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
    N = 10_000_000
    (p, v, e, g) = heterogeneous(N)
    nd1 = Network(g, v, e; execution=ThreadedExecution{true}());
    x0 = randn(dim(nd1))
    dx1 = zeros(dim(nd1))
    @time nd1(dx1, x0, p, 0.0) # call to init caches, we don't want to benchmark this

    TimerOutputs.enable_debug_timings(NDPrototype)
    nd1(dx1, x0, p, 0.0)
    reset_timer!()
    nd1(dx1, x0, p, 0.0)
    print_timer()
    TimerOutputs.disable_debug_timings(NDPrototype)
    @b $nd1($dx1, $x0, $p, 0.0)

    nd2 = Network(g, v, e; execution=ThreadedExecution{false}());
    dx2 = zeros(dim(nd2))
    @time nd2(dx2, x0, p, 0.0) # call to init caches, we don't want to benchmark this

    TimerOutputs.enable_debug_timings(NDPrototype)
    nd2(dx2, x0, p, 0.0);
    reset_timer!()
    nd2(dx2, x0, p, 0.0)
    print_timer()
    TimerOutputs.disable_debug_timings(NDPrototype)
    @b $nd2($dx2, $x0, $p, 0.0)

    @test dx1 ≈ dx2
end
