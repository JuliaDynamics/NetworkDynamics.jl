using NDPrototype
using Graphs
using Random
using OrdinaryDiffEq
using TimerOutputs

#############
include("testnetwork_constructor.jl")

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
    N = 1_000_000
    (p, v, e, g) = heterogeneous(N)

    ####
    nd1 = Network(g, v, e; execution=ThreadedExecution{true}());
    x0 = randn(dim(nd1));
    dx1 = zeros(dim(nd1));
    @time nd1(dx1, x0, p, 0.0) # call to init caches, we don't want to benchmark this

    TimerOutputs.enable_debug_timings(NDPrototype)
    nd1(dx1, x0, p, 0.0)
    reset_timer!()
    nd1(dx1, x0, p, 0.0)
    print_timer()
    TimerOutputs.disable_debug_timings(NDPrototype)
    @b $nd1($dx1, $x0, $p, 0.0)
    # @descend nd1(dx1, x0, p, 0.0)


    ####
    nd2 = Network(g, v, e; execution=ThreadedExecution{false}());
    dx2 = zeros(dim(nd2));
    @time nd2(dx2, x0, p, 0.0) # call to init caches, we don't want to benchmark this
    @test dx1 ≈ dx2

    TimerOutputs.enable_debug_timings(NDPrototype)
    nd2(dx2, x0, p, 0.0);
    reset_timer!()
    nd2(dx2, x0, p, 0.0)
    print_timer()
    TimerOutputs.disable_debug_timings(NDPrototype)
    @b $nd2($dx2, $x0, $p, 0.0)


    ####
    nd3 = Network(g, v, e; execution=SequentialExecution{false}());
    dx3 = zeros(dim(nd3));
    @time nd3(dx3, x0, p, 0.0) # call to init caches, we don't want to benchmark this
    @test dx3 ≈ dx2 ≈ dx1

    TimerOutputs.enable_debug_timings(NDPrototype)
    nd3(dx3, x0, p, 0.0);
    reset_timer!()
    nd3(dx3, x0, p, 0.0)
    print_timer()
    TimerOutputs.disable_debug_timings(NDPrototype)
    @b $nd3($dx3, $x0, $p, 0.0)

    ####
    nd4 = Network(g, v, e; execution=SequentialExecution{true}());
    dx4 = zeros(dim(nd4));
    @time nd4(dx4, x0, p, 0.0) # call to init caches, we don't want to benchmark this
    @test dx4 ≈ dx2 ≈ dx1

    TimerOutputs.enable_debug_timings(NDPrototype)
    nd4(dx4, x0, p, 0.0);
    reset_timer!()
    nd4(dx4, x0, p, 0.0)
    print_timer()
    TimerOutputs.disable_debug_timings(NDPrototype)
    @b $nd4($dx4, $x0, $p, 0.0)
    # @descend nd4(dx4,x0,p,0.0)


    ####
    nd5= Network(g, v, e;
        execution=ThreadedExecution{true}(),
        accumulator=KAAggregator(+));
    dx5 = zeros(dim(nd5));
    @time nd5(dx5, x0, p, 0.0) # call to init caches, we don't want to benchmark this
    @test dx5 ≈ dx5 ≈ dx2 ≈ dx1

    TimerOutputs.enable_debug_timings(NDPrototype)
    nd5(dx5, x0, p, 0.0)
    reset_timer!()
    nd5(dx5, x0, p, 0.0)
    print_timer()
    TimerOutputs.disable_debug_timings(NDPrototype)
    @b $nd5($dx5, $x0, $p, 0.0)

    prob1 = ODEProblem(nd1, x0, (0.0, 1.0), p)
    @time sol1a = solve(prob1, Tsit5());
    @time sol1b = solve(prob1, Rodas5P());
    @time sol1c = solve(prob1, Rosenbrock23());
end
