using NDPrototype
using Graphs
using Random
using OrdinaryDiffEq
using TimerOutputs
using KernelAbstractions

using Metal
using Adapt

#############
include("testnetwork_constructor.jl")

N = 10_000_000
(p, v, e, g) = heterogeneous(N)

nd1 = Network(g, v, e;
    execution=KAExecution{true}(),
    accumulator=NNlibScatter(+));
nd2 = Network(g, v, e; execution=KAExecution{true}(),
    accumulator=KAAggregator(+));

backend = MetalBackend()

x0 = randn(Float32, dim(nd1));
dx1 = zeros(Float32, dim(nd1));
dx2 = copy(dx1);

x0_d = adapt(backend, x0);
dx1_d = adapt(backend, dx1);
dx2_d = adapt(backend, dx2);
p_d  = adapt(backend, convert.(Float32,p));
nd1_d  = adapt(backend, nd1);
nd2_d  = adapt(backend, nd2);

@b $nd1($dx1, $x0, $p, NaN)
@b $nd1($dx2, $x0, $p, NaN)
@b $nd1_d($dx1_d, $x0_d, $p_d, NaN)
@b $nd2_d($dx2_d, $x0_d, $p_d, NaN)
@test dx1 ≈ dx2 ≈ Vector(dx1_d) ≈ Vector(dx2_d)
