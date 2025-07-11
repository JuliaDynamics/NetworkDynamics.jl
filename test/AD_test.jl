using NetworkDynamics
using Graphs
using DifferentiationInterface
using DifferentiationInterfaceTest
using SparseConnectivityTracer

import ForwardDiff, FiniteDiff, ReverseDiff, Enzyme, Mooncake
import Enzyme: EnzymeCore

(isinteractive() && @__MODULE__()==Main ? includet : include)("ComponentLibrary.jl")

g = complete_graph(4)
vf = [Lib.kuramoto_second(), Lib.diffusion_vertex(), Lib.kuramoto_second(), Lib.diffusion_vertex()]
ef = [Lib.diffusion_odeedge(),
      Lib.kuramoto_edge(),
      Lib.kuramoto_edge(),
      Lib.diffusion_edge_fid(),
      Lib.diffusion_odeedge(),
      Lib.diffusion_edge_fid()]
nw = Network(g, vf, ef)

x0 = rand(dim(nw))
p0 = NWParameter(nw)
p0.e[2:3,:K] .= 1.0

fx = function(x)
    dx = similar(x)
    nw(dx, x, pflat(p0), 0.0)
    dx
end
# jacobian(fx, AutoForwardDiff(), x0)
# jacobian(fx, AutoReverseDiff(), x0)
# jacobian(fx, AutoFiniteDiff(), x0)
# jacobian(fx, AutoMooncake(), x0)
@test_broken jacobian(fx, AutoEnzyme(; mode=EnzymeCore.Forward, function_annotation=EnzymeCore.Duplicated), x0)
@test_broken jacobian(fx, AutoEnzyme(; mode=EnzymeCore.Reverse, function_annotation=EnzymeCore.Duplicated), x0)
@test_broken jacobian(fx, AutoEnzyme(; mode=Enzyme.set_runtime_activity(EnzymeCore.Forward), function_annotation=EnzymeCore.Duplicated), x0)
@test_broken jacobian(fx, AutoEnzyme(; mode=Enzyme.set_runtime_activity(EnzymeCore.Reverse), function_annotation=EnzymeCore.Duplicated), x0)

fp = function(p)
    dx = similar(p,length(x0))
    nw(dx, x0, p, 0.0)
    dx
end
# jacobian(fp, AutoForwardDiff(), pflat(p0))
# jacobian(fp, AutoReverseDiff(), pflat(p0))
# jacobian(fp, AutoFiniteDiff(), pflat(p0))
# jacobian(fp, AutoMooncake(), pflat(p0))
@test_broken jacobian(fp, AutoEnzyme(; mode=EnzymeCore.Forward, function_annotation=EnzymeCore.Duplicated), pflat(p0))
@test_broken jacobian(fp, AutoEnzyme(; mode=EnzymeCore.Reverse, function_annotation=EnzymeCore.Duplicated), pflat(p0))
@test_broken jacobian(fp, AutoEnzyme(; mode=Enzyme.set_runtime_activity(EnzymeCore.Forward), function_annotation=EnzymeCore.Duplicated), pflat(p0))
@test_broken jacobian(fp, AutoEnzyme(; mode=Enzyme.set_runtime_activity(EnzymeCore.Reverse), function_annotation=EnzymeCore.Duplicated), pflat(p0))

scenarios = [Scenario{:jacobian, :in}(fx, x0; res1=jacobian(fx, AutoFiniteDiff(), x0)) ,
             Scenario{:jacobian, :in}(fp, pflat(p0); res1=jacobian(fp, AutoFiniteDiff(), pflat(p0)))]
backends = [AutoForwardDiff(), AutoReverseDiff(), AutoMooncake()]
df = test_differentiation(
    backends,             # the backends you want to compare
    scenarios,            # the scenarios you defined,
    correctness=true,     # compares values against the reference
    type_stability=:none, # checks type stability with JET.jl
    detailed=true,        # prints a detailed test set
    benchmark=:full,
)

# for interactive benchmarking
# df = test_differentiation(
#     backends,
#     scenarios[2:2],
#     correctness=true,
#     type_stability=:none,
#     detailed=true,
#     benchmark=:prepared,
# )

# test sparsity tracer
detector = TracerSparsityDetector();
jacobian_sparsity(fx, x0, detector)
jacobian_sparsity(fp, pflat(p0), detector)
