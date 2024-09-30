using NetworkDynamics
using Graphs
using DifferentiationInterface
using DifferentiationInterfaceTest
import ForwardDiff, ReverseDiff, Enzyme, Zygote, FiniteDiff, FiniteDifferences

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


fp = function(p)
    dx = similar(p,length(x0))
    nw(dx, x0, p, 0.0)
    dx
end
# jacobian(fp, AutoEnzyme(), pflat(p0)
# jacobian(fp, AutoZygote(), pflat(p0))
# jacobian(fp, AutoFiniteDifferences(), pflat(p0))
# jacobian(fp, AutoForwardDiff(), pflat(p0))
# jacobian(fp, AutoReverseDiff(), pflat(p0))
# jacobian(fp, AutoFiniteDiff(), pflat(p0))

scenarios = [Scenario{:jacobian, :in}(fx, x0; res1=jacobian(fx, AutoFiniteDiff(), x0)) ,
             Scenario{:jacobian, :in}(fp, pflat(p0); res1=jacobian(fp, AutoFiniteDiff(), pflat(p0)))]
backends = [AutoForwardDiff(), AutoReverseDiff()]
test_differentiation(
    backends,             # the backends you want to compare
    scenarios,            # the scenarios you defined,
    correctness=true,     # compares values against the reference
    type_stability=false, # checks type stability with JET.jl
    detailed=true,        # prints a detailed test set
)
