using DifferentiationInterface
import ForwardDiff, ReverseDiff, Enzyme, Zygote, FiniteDiff, FiniteDifferences

(isinteractive() && @__MODULE__()==Main ? includet : include)("ComponentLibrary.jl")

f(x) = sum(abs2, x)

x = [1.0, 2.0]

value_and_gradient(f, AutoForwardDiff(), x) # returns (5.0, [2.0, 4.0]) with ForwardDiff.jl
value_and_gradient(f, AutoEnzyme(),      x) # returns (5.0, [2.0, 4.0]) with Enzyme.jl
value_and_gradient(f, AutoZygote(),      x) # returns (5.0, [2.0, 4.0]) with Zygote.jl


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

# jacobian(fx, AutoEnzyme(), x0)
# jacobian(fx, AutoZygote(), x0)
# jacobian(fx, AutoFiniteDifferences(), x0)
jacobian(fx, AutoForwardDiff(), x0)
jacobian(fx, AutoReverseDiff(), x0)
jacobian(fx, AutoFiniteDiff(), x0)


fp = function(p)
    dx = similar(x0)
    nw(dx, x0, p, 0.0)
    dx
end

# jacobian(fp, AutoEnzyme(), pflat(p0)
# jacobian(fp, AutoZygote(), pflat(p0))
# jacobian(fp, AutoFiniteDifferences(), pflat(p0))
jacobian(fp, AutoForwardDiff(), pflat(p0))
jacobian(fp, AutoReverseDiff(), pflat(p0))
jacobian(fp, AutoFiniteDiff(), pflat(p0))
