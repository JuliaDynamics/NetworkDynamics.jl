using NetworkDynamics
using SparseConnectivityTracer
using Graphs
using OrdinaryDiffEqRosenbrock

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
_p0 = NWParameter(nw)
_p0.e[2:3,:K] .= 1.0
p0 = pflat(_p0)


get_jac_prototype(nw)
get_jac_prototype(nw; dense=true)
get_jac_prototype(nw; dense=[EIndex(1), VIndex(1)])

# prob = ODEProblem(nw, x0, (0.0, 1.0), p0)
# _nw = ODEFunction(nw; jac_prototype=get_jac_prototype(nw))
# prob_jac = ODEProblem(_nw, x0, (0.0, 1.0), p0)
# @b solve($prob, $Rodas5P())
# @b solve($prob_jac, $Rodas5P())
