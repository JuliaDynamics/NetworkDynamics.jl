# You will find a step-by-step guide to this example in the docs and the
# corresponding jupyter notebook on our github repository.

using DelimitedFiles
using SimpleWeightedGraphs, Graphs
using NetworkDynamics
using OrdinaryDiffEq
using Plots


### Reading the weight matrix (a brain atlas, [N.Tzourio-Mazoyer, 2002])

G = readdlm(joinpath(@__DIR__, "Norm_G_DTI.txt"), ',', Float64, '\n');


### Constructing the graph

g_weighted = SimpleWeightedDiGraph(G)

# for later use we extract the edge.weight attributes with the getfield function
# . is the broadcasting operator and gets the attribute: weight of every edge
edge_weights = getfield.(collect(edges(g_weighted)), :weight);

# we promote the g_weighted graph as a directed graph (weights of the edges are used in parameters)
g_directed = SimpleDiGraph(g_weighted)

### Setting up network_dynamics

```
Fitz-Hugh Nagumo vertex with electrical gap junctions
```
@inline Base.@propagate_inbounds function fhn_electrical_vertex!(dv, v, edges, p, t)
    dv[1] = v[1] - v[1]^3 / 3 - v[2]
    dv[2] = (v[1] - a) * ϵ # x=(u,v)^T
    for e in edges
        dv[1] += e[1]
    end
    nothing
end

@inline Base.@propagate_inbounds function electrical_edge!(e, v_s, v_d, p, t)
    e[1] =  p * (v_s[1] - v_d[1]) # p is the coupling strength
    nothing
end

electricaledge = StaticEdge(f = electrical_edge!, dim = 1, coupling = :directed)
# since the vertex is two dimensional, we specify both symbols u,v
odeelevertex = ODEVertex(f = fhn_electrical_vertex!, dim = 2, sym=[:u, :v])

fhn_network! = network_dynamics(odeelevertex, electricaledge, g_directed)

# global parameters that are accessed several times should be `const` to improve performance

N = 90 # number of nodes in the brain atlas
const ϵ = 0.05 # time separation parameter
const a = .5 # threshold parameter for an oscillatory regime
const σ = .5 # coupling strength

# when passing a tuple of parameters for nodes and edges, then p is passed indi-
# vidually to every edge. For details see the docs on this example

p = (nothing, σ * edge_weights)

x0 = randn(2N) * 5
tspan = (0., 200.)
prob  = ODEProblem(fhn_network!, x0, tspan, p)
sol   = solve(prob, AutoTsit5(TRBDF2()))

### Benchmarking of the simulation

using BenchmarkTools
display(@benchmark sol   = solve(prob, AutoTsit5(TRBDF2())))
display(@benchmark fhn_network!($x0,$x0,$p,0.))

### Plotting

plot(sol, vars = idx_containing(fhn_network!, :u), legend = false, ylim=(-5, 5))
