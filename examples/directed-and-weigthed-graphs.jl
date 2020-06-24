# You will find a step-by-step guide to this example in the docs and the
# corresponding jupyter notebook on our github repository.

using DelimitedFiles
using SimpleWeightedGraphs, LightGraphs
using NetworkDynamics
using OrdinaryDiffEq
using Plots


### Reading the weight matrix (a brain atlas, [N.Tzourio-Mazoyer, 2002])

G = readdlm(joinpath(@__DIR__, "Norm_G_DTI.txt"), ',', Float64, '\n');


### Constructing the graph

g_weighted = SimpleWeightedDiGraph(G)

# For later use we extract the edge.weight attributes
# . is the broadcasting operator and gets the attribute :weight of every edge
edge_weights = getfield.(collect(edges(g_weighted)), :weight);


### Setting up network_dynamics

```
Fitz-Hugh Nagumo vertex with electrical gap junctions
```
@inline Base.@propagate_inbounds function fhn_electrical_vertex!(dv, v, e_s, e_d, p, t)
    dv[1] = v[1] - v[1]^3 / 3 - v[2]
    dv[2] = (v[1] - a) * ϵ
    for e in e_s
        dv[1] -= e[1]
    end
    for e in e_d
        dv[1] += e[1]
    end
    nothing
end
@inline Base.@propagate_inbounds function electrical_edge!(e, v_s, v_d, p, t)
    e[1] =  p * (v_s[1] - v_d[1]) # * σ
    nothing
end

electricaledge = StaticEdge(f! = electrical_edge!, dim = 1)
odeelevertex = ODEVertex(f! = fhn_electrical_vertex!, dim = 2, sym=[:u, :v])

fhn_network! = network_dynamics(odeelevertex, electricaledge, g_weighted)

# Global parameters that are accessed several times should be `const` to improve performance

N = 90 # number of nodes in the brain atlas
const ϵ = 0.05
const a = .5
const σ = .5

# When passing a tuple of parameters for nodes and edges, then p is passed indi-
# vidually to every edge. For details see the docs on this example

p = (nothing, σ * edge_weights)

x0 = randn(2N) * 5
tspan = (0., 200.)
prob  = ODEProblem(fhn_network!, x0, tspan, p)
sol   = solve(prob, AutoTsit5(TRBDF2()))

using BenchmarkTools
display(@benchmark sol   = solve(prob, AutoTsit5(TRBDF2())))
display(@benchmark fhn_network!($x0,$x0,$p,0.))
### Plotting

plot(sol, vars = idx_containing(fhn_network!, :u), legend = false, ylim=(-5, 5))
