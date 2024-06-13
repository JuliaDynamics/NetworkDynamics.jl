#=
# Neurodynamic model of synchronization in the human brain

#### Topics covered in this tutorial include:

  - constructing a directed, weighted graph from data
  - some useful macros
  - parameter handling
  - stiff equations

## The FitzHugh-Nagumo model

Dynamics of spiking neurons have been described in a simplified manner by the [FitzHugh-Nagumo model](https://en.wikipedia.org/wiki/FitzHugh%E2%80%93Nagumo_model).

```math
\begin{aligned}
\varepsilon \dot u &  =  u - \frac{u^3}{3} - v \\
\dot v & =  u + a
\end{aligned}
```

Here $u$ is a fast, excitatory variable corresponding to the membrane potential and $v$ is a slower, inhibitory varibale. $\varepsilon$ is a parameter separating these time-scales, and $a$ is a control parameter.

In simplified models of the brain, such *relaxation oscillators* may be used to model individual neurons, clusters of neurons or even larger areas in the brain. The FitzHugh-Nagumo model has been widely used for studying synchronization in neuronal activity, which in turn has been connected to physiological phenomena such as epileptic seizures.

## Coupling relaxation oscillators

While different coupling schemes for FitzHugh-Nagumo oscillators have been proposed, in this tutorial we focus on coupling of the excitatory variables via electrical gap junctions, as described by the following system of equations.

```math
\begin{aligned}
\varepsilon \dot u_i & =  u_i - \frac{u_i^3}{3} - v_i - \sigma \sum_{j=1}^N G_{ij}(u_i - u_j) \\
\dot v_i & =   u_i + a
\end{aligned}
```

This is a simple diffusive coupling mediated by the difference between activation potentials in pairs of neurons. A similar coupling term was introduced in the "getting started" tutorial.

## The network topology - a brain atlas

In the following we will use a directed and weighted network encoding the strength and directionality of coupling between 90 different areas of the brain [Nathalie Tzourio-Mazoyer et al., 2002, Neuroimage].

The network weight matrix is given as a text file containing 90 lines with 90 numbers representing the coupling strength and separated by commas `,`. The data can be conveniently read into a matrix with the `DelimitedFiles` module.
=#

using DelimitedFiles
## adjust the load path for your filesystem!
G = readdlm(joinpath(@__DIR__, "Norm_G_DTI.txt"), ',', Float64, '\n')
nothing # hide

#=
The data structure for directed, weighted graphs is provided by the package `SimpleWeightedGraphs.jl` which is based on `Graphs.jl`.
=#

using SimpleWeightedGraphs, Graphs

## First we construct a weighted, directed graph
g_weighted = SimpleWeightedDiGraph(G)

## For later use we extract the edge.weight attributes
## . is the broadcasting operator and gets the attribute :weight for every edge
edge_weights = getfield.(collect(edges(g_weighted)), :weight)

## we promote the g_weighted graph as a directed graph (weights of the edges are included in parameters)
g_directed = SimpleDiGraph(g_weighted)

nothing # hide

#=
## Setting up the ODEProblem

Defining `VertexFunction` and `EdgeFunction` is similar to the example before. The macro `Base.@propagate_inbounds` tells the compiler to inline the function and propagate the inbounds context. For more details see the julia [documentation](https://docs.julialang.org/en/v1/devdocs/boundscheck/).

=#

using NetworkDynamics

Base.@propagate_inbounds function fhn_electrical_vertex!(dv, v, edges, p, t)
    dv[1] = v[1] - v[1]^3 / 3 - v[2]
    dv[2] = (v[1] - a) * ϵ
    for e in edges
        dv[1] += e[1]
    end
    nothing
end

Base.@propagate_inbounds function electrical_edge!(e, v_s, v_d, p, t)
    e[1] = p * (v_s[1] - v_d[1]) # * σ
    nothing
end

odeelevertex = ODEVertex(; f=fhn_electrical_vertex!, dim=2, sym=[:u, :v]);
electricaledge = StaticEdge(; f=electrical_edge!, dim=1, coupling=:directed)

fhn_network! = network_dynamics(odeelevertex, electricaledge, g_directed)

nothing # hide

#=
Since this system is a directed one with thus directed edges, the keyword argument `coupling` is used to set the coupling of the edges to `:directed`.

Note that the multiplication with the coupling strength $\sigma$ has been commented out. Since $\sigma$ is the same for every edge we can absorb this multiplication into the edge weight parameter `p`. Since our network has almost 8000 edges, this saves 8000 multiplications at every function call and leads to an 8-fold increase in performance.

## Parameter handling
=#

# Defining global parameters

N = 90         # number of nodes
const ϵ = 0.05 # global variables that are accessed several times should be declared `const`
const a = 0.5
const σ = 0.5

# Tuple of parameters for nodes and edges

p = (nothing, σ * edge_weights)

# Initial conditions

x0 = randn(2N) * 5

nothing # hide

#=
The behaviour of `network_dynamics` changes with the type of parameters `p` being passed. When `p` is an `Array`, the entire Array will be passed to each `VertexFunction` and `EdgeFunction`. When `p` is a tuple of two Arrays with lengths corresponding to the number of nodes and number of edges respectively, then `network_dynamics` passes only the edge or node parameters with the index of the edge or node. When there are no parameters for either edges or nodes the value `nothing` may be used.

## Solving the system

Now we are ready to create an `ODEProblem`. Since for some choices of parameters the FitzHugh-Nagumo model is *stiff* (i.e. numerically unstable), we use a solver with automated stiffness detection. Such a solver switches to a more stable solver only when the solution enters a region of phase space where the problem is numerically unstable. In this case we use `Tsit5` and switch to `TRBDF2` when necessary. `AutoTsit5` is the switching version of the `Tsit5` algorithm.
=#

using OrdinaryDiffEq

tspan = (0.0, 200.0)
prob  = ODEProblem(fhn_network!, x0, tspan, p)
sol = solve(prob, AutoTsit5(TRBDF2()));
nothing # hide

#=
## Plotting

The plot of the excitatory variables shows that they synchronize for this choice of parameters.
=#

using Plots

plot(sol; vars=idx_containing(fhn_network!, :u), legend=false, ylim=(-5, 5), fmt=:png)
