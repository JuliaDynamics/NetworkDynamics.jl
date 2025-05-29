#=
# [Network Diffusion](@id Getting Started)

This introductory example explains the use of the basic types and constructors
in NetworkDynamics.jl by modeling a simple diffusion on an undirected network.

This page will walk you through:
  * the theoretical background of a simple diffusion propagating in an undirected network
  * explain the network dynamics of the system
  * explain how to program the dynamics

!!! note
    An Undirected Network is a network where the connections between the nodes are all bidirectional.

This example can be dowloaded as a normal Julia script [here](@__NAME__.jl). #md
(@Hans: this creates a copy of this Julia file. I am assuming you would only want the commands rather than the text.
I have created a file with just the code and minimal documentation and placed it in the subfolder /examples_code. The
name of the file is simple_diffusion_on_an_undirected_network.jl)
=#

#=
## Theoretical background

Diffusion processes appear in phenomena as diverse as heat conduction, electrical currents, and random walks.
Generally speaking they describe the tendency of systems to evolve towards a state of equally distributed entropy
(that entropy being in the form of e.g. heat, charge or concentration).
If we assume a thermal system, the temperature of a specific spot changes depending on the temperature gradient between
itself and its neighborhood.

We will build a graph $g$ with $N$ nodes and an adjacency matrix $A$. $v = (v_1, \dots, v_n)$ is the vector of (abstract)
temperatures at each node $i = 1, \dots, N$. The rate of change of state $v_i$ will be described by the difference
between the temperature of the node and that of its neighbors. For the above we obtain the following ordinary
differential equation:
```math
\dot v_i = \sum_{j=1}^N A_{ji} (v_j - v_i).
```

The sum on the right hand side plays the role of a (discrete) gradient. If the temperature at node $i$ is higher than
at its neighboring node $j$ it will decrease along that edge.

If a node were to be disconnected from the rest of the network its state would never change, because its adjacency
matrix would be $A_{ji} = 0 \; \forall j$ and hence $\dot v_i = 0$. So because its state would remain unchanged the
model where the nodes have no internal dynamics. This means that the evolution of a node will depend only on its
interaction with its neighbors. In NetworkDynamics.jl, interactions between a node and its neighbors are described
using edge equations.

In order to bring this equation into the form required by NetworkDynamics.jl we need split the dynamics into edge and
vertex parts and bring them into the correct input-output formulation.
=#

#=
#### Vertex Dynamics:
All vertices have one internal state $v$. This state is also the vertex output. It is the sum over all incoming flows
of edges connected to the vertex. This directly corresponds to the component model definition outlined in
[Mathematical Model](@ref):
```math
\begin{aligned}
\dot x^\mathrm{v} &= f^{\mathrm v}(u^{\mathrm v}, \sum_k^{\text{incident}} y^{\mathrm e}_k, p^{\mathrm v}, t) &&= \sum_k^\mathrm{incident} y^\mathrm{e}_k \\
y^\mathrm{v} &= g^{\mathrm v}(u^{\mathrm v}, \sum_k^{\text{incident}} y^{\mathrm e}_k, p^{\mathrm v}, t) &&= x^\mathrm{v}
\end{aligned}
```
=#

#=
#### Edge Dynamics:
The edge dynamics on the other hand do not have any internal states. Thus we can define the edge output as the
difference between the source and destination vertex:
```math
\begin{aligned}
y^{\mathrm e}_{\mathrm{dst}} &= g_\mathrm{dst}^{\mathrm e}(u^{\mathrm e}, y^{\mathrm v}_{\mathrm{src}}, y^{\mathrm v}_{\mathrm{dst}}, p^{\mathrm e}, t) &&= y^{\mathrm v}_{\mathrm{src}} - y^{\mathrm v}_{\mathrm{dst}}\\
y^{\mathrm e}_{\mathrm{src}} &= g_\mathrm{src}^{\mathrm e}(u^{\mathrm e}, y^{\mathrm v}_{\mathrm{src}}, y^{\mathrm v}_{\mathrm{dst}}, p^{\mathrm e}, t) &&= y^{\mathrm v}_{\mathrm{dst}} - y^{\mathrm v}_{\mathrm{src}}
\end{aligned}
```
=#

#=
### Modelling dynamics in NetworkDynamics.jl
To model the vertex dynamics we need to create an `VertexModel` and to model the edge dynamics we need to create an
`EdgeModel`.

First we need to have Julia and the necessary packages installed
Then we need to load the necessary packages
=#

using Graphs
using NetworkDynamics
using OrdinaryDiffEqTsit5
using StableRNGs
using Plots
nothing #hide

#=
#### Defining an `EdgeModel`

Then we define the `EdgeModel`.
To it we use the function `diffusionedge_g!`. It takes as as inputs the current state of the edge `e`,
its source vertex `v_src`, its destination vertex `v_dst`, a vector of parameters `p` and the time `t`.
In order to comply with the syntax of NetworkDynamics.jl we must always define functions for edges with exactly
these arguments. (In the case of this example, the values for `p` and `t` are not used). (@Hans Why are they not used?)
After the function call the edge's output value `e` equals the difference between its source and its destination vertex
(i.e. the discrete gradient along that edge).

!!! note
    `diffusionedge_g!` is called a **mutating** function, because it mutates (modifies) the edge state `e` (which is the
    first of its inputs). (I Julia names of mutating functions end with an `!`.)
    We use mutating functions because they reduce allocations and as a result speed up computations because .
=#

function diffusionedge_g!(e_dst, v_src, v_dst, p, t)
    ## e_dst, v_src, v_dst are arrays, hence we use the broadcasting operator
    e_dst .= v_src .- v_dst
    nothing
end
nothing #hide #md

#=
Notably, this function only models $g_\mathrm{dst}$. However we can wrap this single-sided output function in an
[`AntiSymmetric`](@ref) output wrapper to construct the [`EdgeModel`](@ref):
=#

nd_diffusion_edge = EdgeModel(; g=AntiSymmetric(diffusionedge_g!), outsym=[:flow])


#=
#### Definition of `VertexModel`
Next we need to define the `VertexModel`
For undirected graphs, the `edgefunction!` specifies the coupling from a source- to a destination vertex.
The contributions of the connected edges to a single vertex are "aggregated". Default aggregation is the summation of
all incident edge states. The aggregated edge state is made available via the `esum` argument of the vertex function.
=#

function diffusionvertex_f!(dv, v, esum, p, t)
    ## dv, v and esum are arrays, hence we use the broadcasting operator .
    dv .= esum
    nothing
end
nothing #hide #md

#=
Just like above the input arguments `v, esum, p, t` are mandatory for the syntax of vertex functions. The additional
    input `dv` corresponding to the derivative of the vertex' state is mandatory for vertices described by ordinary
    differential equations.

The output function `g` is just taking part of the internal states. For that we can use the [`StateMask`](@ref) helper
    function `g = StateMaks(1:1)`
=#

nd_diffusion_vertex = VertexModel(; f=diffusionvertex_f!, g=StateMask(1:1), dim=1)

#=
## Constructing the network

With the components defined, we can define the topology and assemble the network dynamics.
=#

N = 20 # number of nodes
k = 4  # average degree
g = barabasi_albert(N, k) # a little more exciting than a bare random graph

nothing #hide #md

#= The [Barabási–Albert model](https://en.wikipedia.org/wiki/Barab%C3%A1si%E2%80%93Albert_model) generates a scale-free
random graph.
=#

nd = Network(g, nd_diffusion_vertex, nd_diffusion_edge)

#=
The constructor `Network` combines the component model with the topological information contained in the graph **`g`**
and returns an `Network` compatible with the solvers of `DifferentialEquations.jl`.
=#

rng = StableRNG(1)
x0 = randn(rng, N) # random initial conditions
ode_prob = ODEProblem(nd, x0, (0.0, 2.0))
Main.test_execution_styles(ode_prob) # testing all ex styles #src
sol = solve(ode_prob, Tsit5());
nothing #hide #md

#=
We are solving the diffusion problem on the time interval $[0, 2]$ with the `Tsit5()` algorithm, which is recommended
by the authors of `DifferentialEquations.jl` for most non-stiff problems.
=#

plot(sol; idxs=vidxs(nd, :, :), fmt=:png)

#=
The plotting is straightforward. The **`idxs`** keyword allows us to pass a list of indices. Indices can be also
"symbolic" indices which specify components and their symbols directly. For example `idxs = VIndex(1, :v)` acesses
state `:v` of vertex 1. See [Symbolic Indexing](@ref) for more details.

In oder to collect multiple indices we can use the helper function [`vidxs`](@ref) and [`eidxs`](@ref), which help to
    collect all symbolic indices matching a certain criteria.

## Two Dimensional Extension

To illustrate a very simple multi-dimensional case, in the following we simulate two independent diffusions on an
identical graph. The first uses the symbol `x` and is started with initial conditions drawn from the standard normal
distribution $N(0,1)$, the second uses the symbol `ϕ` with squared standard normal inital conditions.

The symbols have to be passed with the keyword **`sym`** to `VertexModel`.
=#

N = 10 # number of nodes
k = 4  # average degree
g = barabasi_albert(N, k) # a little more exciting than a bare random graph


# We will have two independent diffusions on the network, hence dim = 2

nd_diffusion_vertex_2 = VertexModel(; f=diffusionvertex_f!, g=1:2, dim=2, sym=[:x, :ϕ])
nd_diffusion_edge_2 = EdgeModel(; g=AntiSymmetric(diffusionedge_g!), outsym=[:flow_x, :flow_ϕ])
nd_2 = Network(g, nd_diffusion_vertex_2, nd_diffusion_edge_2)

x0_2 = vec(transpose([randn(rng, N) .^ 2 randn(rng, N)])) # x ~ N(0,1)^2; ϕ ~ N(0,1)
ode_prob_2 = ODEProblem(nd_2, x0_2, (0.0, 3.0))
Main.test_execution_styles(ode_prob_2) # testing all ex styles #src
sol_2 = solve(ode_prob_2, Tsit5());


# Try plotting the variables ϕ_i yourself. [To write ϕ type \phi and press TAB]
plot(sol_2; idxs=vidxs(nd_2, :, :x), fmt=:png)

# Using the `eidxs` helper function we can also plot the flow variables
plot(sol_2; idxs=eidxs(nd_2, :, :flow_x), fmt=:png)

#=
## Appendix: The network Laplacian $L$

The diffusion equation on a network can be rewritten as

```math
\dot v_i  = \sum_{j=1}^N A_{ji} v_j - d_i v_i =  e_i^T A v - d_i v_i
```

where $d_i$ is the degree of node $i$ and $e_i^T$ is the $i$-th standard basis vector. Introducing the diagonal matrix
$D$ that has the degree of node $i$ in its $i$-th row and the Laplacian matrix $L = D - A$ we arrive at

```math
\dot v = e_i^T(A - D) v
```

and finally

```math
\dot v = - L v
```

This is a linear system of ODEs and its solution is a matrix exponential. To study the asymptotic behaviour of the
system it suffices to analyze the eigenspectrum of $L$. For this reason $L$ is an important construction in network
science.
=#

#=
# Putting it all together
We begin by loading the necessary packages.
=#

using Graphs
using NetworkDynamics
using OrdinaryDiffEqTsit5
using StableRNGs
using Plots
nothing #hide

using Graphs
using NetworkDynamics
using OrdinaryDiffEqTsit5
using StableRNGs
using Plots
nothing #hide

function diffusionedge_g!(e_dst, v_src, v_dst, p, t)
    ## e_dst, v_src, v_dst are arrays, hence we use the broadcasting operator
    e_dst .= v_src .- v_dst
    nothing
end
nothing #hide #md

nd_diffusion_edge = EdgeModel(; g=AntiSymmetric(diffusionedge_g!), outsym=[:flow])

function diffusionvertex_f!(dv, v, esum, p, t)
    ## dv, v and esum are arrays, hence we use the broadcasting operator .
    dv .= esum
    nothing
end
nothing #hide #md

nd_diffusion_vertex = VertexModel(; f=diffusionvertex_f!, g=StateMask(1:1), dim=1)

N = 20 # number of nodes
k = 4  # average degree
g = barabasi_albert(N, k) # a little more exciting than a bare random graph
nothing #hide #md

nd = Network(g, nd_diffusion_vertex, nd_diffusion_edge)
rng = StableRNG(1)
x0 = randn(rng, N) # random initial conditions
ode_prob = ODEProblem(nd, x0, (0.0, 2.0))
Main.test_execution_styles(ode_prob) # testing all ex styles #src
sol = solve(ode_prob, Tsit5());
nothing #hide #md
plot(sol; idxs=vidxs(nd, :, :), fmt=:png)
N = 10 # number of nodes
k = 4  # average degree
g = barabasi_albert(N, k) # a little more exciting than a bare random graph
nd_diffusion_vertex_2 = VertexModel(; f=diffusionvertex_f!, g=1:2, dim=2, sym=[:x, :ϕ])
nd_diffusion_edge_2 = EdgeModel(; g=AntiSymmetric(diffusionedge_g!), outsym=[:flow_x, :flow_ϕ])
nd_2 = Network(g, nd_diffusion_vertex_2, nd_diffusion_edge_2)
x0_2 = vec(transpose([randn(rng, N) .^ 2 randn(rng, N)])) # x ~ N(0,1)^2; ϕ ~ N(0,1)
ode_prob_2 = ODEProblem(nd_2, x0_2, (0.0, 3.0))
Main.test_execution_styles(ode_prob_2) # testing all ex styles #src
sol_2 = solve(ode_prob_2, Tsit5());
plot(sol_2; idxs=vidxs(nd_2, :, :x), fmt=:png)
plot(sol_2; idxs=eidxs(nd_2, :, :flow_x), fmt=:png)