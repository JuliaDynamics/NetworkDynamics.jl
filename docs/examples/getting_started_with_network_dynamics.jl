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
(@Hans: this creates a copy of this Julia file. I am assuming you would only want the commands rather than the text)
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

#### Defining an `EdgeModel`
To define an Edgemodel we use the function `diffusionedge_g!`. It takes as as inputs the current state of the edge `e`,
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