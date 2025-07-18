# Mathematical Model

The core of the `NetworkDynamics.jl` package is the [`Network`](@ref) function. It accepts the functions describing the 
local dynamics on the edges and nodes of the graph `g` as inputs, and returns a composite function compatible with the 
DifferentialEquations.jl syntax as output.

```julia
nd = Network(g, vertex_dynamics,  edge_dynamics)
nd(dx, x, p, t)
```

In general the local dynamics on the edges and nodes of a graph can be described through the use of (a) algebraic 
equations, (b) differential algebraic equation (DAEs) in mass matrix form or (c) ordinary differential equations (ODE).
The `NetworkDynamics.jl` package uses 
[Differential-Algebraic-Equation (DAE)](https://mathworld.wolfram.com/Differential-AlgebraicEquation.html) to express 
the overall network dynamics:
```math
M\,\frac{\mathrm{d}}{\mathrm{d}t}u = f^{\mathrm{nw}}(u, p, t)
```
where $M$ is a (possibly singular) mass matrix, $u$ is the internal state vector of the system, $p$ are the parameters
and $t$ is the time. To make this compatible with the solvers used in `OrdinaryDiffEq.jl`, the generated
[`Network`](@ref) object is callable
```
nw(du, u, p, t) # mutates du as an "output"
```
and represents the right-hand-side (RHS) of the equation above. The mass-matrix $m$ is stored in the `Network` object 
as well.

## Modelling the Dynamics of the System
Each component model $\mathrm c$ is modeled as a general input-output-system:
```math
\begin{aligned}
M_{\mathrm c}\,\frac{\mathrm{d}}{\mathrm{d}t}x_{\mathrm c} &= f^{\mathrm c}(x^{\mathrm c}, i_{\mathrm c}, p_{\mathrm c}, t)\\
y^{\mathrm c} &= g^{\mathrm c}(x^\mathrm{c}, i_{\mathrm c}, p_{\mathrm c}, t)
\end{aligned}
```
where $M_{\mathrm{c}}$ is the component mass matrix, $x^{\mathrm c}$ are the component states, $i^{\mathrm c}$ are the
inputs of the component and $y^{\mathrm c}$ is the output of the component. If 
$\mathrm{dim}(x^{\mathrm{c}}) = 0$, the number of internal states is 0.

The mathematical model of `NetworkDynamics.jl` splits the network system in two parts: the vertex and
the edge components (the nodes and edges, respectively). Instead of defining the $f^{\mathrm{nw}}$ by hand, `ND.jl`
builds it automatically based on a list of decentralized nodal and edge dynamics that the user provides (the 
`VertexModel` and `EdgeModel` objects).

In the context of the network, the **output of the edges are flow variables** and the **outputs of vertices are 
potential variables**. When the node and edge models are placed on a graph, the inputs and outputs are connected: 
the nodes receive the output of the adjacent edges as inputs and the edges receive the output of the adjacent nodes as 
inputs. Thus, the *flow* on the edges depends on the *potentials* at both ends as inputs. The *potentials* of the nodes 
depend on the incoming *flows* from all connected edges as an input. (Here, flow and potentials are meant in a 
conceptional and not necessarily physical way.)

In this graphical representation of a partial network graph
```@raw html
<img src="../assets/mathmodel.svg" width="100% height="100%"/>
```
three nodes are visible (node 1, node 2 and node 3) as well as the edges connecting node 1 and node 2 ($e_{\mathrm 12}$, 
$e_{\mathrm 21}$). Above the network, the mass matrix equations on node 1 and node 2 ($M_{\mathrm{c}} x_{\mathrm c}$), 
the equations on the connecting edges ($e_{\mathrm 12}$, $e_{\mathrm 21}$), as well as the internal state vector 
equations of node1 and node2 ($u_1$ and $u_2$) are also shown.

(@Hans: the .svg graphics in this page looks black with dark letters for the most part in some browsers (e.g. firefox, 
duckduckgo). I tried to add a white background to them, but it does not seem to have worked)

### Vertex Models
The equations of a (single-layer) full vertex model are: 
```math
\begin{aligned}
M^{\mathrm v}\,\frac{\mathrm{d}}{\mathrm{d}t}x^{\mathrm v} &= f^{\mathrm v}(x^{\mathrm v}, i^{\mathrm v}, p^{\mathrm v}, t)\\
y^{\mathrm v} &= g^{\mathrm v}(x^{\mathrm v}, i^{\mathrm v}, p^{\mathrm v}, t)
\end{aligned}
```
and they correspond to the Julia functions:
```julia
function fᵥ(dxᵥ, xᵥ, e_aggr, pᵥ, t)
    # mutate dxᵥ
    nothing
end
function gᵥ(yᵥ, xᵥ, e_aggr, pᵥ, t)
    # mutate yᵥ
    nothing
end
vertf = VertexModel(; f=fᵥ, g=gᵥ, mass_matrix=Mᵥ, ...)
```

A (single-layer) full vertex model has one input, and one output. Its input is an aggregation/reduction over all the 
*incident edge outputs* which is calculated using:
```math
i^{\mathrm v} = \mathop{\mathrm{agg}}\limits_k^{\text{incident}} y^{\mathrm e}_k \qquad\text{often}\qquad
i^{\mathrm v} = \sum_k^{\text{incident}} y^{\mathrm e}_k
```

The graphical representation of such a model is:
```@raw html
<img src="../assets/nodemodel.svg" width="70% height="70%"/>
```
where $y^e_i$ and $y^e_j$ are two of the $n$ incident edge outputs that are aggregated to produce the model input 
$i^u$ and the model output $y^v$ (the vertex model output).

 
### Edge Models
In contrast to the vertex models, edge models in general have *two* inputs and *two* outputs, for both the source and 
the destination end of the edge. We commonly use `src` and `dst` to describe the source and destination end of an edge,
respectively.

!!! note "Source and Destination definitions                                                                                                      "
    A source node (`src`) is the node from which the flow exits.
    A destination node (`dst`) is the node into which the flow enters.

!!! note "On the directionality of edges"
    Mathematically, in a system defined on an undirected graph there is no difference between edge $(1,2)$ and
    edge $(2,1)$, because the edge has no direction. However, from an implementation point of view we always need to
    have some kind of ordering, which is why we introduce the source (`src`) and destination (`dst`) terminology.

The full edge model equations are:
```math
\begin{aligned}
M^{\mathrm e}\,\frac{\mathrm{d}}{\mathrm{d}t}x^{\mathrm e} &= f^{\mathrm e}(u^{\mathrm e}, y^{\mathrm v}_{\mathrm{src}}, y^{\mathrm v}_{\mathrm{dst}}, p^{\mathrm e}, t)\\
y^{\mathrm e}_{\mathrm{dst}} &= g_\mathrm{dst}^{\mathrm e}(u^{\mathrm e}, y^{\mathrm v}_{\mathrm{src}}, y^{\mathrm v}_{\mathrm{dst}}, p^{\mathrm e}, t)\\
y^{\mathrm e}_{\mathrm{src}} &= g_\mathrm{src}^{\mathrm e}(u^{\mathrm e}, y^{\mathrm v}_{\mathrm{src}}, y^{\mathrm v}_{\mathrm{dst}}, p^{\mathrm e}, t)
\end{aligned}
```
and they correspond to the Julia functions:
```julia
function fₑ(dxₑ, xₑ, v_src, v_dst, pₑ, t)
    # mutate dxᵥ
    nothing
end
function gₑ(y_src, y_dst, xᵥ, v_src, v_dst, pₑ, t)
    # mutate y_src and y_dst
    nothing
end
vertf = EdgeModel(; f=fₑ, g=gₑ, mass_matrix=Mₑ, ...)
```

The *inputs* of an edge connecting two nodes are the *outputs* of the nodes at both their ends. The output of each node 
is split into two parts:
1. the *`dst`* output which is used as the input of the vertex at the destination end
2. the `src` output which is used as the input of the vertex at the `src` end.

Because a Vertex only receives flows from the edges connected to it, it does not know whether the flow it
received was produced by the source or the destination end of an edge. To solve this problem a sign convention has been
introduced. A positive flow represents a flow *into* the connected vertex, while a negative flow represents a flow *out*
of the connected vertex.

When we consider $n_1$ as the source and $n_2$ as the destination, flows out of $n_1$ and into $n_2$ are negative,
so $y_dst < 0$. Whereas flows out of $n_2$ and into $n_1$ are positive, so $y_src > 0$.
```
                 y_dst < 0 
  n_1 (src) o───────→────────o n_2 (dst)
                 y_src > 0 
  n_1 (src) o───────←────────o n_2 (dst)
```  
But when we consider $n_2$ as the source and $n_1$ as the destination, flows out of $n_2$ and into $n_1$ are negative, 
so $y_dst < 0$. Whereas flows out of $n_1$ and into $n_2$ are positive, so $y_src > 0$.
```
                y_src > 0
  $n_1$ (dst) o───────→────────o $n_2$ (src)
                y_dst < 0
  $n_1$ (dst) o───────→────────o $n_2$ (src)
```

For undirected graphs, `Graphs.jl` chooses the direction of an edge (so which node is the `src` and which the `dst`) 
$v_1->v_2$ such that $v_1 < v_2$, so the edge between vertices 16 and 12 will always be an edge with source 
`src=12` and destination `dst=16`.



### Single Sided Edge Outputs
To simplify the calculations, we split the output flows of edges into "single sided" edge output functions. That way we
can split the output flows of edges into "single sided" edge output functions which simplifies the calculations:
```julia
function g_single(y, xᵥ, v_src, v_dst, pₑ, t)
    # mutate y
    nothing
end
```

There are multiple wrappers available to automatically convert them into double-sided edge output functions:
- `Directed(g_single)` builds a double-sided function *which only couples* to the destination side.
- `Symmetric(g_single)` builds a double-sided function in which both ends receive `y`.
- `AntiSymmetric(g_single)` builds a double-sided function where the destination receives `y` and the source receives `-y`.
- `Fiducial(g_single_src, g_singl_dst)` builds a double-sided edge output function based on two single sided functions.


## Feed Forward Behavior
!!! warning "Feed Forward Vertices"
    As of 11/2024, vertices with feed forward behaviour (FF) are not supported at all. Use [`ff_to_constraint`](@ref) to 
    transform them into vertex models without FF.

Component models can have a so-called Feed Forward behaviour, which provides a direct link between the input and the 
output.

The most generic version of the component models can contain direct FFs from the input to the output. This means that 
the output function $g$ depends directly on the component inputs $i$ rather than just on the component state $x$.

Whenever possible, you should define output functions without FFs in the following way:
```julia
gᵥ_noff(yᵥ, xᵥ, pᵥ, t)
gₑ_noff([y_src,] y_dst, xᵥ, pₑ, t)
```
instead of the more general
```julia
gᵥ(yᵥ, xᵥ, e_aggr, pᵥ, t)
gₑ([y_src], y_dst, xᵥ, v_src, v_dst, pₑ, t)
```

NetworkDynamics cannot couple two components with FFs to each other. But, it is always possible to transform 
feed forward behaviour to an internal state `x` with mass matrix entry zero to circumvent this problem. This 
transformation can be performed automatically using [`ff_to_constraint`](@ref).

Concretely, NetworkDynamics distinguishes between 4 types of feed forward behaviours of `g` functions based on the 
[`FeedForwardType`](@ref) trait.
The feed forward type is inferred automatically based on the provided function `g` (by calculating the signature of the 
methods used over the number of arguments it is given.)

The code bloke below presents the different `g` signatures for the different feed forward types:

**[`PureFeedForward()`](@ref)**
```julia
g!(outs...,          ins...,       p, t) # abstractly
g!(out_dst,          v_src, v_dst, p, t) # single-sided edge
g!(out_src, out_dst, v_src, v_dst, p, t) # double-sided edge
g!(v_out,            e_aggr,       p, t) # single layer vertex
```
**[`FeedForward()`](@ref)**
```julia
g!(outs...,          x, ins...,       p, t) # abstractly
g!(out_dst,          x, v_src, v_dst, p, t) # single-sided edge
g!(out_src, out_dst, x, v_src, v_dst, p, t) # double-sided edge
g!(v_out,            x, e_aggr,       p, t) # single layer vertex
```
**[`NoFeedForward()`](@ref)**
```julia
g!(outs...,          x, p, t) # abstractly
g!(out_dst,          x, p, t) # single-sided edge
g!(out_src, out_dst, x, p, t) # double-sided edge
g!(v_out,            x, p, t) # single layer vertex
```
**[`PureStateMap()`](@ref)**
```julia
g!(outs...,          x) # abstractly
g!(out_dst,          x) # single-sided edge
g!(out_src, out_dst, x) # double-sided edge
g!(v_out,            x) # single layer vertex
```

If the automatic inference of feed forward type fails, the user may specify it explicitly using the `ff` keyword 
argument of the Edge/VertexModel constructor.
