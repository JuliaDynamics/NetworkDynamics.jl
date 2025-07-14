# Mathematical Model

The core of the `NetworkDynamics.jl` package is the [`Network`](@ref) function. It accepts the functions describing the 
local dynamics on the edges and nodes of the graph `g` as inputs, and returns a composite function compatible with the 
DifferentialEquations.jl syntax.

```julia
nd = Network(g, vertex_dynamics,  edge_dynamics)
nd(dx, x, p, t)
```

The local dynamics on the edges and nodes of the graph can be described through the use of (a) algebraic equations, 
(b) differential algebraic equation (DAEs) in mass matrix form or (c) ordinary differential equations (ODE). 

The `NetworkDynamics.jl` package uses 
[Differential-Algebraic-Equation (DAE)](https://mathworld.wolfram.com/Differential-AlgebraicEquation.html)
to express the overall network dynamics:
```math
M\,\frac{\mathrm{d}}{\mathrm{d}t}u = f^{\mathrm{nw}}(u, p, t)
```
where $M$ is a (possibly singular) mass matrix, $u$ is the internal state vector of the system, $p$ are the parameters
and $t$ is the time. To make this compatible with the solvers used in `OrdinaryDiffEq.jl`, the generated
[`Network`](@ref) object is a callable object:
```
nw(du, u, p, t) # mutates du as an "output"
```
which represents the right-hand-side (RHS) of the equation above. The mass-matrix $m$ is stored in the `Network` object 
as well.

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
builds it automatically based on a list of decentralized nodal and edge dynamics that the user provides ( 
the so-called `VertexModel` and `EdgeModel` objects).

In the context of the network, the **output of the edges are flow variables** and the **outputs of vertices are 
potential variables**. When the node and edge models are placed on a graph, the inputs and outputs are connected: 
the nodes receive the output of the adjacent edges as inputs and the edges receive the output of the adjacent nodes as 
inputs. Thus, the *flow* on the edges depends on the *potentials* at both ends as inputs. The *potentials* of the nodes 
depend on the incoming *flows* from all connected edges as an input. (Here, flow and potentials are meant in a 
conceptional and not necessarily physical way.)

```@raw html
<img src="../assets/mathmodel.svg" width="100% height="100%"/>
```
Figure 2: Here part of a network is shown. Three nodes are visible (node1, node 2 and node3) as well as the edges
connecting node1 and node2 ($e_12$, $e_21$). Above the network, the mass matrix equations on node1 and node2 
($M_{\mathrm{c}}x^{\mathrm c}$), the equations on the connecting edges ($e_12$, $e_21$), as well as the internal 
state vector equations of node1 and node2($u_1$ and $u_2$) are also shown.
(@Hans: this graphic looks black with dark letter for the most part in some browsers. It may be a better idea to add
a white background to it)

## Vertex Models
A (single-layer) vertex model has one input, and one output.
The input is an aggregation/reduction over all of the *incident edge outputs*:
```math
i^{\mathrm v} = \mathop{\mathrm{agg}}\limits_k^{\text{incident}} y^{\mathrm e}_k \qquad\text{often}\qquad
i^{\mathrm v} = \sum_k^{\text{incident}} y^{\mathrm e}_k
```
```@raw html
<img src="../assets/nodemodel.svg" width="100% style="red"/>
``` 

The full vertex model
```math
\begin{aligned}
M^{\mathrm v}\,\frac{\mathrm{d}}{\mathrm{d}t}x^{\mathrm v} &= f^{\mathrm v}(x^{\mathrm v}, i^{\mathrm v}, p^{\mathrm v}, t)\\
y^{\mathrm v} &= g^{\mathrm v}(x^{\mathrm v}, i^{\mathrm v}, p^{\mathrm v}, t)
\end{aligned}
```
corresponds to the Julia functions
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

## Edge Models
In contrast to vertex models, edge models in general have *two* inputs and *two* outputs, for both the source and the
destination end of the edge. We commonly use `src` and `dst` to describe the source and destination end of an edge,
respectively.

The *inputs* of the edge are the outputs of the two nodes at both their ends. The output is split into two parts:
the `dst` output goes to the input of the vertex at the destination end, the `src` output goes to the input of the
vertex at the `src` end.

For undirected graphs, `Graphs.jl` chooses the direction of an edge `v1->v2` such that `v1 < v2`, i.e. the edge between
vertices 16 and 12 will be always an edge with source `src=12` and destination `dst=16`.

!!! note "On the directionality of edges"
    Mathematically, in a system defined on an undirected graph there is no difference between edge $(1,2)$ and 
    edge $(2,1)$, because the edge has no direction. However, from an implementation point of view we always need to have 
    some kind of ordering, which is why we introduce the source and destination terminology. 

```@raw html
<img src="../assets/edgemodel.svg" width="100%"/>
```

The full model of an edge
```math
\begin{aligned}
M^{\mathrm e}\,\frac{\mathrm{d}}{\mathrm{d}t}x^{\mathrm e} &= f^{\mathrm e}(u^{\mathrm e}, y^{\mathrm v}_{\mathrm{src}}, y^{\mathrm v}_{\mathrm{dst}}, p^{\mathrm e}, t)\\
y^{\mathrm e}_{\mathrm{dst}} &= g_\mathrm{dst}^{\mathrm e}(u^{\mathrm e}, y^{\mathrm v}_{\mathrm{src}}, y^{\mathrm v}_{\mathrm{dst}}, p^{\mathrm e}, t)\\
y^{\mathrm e}_{\mathrm{src}} &= g_\mathrm{src}^{\mathrm e}(u^{\mathrm e}, y^{\mathrm v}_{\mathrm{src}}, y^{\mathrm v}_{\mathrm{dst}}, p^{\mathrm e}, t)
\end{aligned}
```
which corresponds to the Julia functions:
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

The sign convention for both outputs of an edge must be identical, so typically, a positive flow represents a flow 
*into* the connected vertex. This is important, because the vertex only receives the flows, it does not know whether 
the flow was produced by the source or the destination end of an edge.
```
          y_src     y_dst 
  V_src o───←─────────→───o V_dst

```


### Single Sided Edge Outputs
Often, edge outputs will possess some symmetry. This makes it more convenient to define
"single sided" edge output functions:
```julia
function g_single(y, xᵥ, v_src, v_dst, pₑ, t)
    # mutate y
    nothing
end
```
There are multiple wrappers available to automatically convert them into double-sided edge
output functions:

- `Directed(g_single)` builds a double-sided function *which only couples* to the destination side.
- `Symmetric(g_single)` builds a double-sided function in which both ends receive `y`.
- `AntiSymmetric(g_single)` builds a double-sided function where the destination receives `y` and the source receives `-y`.
- `Fiducial(g_single_src, g_singl_dst)` builds a double-sided edge output function based on two single sided functions.

## Feed Forward Behavior
Component models can show have so-called feed forward behavior. Feed forward means, that there is a direct link from input to output.

The most generic version of the component models can contain direct feed forwards from the input to the output.
This means, the output function depends $g$ directly depends on the component inputs $i$ and not only on the component state $x$.

Whenever possible, you should define output functions without feed forwards, i.e.
```julia
gᵥ_noff(yᵥ, xᵥ, pᵥ, t)
gₑ_noff([y_src,] y_dst, xᵥ, pₑ, t)
```
instead of the more general
```julia
gᵥ(yᵥ, xᵥ, e_aggr, pᵥ, t)
gₑ([y_src], y_dst, xᵥ, v_src, v_dst, pₑ, t)
```

NetworkDynamics cannot couple two components with feed forward to each other.
But, it is always possible to transform feed forward behavior to an internal state `x` with mass matrix entry zero to 
circumvent this problem. This transformation can be performed automatically by using [`ff_to_constraint`](@ref).


!!! warning "Feed Forward Vertices"
As of 11/2024, vertices with feed forward are not supported at all. Use [`ff_to_constraint`](@ref) to transform them 
into vertex model without FF.

Concretely, NetworkDynamics distinguishes between 4 types of feed forward behaviors of `g` functions based on the 
[`FeedForwardType`](@ref) trait.
The feed forward type is inferred automatically based on the provided function `g`, more concretely it is determined by the signature of that method / the number of arguments it takes.

The code blow shows the different `g` signatures for the different feed forward types
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

If the automatic inference of feed forward type fails, the user may specify it explicitly using the `ff` keyword argument of the Edge/VertexModel constructor.
