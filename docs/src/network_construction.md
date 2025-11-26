# Network Construction

## Building a Network
The main object type of `NetworkDynamics.jl` is the [`Network`](@ref). It bundles together the various component models 
(edge and vertex models) with a graph to form a callable object. This object represents the right hand side (RHS) of the 
overall dynamical system (for more information see [Mathematical Model](@ref)).

A `Network` is build by passing a graph `g`, as well as the vertex models `vertexm` and the edge models `edgem` to the 
[`Network`](@ref) network constructor:
```julia
nw = Network(g, vertexm, edgem; kwargs...)
```

In order to perform the simulation we need to tell the backend how to parallelize the system and how to aggregate the
inputs and outputs at the node and edge level. To achieve this, we use the `execution` and `aggregator` keywords of 
the [`Network`](@ref) constructor:

- `execution`:
    Defines the [`ExecutionStyle`](@ref) of the coreloop, e.g. `SequentialExecution{true}()`. 
    (@Hans: the correloop is mentioned here for the first time, please add a definition somewhere and reference it here)
    The [`ExecutionStyle`](@ref) is a special Julia object that tells the backend how to parallelize the system. 
    For example `ThreadedExecution{true}()` [ThreadedExecution](@ref) will use native Julia threads to parallelize the 
    RHS call of the system.
    A list of available executions styles can be found under [Execution Types](@ref) in the API.

- `aggregator`:
    Instructs the backend on how to perform the aggregation and on which aggregation function to use.

!!! note "Regarding Aggregation"
    Aggregation is the process of creating a single vertex input by reducing over the outputs of adjacent edges of 
    said vertex. The `aggregator` contains both the function and the algorithm. E.g. `SequentialAggregator(+)` is a 
    sequential aggregation by summation. A list of available Aggregators can be found under [`Aggregators`](@ref) in the 
    API chapter.

(@Hans: the novice user has no idea how to apply the above. An example is necessary)

### Graphless Constructor
If each of the network components has a "graphelement" `g` [metadata](@ref Metadata) predefined, we can omit the 
explicit graph `g` in the network definition `nw` (for more details on the "graphelement" see  [Metadata](@ref)).

To create a graphless constructor, we first have to define a `Network` object (while omiting the "graphelement" `g`):
```julia
nw = Network(vertexm, edgem)
```
The "graphelement" `g` metadata can then be set using the following syntax options:
```julia
graphName_vertex = VertexModel(; ..., vidx=1)       # places a vertex at position 1
graphName_edge = EdgeModel(; ..., src=1, dst=2)     # places an edge between 1 and 2
graphName_edge = EdgeModel(; ..., src=:v1, dst=:v2) # places an edge between vertices with names `:v1` and `:v2`
```
where "graphName" is the name given to the graph.

## Building `VertexModel`s
This section will walk you through the most important aspects of defining a custom vertex model. 
For a list of all keyword arguments check out [`VertexModel`](@ref).

Next we will demonstrate how to construct a second order Kuramoto model.
```@example construction
using NetworkDynamics #hide
function kuramoto_f!(dv, v, esum, p, t)
    M, P, D = p
    dv[1] = v[2]
    dv[2] = (P - D*v[2] + esum[1])/M
    nothing
end
function kuramoto_g!(y, v, esum, p, t)
    y[1] = v[1]
    nothing
end
kuramoto_vertex = VertexModel(; f=kuramoto_f!, g=kuramoto_g!, dim=2, pdim=3, outdim=1)
```
Those keywords are the minimum amount of metadata we need to provide.

However, there is a problem: the vertex is classified as a `FeedForward` vertex, which is unnecessary. 
We can improve the implementation of `g` according to the [Feed Forward Behavior](@ref) section. 
(@Hans: the Feed Forward has been deprecated according to the [Mathematical Mode](@ref).)
```@example construction
function kuramoto_g_noff!(y, v, p, t)
    y[1] = v[1]
    nothing
end
kuramoto_vertex = VertexModel(; f=kuramoto_f!, g=kuramoto_g_noff!, dim=2, pdim=3, outdim=1)
```

To simplify the programming experience we can use [`StateMask`](@ref) and replace `g=kuramoto_g_noff!` with 
`g=StateMask(1:1)` by writing:
```@example construction
kuramoto_vertex = VertexModel(; f=kuramoto_f!, g=StateMask(1:1), dim=2, pdim=3)
```
Here we are instructing the vertex model, that the output is part of the states `x[1:1]`.
This results in the following changes:
- `outdim=1` is removed because it can be inferred from `StateMask`
- `outsym` is not a generic `:o` anymore but inferred from the state symbols. 

(@Hans: outsym has not been referenced before this line, so the user has no idea what it is. Also, the "generic `:o`" 
has not been explained anywhere before)

We can be even less verbose by replacing `g=StateMask(1:1)` with `g=1:1` or even just `g=1`.
```@example construction
kuramoto_vertex = VertexModel(; f=kuramoto_f!, g=1:1, dim=2, pdim=3)
```
or
```@example construction
kuramoto_vertex = VertexModel(; f=kuramoto_f!, g=1, dim=2, pdim=3)
```

Lastly, we define improved names for our states by using the input parameter `name=:swing` and assign a position in the 
graph to enable the graphless network construction. Whenever a `sym` keyword 
(@Hans: are the sym keywords `psym`, `insym` and `sym` or just the `sym`) is provided, the
corresponding `dim` keyword (@Hans: does the `dim` mentioned here include `pdim` too?) stops being necessary.
So, we end up with the relatively short definition:
```@example construction
VertexModel(; f=kuramoto_f!, g=1,
              sym=[:θ, :ω], psym=[:M=>1, :P=>0.1, :D=>0],
              insym=[:P_nw], name=:swing, vidx=1)
```
(@Hans: this definition is actually longer than the one above so either the above sentence needs to change or the
example needs to change)

## Building `EdgeModel`s
This chapter walks us through the most important aspects of defining a custom Edge Model. For a list of all keyword 
arguments available, check the docstring of [`EdgeModel`](@ref).

To begin with, we define a standard edge model which uses sinusoidal coupling between the vertices in its network.
```@example construction
function edge_f!(de, e, vsrc, vdst, p, t)
    nothing
end
function edge_g!(ysrc, ydst, e, vsrc, vdst, p, t)
    ydst[1] = p[1] * sin(vsrc[1] - vdst[1])
    ysrc[1] = -ydst[1]
end
EdgeModel(; f=edge_f!, g=edge_g!, dim=0, pdim=1, outdim=1)
```
This is a purely "static" edge model, where the nodes have no internal states. This means we can omit the `f` and `dim` 
input parameters entirely. We can also define a variant of `g` without the `e` input for the same reason.
```@example construction
function edge_g_ff!(ysrc, ydst, vsrc, vdst, p, t)
    ydst[1] = p[1] * sin(vsrc[1] - vdst[1])
    ysrc[1] = -ydst[1]
end
EdgeModel(;g=edge_g_ff!, pdim=1, outdim=1)
```
This classifies the edges of the model as `PureFeedForward` edges.
In cases like this, where the edge is actually anti-symmetrical (because ysrc=-ydst) we can define a single sided output
function and wrap it in an `AntiSymmetric` object:
```@example construction
function edge_g_s!(ydst, vsrc, vdst, p, t)
    ydst[1] = p[1] * sin(vsrc[1] - vdst[1])
end
EdgeModel(;g=AntiSymmetric(edge_g_ff!), pdim=1, outdim=1)
```
This can also lead to briefer output naming. Available single sided wrappers are:
- [`Directed`](@ref) (no coupling at `src`),
- [`AntiSymmetric`](@ref) (same coupling at `src` and `dst`),
- [`Symmetric`](@ref) (inverse coupling at `dst`) and
- [`Fiducial`](@ref) (define separate `g` for both ends).

Once again we can add additional data like defining a `src` and `dst` index using the following commands:
```@example construction
function edge_g_s!(ydst, vsrc, vdst, p, t)
    ydst[1] = p[1] * sin(vsrc[1] - vdst[1])
end
EdgeModel(;g=AntiSymmetric(edge_g_ff!), psym=:K=>1, outsym=:P, insym=:θ, src=1, dst=4)
```

So this is how you can construct a second order Kuramoto model using the NetworkDynamics.jl package.
