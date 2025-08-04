# Network Construction

## Building a Network
The main type of `NetworkDynamics.jl` is a [`Network`](@ref).
A network bundles various component models (edge and vertex models) together with a graph to form a callable object which represents the right hand side (RHS) of the overall dynamical system, see [Mathematical Model](@ref).

A `Network` is build by passing a graph `g`, vertex models `vertexm` and edge models `edgem` to the [`Network`](@ref) constructor:.
```julia
nw = Network(g, vertexm, edgem; kwargs...)
```

Two important keywords for the [`Network`](@ref) constructor are:

- `execution`:
    Defines the [`ExecutionStyle`](@ref) of the coreloop, e.g. `SequentialExecution{true}()`.
    A execution style is a special Julia object, which tells the backend how to parallelize (e.g. `ThreadedExecution{true}()` will use native Julia threads to parallelize the RHS call).
    A list of available executions styles can be found under [Execution Types](@ref) in the API.

- `aggregator`:
    Instructs the backend how to perform the aggregation and which aggregation function to use.
    Aggregation is the process of creating a single vertex input by reducing over the outputs of adjecent edges of said vertex. The `aggregator` contains both the function and the algorithm. E.g. `SequentialAggregator(+)` is a sequential aggregation by summation. A list of availabe Aggregators can be found under [`Aggregators`](@ref) in the API.

### Graphless Constructor
If each of the network components has a "graphelement" [metadata](@ref Metadata), we may omit the explicit graph.
```julia
nw = Network(vertexm, edgem)
```
The graphelement metadata can be set using the following syntax:
```julia
VertexModel(; ..., vidx=1)         # places vertex at position 1
EdgeModel(; ..., src=1, dst=2)     # places edge between 1 and 2
EdgeModel(; ..., src=:v1, dst=:v2) # places edge between vertices with names `:v1` and `:v2`
```

## Building `VertexModel`s
This chapter will walk you through the most important aspects of defining a custom vertex model. For a list of all keyword arguments please check out the docstring of [`VertexModel`](@ref).

As an example, we'll construct an second order kuramoto model.
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
VertexModel(; f=kuramoto_f!, g=kuramoto_g!, dim=2, pdim=3, outdim=1)
```
Those keywords are the minimum metadata we need to provide.

However there is a problem: the vertex is classified as a `FeedForward` vertex, which is unnecessary. We can improve the implementation of `g` according to the [Feed Forward Behavior](@ref) section.
```@example construction
function kuramoto_g_noff!(y, v, p, t)
    y[1] = v[1]
    nothing
end
VertexModel(; f=kuramoto_f!, g=kuramoto_g_noff!, dim=2, pdim=3, outdim=1)
```

To simplify your programming and avoid explicitly writing the above trivial output function you can use [`StateMask`](@ref).
By writing
```@example construction
VertexModel(; f=kuramoto_f!, g=StateMask(1:1), dim=2, pdim=3)
```
we are instructing the vertex model, that the output is part of the states `x[1:1]`.
This results in the following changes:
- `outdim` is removed because it can be inferred from `StateMask`
- `outsym` is not a generic `:o` any more but inferred from the state symbols.

We can be even less verbose by writing `g=1:1` or just `g=1`.

Lastly, we define improved names for our states and parameters as well as assigning a position in the graph to enable the graphless network construction.
Whenever you provide a `sym` keyword the corresponding `dim` keyword stops being neccessary. So, we end up with a relatively short definition
```@example construction
VertexModel(; f=kuramoto_f!, g=1,
              sym=[:θ, :ω], psym=[:M=>1, :P=>0.1, :D=>0],
              insym=[:P_nw], name=:swing, vidx=1)
```

## Building `EdgeModel`s
This chapter walks you through the most important aspects when defining custom edge models. For a list of all keyword arguments please check the docstring of [`EdgeModel`](@ref).

As an example edge model we define a standard sinusoidal coupling between the vertices in our network. The full definition is:

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
This is a purely "static" edge without internal states. This means we can omit `f` and `dim` entirely.
Also, we can define a variant of `g` without the `e` input
```@example construction
function edge_g_ff!(ysrc, ydst, vsrc, vdst, p, t)
    ydst[1] = p[1] * sin(vsrc[1] - vdst[1])
    ysrc[1] = -ydst[1]
end
EdgeModel(;g=edge_g_ff!, pdim=1, outdim=1)
```
which classifies as a `PureFeedForward` edge.
In cases like this, where the edge is actually anti-symmetrical we can define a single sided output function and wrap it in an `AntiSymmetric` object:
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

Once again we can add additonal data like defining a `src` and `dst` index
```@example construction
function edge_g_s!(ydst, vsrc, vdst, p, t)
    ydst[1] = p[1] * sin(vsrc[1] - vdst[1])
end
EdgeModel(;g=AntiSymmetric(edge_g_ff!), psym=:K=>1, outsym=:P, insym=:θ, src=1, dst=4)
```
