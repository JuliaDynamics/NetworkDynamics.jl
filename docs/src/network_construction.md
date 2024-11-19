# Network Construction

## Building a Network
The main type of `NetworkDyanmics.jl` is a [`Network`](@ref).
A network bundles various component models (edge and vertex models) together with a graph to form a callable object which represents the RHS of the overall dynamical system, see [Mathematical Model](@ref).

A `Network` is build by passing a graph `g`, vertex models `vertexm` and edge models `edgem`.
```julia
nw = Network(g, vertexm, edgem; kwargs...)
```

Two important keywords for the [`Network`](@ref) constructor are:

- `execution`: 
    Defines the [`ExecutionStyle`](@ref) of the coreloop, e.g. `SequentialExecution{true}()`.
    A execution style is a special struct which tells the backend how to parallelize for example.
    A list of available executions styles can be found under [Execution Types](@ref) in the API.

- `aggregator`: 
    Tells the backend how to aggregate and which aggregation function to use.
    Aggregation is the process of creating a single vertex input by reducing over
    the outputs of adjecent edges of said vertex. The `aggregator` contains both the
    function and the algorith. E.g. `SequentialAggregator(+)` is a sequential
    aggregation by summation. A list of availabe Aggregators can be found under
    [`Aggregators`](@ref) in the API.

## Building `VertexModel`s
This chapter walks through the most important aspects when defining custom vertex model. For a list of all keyword arguments please check out the docstring of [`VertexModel`](@ref).
As an example, we'll construct an second order kuramoto model, because that's what we do.
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
It is still annoying to explicitly write this trivial output function. You can prevent this by using [`StateMask`](@ref).
By writing
```@example construction
VertexModel(; f=kuramoto_f!, g=StateMask(1:1), dim=2, pdim=3)
```
we told the vertex model, that the output is part of the states `x[1:1]`.
This enables a few things:
- `outdim` is not needed anymore, can be inferred from `StateMask`
- `outsym` is not a generic `:o` any more but inferred from the state symbols.

We can be even less verbose by writing `g=1:1` or just `g=1`.

In a last we define better names for our states and parameters as well as assigning a position in the graph to enable the graphless network construction.
Whenever you provide `sym` keyword the corresponding `dim` keyword is not neccessary anymore. We end up with a relatively short definition
```@example construction
VertexModel(; f=kuramoto_f!, g=1,
              sym=[:θ, :ω], psym=[:M=>1, :P=>0.1, :D=>0], 
              insym=[:P_nw], name=:swing, vidx=1)
```

## Building `EdgeModel`s
This chapter walks through the most important aspects when defining custom edge models. For a list of all keyword arguments please check out the docstring of [`EdgeModel`](@ref).

As an example edge model we want to define standard sinusoidal coupling between the vertices in our network. The full definition looks like this:

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
which no classifies as a `PureFeedForward` edge.
In cases like this, where the edge is actually anti symmetric we can alternatively define a single sided output function and wrapping it in an `AntiSymmetric` object
```@example construction
function edge_g_s!(ydst, vsrc, vdst, p, t)
    ydst[1] = p[1] * sin(vsrc[1] - vdst[1])
end
EdgeModel(;g=AntiSymmetric(edge_g_ff!), pdim=1, outdim=1)
```
which can also lead to briefer output naming. Available single sided wrappers are
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


