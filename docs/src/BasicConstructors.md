# Functions

The key constructor [`network_dynamics`](@ref) assembles the dynamics of the whole network from functions for the single vertices and edges of the graph `g`.
Since the equations describing the local dynamics may differ strongly from each
other, the types `VertexFunction` and `EdgeFunction` are introduced. They
provide a unifying interface between different classes of nodes and edges. Both
have several subtypes that account for the different classes of equations that may
represent the local dynamics. At the moment, algebraic (static) equations, ordinary differential equations (ODEs) and delay differential equations (DDEs) are supported:


```julia
# VertexFunctions
StaticVertex(vertexfunction!, dimension, symbol)
ODEVertex(vertexfunction!, dimension, mass_matrix, symbol)
DDEVertex(vertexfunction!, dimension, mass_matrix, symbol)

# EdgeFunctions
StaticEdge(edgefunction!, dimension, coupling, symbol)
ODEEdge(edgefunction!, dimension, mass_matrix, coupling, symbol)
StaticDelayEdge(edgefunction!, dimension, mass_matrix, coupling, symbol)
```


# VertexFunctions

Given a set of (algebraic or differential) equations describing a node or an edge
the first step is to turn them into a **mutating** function `vertexfunction!`. Depending on the class of the function `vertexfunction!`, the constructors `StaticVertex`, `ODEVertex` or `DDEVertex` are called in order to turn `vertexfunction!` into a `VertexFunction` object compatible with [`network_dynamics`](@ref).


Since in general the state of a vertex depends on the vertex value itself as well as on the in- and outgoing edges, the function `vertexfunction!`
has to respect one of the following calling syntaxes.

```julia
# For static nodes
function vertexfunction!(v, edges, p, t) end
# For dynamic nodes
function vertexfunction!(dv, v, edges, p, t) end
# For delay nodes
function vertexfunction!(dv, v, edges, h_v, p, t) end
```

Here `dv`, `v`, `p`, and `t` are the usual ODE arguments, while `edges` is an array containing the edges for which the described vertex is the destination. In the delay case, the array `h_v` denotes the vertex history. The typical case of diffusive coupling on a directed graph could be described as

```julia
function vertex!(dv, v, edges, p, t)
    dv .= 0.
    for e in edges
        dv .+= e
    end
    nothing
end
```

### [`StaticVertex`](@ref)

If a vertex is described by an algebraic equation  `vertexfunction!(v, edges, p, t)`, i.e. `dv = 0` the `VertexFunction` is constructed as

```julia
StaticVertex(vertexfunction!, dim, sym)
```

Here, **dim** is the number of independent variables in the vertex equations and **sym** is an array of symbols for these variables. For example, if a node
models a constant input ``I = p``, then `dim = 1` and `sym = [:I]`. For more details on the use of symbols, check out the [Getting started](@ref getting_started) tutorial and the Julia [documentation](https://docs.julialang.org/en/v1/manual/metaprogramming/). The use of symbols makes it easier to later fish out the interesting variables one wants to look at.


### [`ODEVertex`](@ref)

If a vertex has local dynamics `vertexfunction!(dv, v, edges, p, t)` described by an ODE,
the `VertexFunction` is contructed as

```julia
ODEVertex(vertexfunction!, dim, mass_matrix, sym)
```

As above, **dim** is the number of independent variables in the vertex equations and **sym** corresponds to the symbols of these variables.

**mass_matrix** is an optional argument that defaults to the identity matrix `I`. If a mass matrix M is given, the local system `M * dv = vertexfunction!` will be solved. `network_dynamics` assembles all local mass matrices into one global mass matrix that can be passed to a differential equation solver like `Rodas4`.


One may also call ODEVertex with keyword arguments, omitting optional arguments:

```julia
ODEVertex(f! = vertexfunction!, dim = dim)
```

The function then defaults to using the identity as mass matrix and `[:v for i in 1:dimension]` as symbols.

### [`DDEVertex`](@ref)

If a vertex has local dynamics described by a delay differential equation (DDE) the local dynamics need to have the signature `vertexfunction!(dv, v, edges, h_v, p, t)`, where `h_v` are the history values of `v`. Then the `VertexFunction` is constructed as

```julia
DDEVertex(vertexfunction!, dim, mass_matrix, sym)
```
As before, **dim** is the number of independent variables in the vertex function and **sym** describes the symbols of the variables. The **mass matrix** is optional and defaults to the identity matrix `I`.

As above, the DDEVertex function can also be called without the optional arguments:

```julia
DDEVertex(f! = vertexfunction!, dim = dim)
```

## EdgeFunctions

Similar to the case of vertices, an edge is described by **mutating** function `edgefunction!`. At the moment the constructors `StaticEdge`, `ODEEdge`, and `StaticDelayEdge` are available. `edgefunction!` has to respect one of the following syntaxes:

```julia
# For static edges
function edgefunction!(e, v_s, v_d, p, t) end
# For dynamics edges
function edgefunction!(de, e, v_s, v_d, p, t) end
# For static delay edges
function edgefunction!(e, v_s, v_d, h_v_s, h_v_d, p, t) end
```
Just like above, `de`, `e`, `p` and `t` are the usual ODE arguments, while `v_s`
and `v_d` are the source and destination vertices respectively and `h_v_s` and `h_v_d` the corresponding history arrays.


### [`StaticEdge`](@ref)

Static here means, that the edge value described by `edgefunction!` only depends on the values of the vertices the edge connects to and that no derivative of the edge's internal state is involved. One very simple and natural example is a diffusive edge:

```julia
edgefunction! = (e, v_s, v_d, p, t) -> e .= v_s .- v_d
```

In this case the `EdgeFunction` is constructed by

```julia
StaticEdge(edgefunction!, dim, coupling, sym)
```
The keywords are the same as for the vertices, except from the additional keyword **coupling**, which describes if the `EdgeFunction` is intended for a (`:directed`) or for an undirected graph (`{:undirected, :symmetric, :antisymmetric, :fiducial}`). `:directed` is intended for directed graphs. `:undirected` is the default option and is only compatible with SimpleGraph. In this case the `edgefunction!` should specify the coupling from a source vertex to a destination vertex. `:symmetric` and `:antisymmetric` trigger performance optimizations, if `edgefunction!` has that symmetry property. `:fiducial` lets the user specify both the coupling from source to destination, as well as the coupling from destination to source and is intended for advanced users.

### [`ODEEdge`](@ref)

For problems where `edgefunction!` describes the differential of an edge value, we use the `ODEEdge` function. An example for such a system is given by:

```julia
edgefunction! = (de, e, v_s, v_d, p, t) -> de .= 1000 * (v_s .- v_d .- e)
```
The `EdgeFunction` object is constructed as

```julia
ODEEdge(edgefunction!, dim, coupling, mass_matrix, sym)
```

The keywords are the same as above for the `StaticEdge`. For `ODEEdge` the same simplified construction rules apply when keyword arguments are used.

```julia
ODEEdge(f! = edgefunction!, dim = n)
```

In this case the function defaults to using undirected as coupling, the identity as mass matrix and `[:e for in 1:dimension]` as symbols.

### [`StaticDelayEdge`](@ref)
This constructor is used when edge variables depend on past values of the vertex variables. In this case the `edgefunction!` has to accept two additional arguments `h_v_s` and `h_v_d` that hold the history of `v_s` and `v_d`. *Static* means that the edge depends only on the dynamics of the vertices the edge is connected to and not on an internal derivative of the edge variables itself.

As an example for such system, we show a diffusive coupling with delay:

```julia
edgefunction! = (e, v_s, v_d, h_v_s, h_v_d, p, t) -> e .=.1 * (h_v_s .- v_d)
```
The `EdgeFunction` object is constructed as

```julia
StaticDelayEdge(edgefunction!, dim, coupling, mass_matrix, sym)
```
Again, we can also leave out the optional keywords **coupling**, **sym** and **mass_matrix**.

```julia
StaticDelayEdge(f! = edgefunction!, dim = n)
```

## Constructor

The key constructor is the function [`network_dynamics`](@ref) that takes in
two arrays of `EdgeFunctions` and `VertexFunctions` describing the local dynamics on the edges and nodes of a graph `g`, given as a [LightGraphs.jl](https://github.com/JuliaGraphs/LightGraphs.jl) object. It returns a composite function compatible with the [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl) calling syntax.

```julia
nd = network_dynamics(vertices!::Array{VertexFunction},
                      edges!::Array{EdgeFunction}, g)
nd(dx, x, p, t)
```

If all vertices, respectively edges share the same `VertexFunction` or `EdgeFunction`, than `network_dynamics` can be conveniently called with these functions as arguments.

```julia
nd = network_dynamics(vertexfunction!::VertexFunction,
                      edgefunction!::EdgeFunction, g)

```

The optional keyword argument `parallel` is `false` by default and can be set to `true` if a multi-threaded `ODEFunction` should be used. This may significantly improve performance on multi-core machines, for more details see section [Multi-Threading](@ref).

```julia
nd = network_dynamics(vertexfunction!::VertexFunction,
                      edgefunction!::EdgeFunction, g, parallel=true)

```

### Example

 Let's look at an example. First, we define our graph as well as the differential systems connected to its vertices and edges:

```@example
using NetworkDynamics, LightGraphs

g = erdos_renyi(10, 25) # random graph with 10 vertices and 25 edges

function vertexfunction!(dv, v, edges, p, t)
  dv .= 0
  for e in edges
    dv .+= e
  end
end

function edgefunction!(de, e, v_s, v_d, p, t)
     de .= 1000 .*(v_s .- v_d .- e)
     nothing
end

vertex = ODEVertex(f! = vertexfunction!, dim = 1)
vertexarr = [vertex for v in vertices(g)]

edge = ODEEdge(f! = edgefunction!, dim = 1)
edgearr = [edge for e in edges(g)]

nd = network_dynamics(vertexarr, edgearr, g)

nothing # hide
```

Now we have an `ODEFunction nd` that can be solved with the tools provided by
[DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl).

For more details check out the Tutorials section.
