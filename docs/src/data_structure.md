# Data Structure

A [`Network`](@ref) contains a list of vertex and edge models along with a graph.
However, in tight numerical loops, it will never access these lists of models directly.
Instead, the network maintains an internal representation that tracks all symbolic indices, defining the precise ordering of states and parameters in a flat array representation. To optimize performance, especially for heterogeneous networks, the network employs specialized data structures that batch identical models together.

This disconnect between the explicit lists and the internal data structures can be confusing.

## Flat Parameter and State Arrays

The vertex and edge models may contain metadata, such as the initial values for states and parameters.
Crucially, this metadata is **only** for the building and initialization of the simulation.
During actual simulation, the state and parameters are handled as **flat arrays**, i.e., plain `Vector{Float64}` objects.

[`NWState`](@ref) and [`NWParameter`](@ref) serve as wrappers around flat arrays and the [`Network`](@ref) objects, allowing you to inspect and modify those flat arrays by addressing vertices and edges directly.

A typical workflow is the following:

1. Set default values in the models using the metadata (see [Metadata](@ref)).
2. Create a network (see [Network Construction](@ref)).
3. Generate a state `s = NWState(nw)` which will be prefilled with the default values from the component metadata (see [Symbolic Indexing](@ref)).
4. Change the values of `s`, i.e., `s.v[1,:x] = 1.0`: This changes the **underlying flat array** but not the metadata of the models.
5. Build a problem with the updated flat arrays using `uflat(s)` and `pflat(s)`.

## Accessing Components
Per default, the models are not copied on Network construction:

```@example data_structure
using NetworkDynamics # hide
using Graphs #hide
include(joinpath(pkgdir(NetworkDynamics), "test", "ComponentLibrary.jl")) # hide
kuramoto_first = Lib.kuramoto_vertex! # hide
kuramoto_secnd = Lib.kuramoto_inertia! # hide
kuramoto_edge = Lib.kuramoto_edge! # hide

v1 = VertexModel(f=kuramoto_first, sym=[:θ], psym=[:ω], g=1)
v2 = VertexModel(f=kuramoto_secnd, sym=[:δ, :ω], psym=[:M, :D, :Pm], g=1)
e = EdgeModel(;g=AntiSymmetric(kuramoto_edge), outsym=[:P], psym=[:K])
nw = Network(complete_graph(2), [v1, v2], e)
```
You can access the models using `getindex`/`[]` with `VIndex` or `EIndex`:
```@example data_structure
v1 === nw[VIndex(1)]
```
This can be important when changing the metadata of components. i.e., both lines below are equivalent:
```@example data_structure
set_position!(v1, (1,0))
set_position!(nw[VIndex(1)], (1,0))
nothing #hide
```

!!! note "Aliasing of component models"
    Since components are not copied, multiple entries in the vertex and edge lists might point to the same instance of a model. 
    ```@example data_structure
    nw = Network(complete_graph(3), [v1,v2,v1], e)
    v1 === nw[VIndex(1)] === nw[VIndex(3)]
    ```
    Consequently, metadata set for one model might affect another model.
    This behavior can be beneficial for performance reasons.
    To force the copying of components, use the `dealias` keyword:
    ```@example data_structure
    nw = Network(complete_graph(3), [v1,v2,v1], e; dealias=true)
    nw[VIndex(1)] === nw[VIndex(3)] # neither of them === v1
    ```

## Extracting a `Network`-object from Containers
`NetworkDynamics.jl` provides a [`extract_nw`](@ref) function, to get a reference to the wrapped `Network` object from different containers, such as solution objects or [integrator objects](https://docs.sciml.ai/DiffEqDocs/stable/basics/integrator). 
