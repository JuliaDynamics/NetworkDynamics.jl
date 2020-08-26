# Accessing edge values

 In this tutorial we explain how to access the edge values at any time $t$ of the simulation. `NetworkDynamics.jl` offers several possibilities for this purpose. In the following, we will present two different ways how the values of the edges can be retrieved at a time $t$ to be able to access them later.

## Accessing edge values via helper function GetGD

`NetworkDynamics.jl` includes the GetGD type, which provides access to the underlying GraphData and GraphStruct objects by multiple dispatching of the network dynamics functions. Subsequently, we use the example XXX to demonstrate the usage of `GetGD`. In this example a barabasi-albert graph 'g' with $N = 20$ nodes was created. Then, using the functions `diffusionedge!` and `diffusionvertex!` and the functions wrappers, objects of the type `VertexFunction` and `EdgeFunction` were formed. The key constructor `network_dynamics` combines these objects with the graph `g` to the object `nd` (for more information see `getting-started.jl`).

To now get access to the edge values at any time t, the object nd can be called as follows:

```@example accessing_edge_values
gd_nd = nd(x, p, t, GetGD) # exposes underlying graph data struct
e_values_1 = gd_nd.e_array
nothing # hide
```

By calling GetGD in the arguments of object X, the newly created object gd_nd contains information about the underlying graph data struct. For t a time point is chosen at which the information about the graph should be retrieved. The argument x has to be the precomputed state of the system at time t, e.g. from a solution object.

## Accessing edge values using SavingCallback
