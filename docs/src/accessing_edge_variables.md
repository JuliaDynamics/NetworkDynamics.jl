# Accessing edge values

 In this tutorial we explain how to access the edge values at any time $t$ of the simulation. `NetworkDynamics.jl` offers several possibilities for this purpose. In the following, we will present two different ways how the values of the edges can be retrieved at a time $t$ to be able to access them later.

## Accessing edge values via helper function GetGD

`NetworkDynamics.jl` includes the GetGD type, which provides access to the underlying GraphData and GraphStruct objects by multiple dispatching of the network dynamics functions. Subsequently, we demonstrate the usage of `GetGD`. For the use of network dynamics we combine functions that describe the dynamics of vertices and edges with the information of the graph created before. For instance, in the example `getting-started.jl` the object `nd` has been formed (for more information see `getting-started.jl`).

To now get access to the edge values at any time t, the object nd can be called as follows:

```@example accessing_edge_values
gd_nd = nd(x, p, t, GetGD) # exposes underlying graph data struct
e_values_1 = gd_nd.e_array
nothing # hide
```

By calling GetGD in the arguments of object X, the newly created object gd_nd contains information about the underlying graph data struct. For t a time point is chosen at which the information about the graph should be retrieved. The argument x has to be the precomputed state of the system at time t, e.g. from a solution object.

## Accessing edge values using SavingCallback
