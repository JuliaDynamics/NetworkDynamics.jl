# NetworkDynamics

# Overview

This package implements functions for defining and studying dynamics on networks.
The key construction is a callable function compatible with the
DifferentialEquations.jl calling syntax.

```julia
nd = network_dynamics(vertices!::Array{VertexFunction}, edges!::Array{EdgeFunction}, g)
nd(dx, x, p, t)
```

The first two parameters are the functions, or function arrays from which a network dynamics is
built. The types VertexFunction and EdgeFunction are specified on the next page.
The last parameter g is a graph encoding the network constructed with
the LightGraphs.jl package.

This page is still in development. :)
