# NetworkDynamics

A package for working with dynamical systems on complex networks. NetworkDynamics.jl provides an interface between [LightGraphs.jl](https://github.com/JuliaGraphs/LightGraphs.jl) and [DifferentialEquations.jl](https://github.com/JuliaDiffEq/DifferentialEquations.jl). It allows to define several types of dynamic and static nodes and edges and to link them up in order to create complex network dynamics.

At the moment the behaviour of a node or an edge can be described by algebraic equations or by ordinary differential equations (ODE). Support for stochastic differential equations (SDE) and delay differential equations (DDE) will be added in future releases.

## Installation

Installation is straightforward with Julia's package manager.

```julia-repl
(v1.3) pkg> add NetworkDynamics
```

## PowerDynamics

[PowerDynamics.jl](https://juliaenergy.github.io/PowerDynamics.jl/stable/) is an open-source framework for dynamic power grid modeling and analysis build on top of NetworkDynamics.jl.


# Overview

The key construction is a callable function compatible with the
DifferentialEquations.jl calling syntax.

```julia
nd = network_dynamics(vertices!::Array{VertexFunction}, edges!::Array{EdgeFunction}, g)
nd(dx, x, p, t)
```

The first two parameters are the functions, or arrays of functions that define the dynamics on the network. The types VertexFunction and EdgeFunction are specified on the next page. The last parameter g is a graph encoding the network in the LightGraphs.jl format.

This page is still in development. :)
