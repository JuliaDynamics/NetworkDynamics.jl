# NetworkDynamics

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://fhell.github.io/NetworkDynamics.jl/stable)


A package for working with dynamical systems on complex networks. NetworkDynamics.jl provides an interface between [LightGraphs.jl](https://github.com/JuliaGraphs/LightGraphs.jl) and [DifferentialEquations.jl](https://github.com/JuliaDiffEq/DifferentialEquations.jl). It allows to define several types of dynamic and static nodes and edges and to link them up in order to create complex network dynamics.

At the moment the behaviour of a node or an edges can be described by algebraic equations or by ordinary differential equations (ODE). Support for stochastic differential equations (SDE) and delay differential equations (DDE) will be added in future releases.

## Installation

Installation is straightforward with Julia's package manager.

```julia-repl
(v1.3) pkg> add NetworkDynamics
```

## PowerDynamics

[PowerDynamics.jl](https://juliaenergy.github.io/PowerDynamics.jl/stable/) is an open-source framework for dynamic power grid modeling and analysis build on top of NetworkDynamics.jl.

## Documentation

The documentation is currently receiving a major rework. The old version can be found [here](https://fhell.github.io/NetworkDynamics.jl/stable).
