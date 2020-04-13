# NetworkDynamics

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://fhell.github.io/NetworkDynamics.jl/dev)


A package for working with dynamical systems on complex networks. NetworkDynamics.jl provides an interface between [LightGraphs.jl](https://github.com/JuliaGraphs/LightGraphs.jl) and [DifferentialEquations.jl](https://github.com/JuliaDiffEq/DifferentialEquations.jl). It allows to define several types of dynamic and static nodes and edges and to link them up in order to create complex network dynamics.

At the moment the behaviour of a node or an edge can be described by algebraic equations or by ordinary differential equations (ODE). Stochastic ordinary differential equations (SDE) can be implemented as a [two-layer network](https://github.com/FHell/NetworkDynamics.jl/blob/master/examples/sde.jl) - this is still an experimental feature, though. Support for delay differential equations (DDE) is planned for a future release.

## Installation

Installation is straightforward with Julia's package manager.

```julia-repl
(v1.3) pkg> add NetworkDynamics
```

## Getting started

Check out our step-by-step tutorial as a [jupyter notebook](https://github.com/FHell/NetworkDynamics.jl/blob/master/examples/getting_started_with_network_dynamics.ipynb) or [in the docs](https://fhell.github.io/NetworkDynamics.jl/dev/getting_started_with_network_dynamics/).

## PowerDynamics

[PowerDynamics.jl](https://juliaenergy.github.io/PowerDynamics.jl/stable/) is an open-source framework for dynamic power grid modeling and analysis build on top of NetworkDynamics.jl.
