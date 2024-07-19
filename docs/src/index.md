# NetworkDynamics

A package for working with dynamical systems on complex networks. NetworkDynamics.jl provides an interface between [Graphs.jl](https://github.com/JuliaGraphs/Graphs.jl) and [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl). It allows to define several types of dynamic and static nodes and edges and to link them up in order to create complex network dynamics.

The behavior of a node or an edge can be described by algebraic equations, by differential algebraic equation (DAEs) in mass matrix form or by ordinary differential equations (ODE). Stochastic ordinary differential equations (SDE) can e.g. be implemented as a [two-layer network](https://github.com/pik-icone/NetworkDynamics.jl/blob/master/examples/sde.jl). For detail see the tutorials section.

## Installation

Installation is straightforward with Julia's package manager.

```julia-repl
(v1.10) pkg> add NetworkDynamics
```

## PowerDynamics

[PowerDynamics.jl](https://juliaenergy.github.io/PowerDynamics.jl/stable/) is an open-source framework for dynamic power grid modeling and analysis build on top of NetworkDynamics.jl.


# Overview

The central construction is the function [`network_dynamics`](@ref) that receives functions describing the local dynamics on the edges and nodes of
a graph `g` as inputs, and returns a composite function compatible with the
DifferentialEquations.jl calling syntax.

```julia
nd = network_dynamics(vertices!,  edges!, g)
nd(dx, x, p, t)
```

## Reproducibility

```@raw html
<details><summary>Direct dependencies used for this documentation:</summary>
```

```@example
using Pkg # hide
Pkg.status() # hide
```

```@raw html
</details>
```

```@raw html
<details><summary>Julia Version:</summary>
```

```@example
using InteractiveUtils # hide
versioninfo() # hide
```

```@raw html
</details>
```

```@raw html
<details><summary>Full Manifest:</summary>
```

```@example
using Pkg # hide
Pkg.status(; mode = PKGMODE_MANIFEST) # hide
```

```@raw html
</details>
```
