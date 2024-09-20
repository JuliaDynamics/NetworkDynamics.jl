# NetworkDynamics

A package for working with dynamical systems on complex networks. NetworkDynamics.jl provides an interface between [Graphs.jl](https://github.com/JuliaGraphs/Graphs.jl) and [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl). It allows to define several types of dynamic and static nodes and edges and to link them up in order to create complex network dynamics.

The behavior of a node or an edge can be described by algebraic equations, by differential algebraic equation (DAEs) in mass matrix form or by ordinary differential equations (ODE). Stochastic ordinary differential equations (SDE) can e.g. be implemented as a two-layer network. For detail see the [SDE Tutorial](@ref) in the tutorials section.

## Installation

Installation is straightforward with Julia's package manager.

```julia-repl
(v1.10) pkg> add NetworkDynamics
```

## PowerDynamics

[PowerDynamics.jl](https://juliaenergy.github.io/PowerDynamics.jl/stable/) is an open-source framework for dynamic power grid modeling and analysis build on top of NetworkDynamics.jl.


# Overview

The central construction is the function [`Network`](@ref) that receives functions describing the local dynamics on the edges and nodes of
a graph `g` as inputs, and returns a composite function compatible with the
DifferentialEquations.jl calling syntax.

```julia
nd = Network(g, vertices!,  edges!)
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

## Funding
Development of this project was in part funded by the *German Federal Ministry for Economic Affairs and Climate Action* as part of the *OpPoDyn*-Project (Project ID 01258425/1, 2024-2027).

```@raw html
<img src="assets/bmwk_logo_en.svg" width="300"/>
```
