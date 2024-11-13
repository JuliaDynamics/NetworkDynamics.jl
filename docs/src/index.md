# NetworkDynamics

A package for working with dynamical systems on complex networks. NetworkDynamics.jl provides an interface between [Graphs.jl](https://github.com/JuliaGraphs/Graphs.jl) and [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl).
It allows for high performant simulation of dynamic networks by describing the local dynamics on the edges and vertices of the graph.

The behavior of a node or an edge can be described by algebraic equations, by differential algebraic equation (DAEs) in mass matrix form or by ordinary differential equations (ODE). 

The central construction is the function [`Network`](@ref) that receives functions describing the local dynamics on the edges and nodes of
a graph `g` as inputs, and returns a composite function compatible with the
DifferentialEquations.jl calling syntax.

```julia
nd = Network(g, vertex_dynamics,  edge_dynamics)
nd(dx, x, p, t)
```

Main features:
- Clear separation of local dynamics and topology: you can easily change topology of your system or switching out dynamical components.
- High performance when working with heterogeneous models (which means heaving different local dynamics in different parts of your network).
- [Symbolic Indexing](@ref) into solutions and states: NetworkDynamics keeps track of the states of the individual subsystems.
- Different execution schemes: NetworkDynamics exploits the known inter-dependencies between components to auto parallelize execution, even on GPUs!
- Equation based models: use [ModelingToolkit.jl](https://docs.sciml.ai/ModelingToolkit/dev/) to define local dynamics, use `NetworkDynamics.jl` to combine them into large networks!

## Where to begin?
Check out the [Mathematical Model](@ref) to understand the underlying modelling ideas of NetworkDynamics followed by the page on [Network Construction](@ref) to learn how to implement you own models.

If you prefer to look at some concrete code first check out the [Getting Started](@ref) tutorial!


## Installation

Installation is straightforward with Julia's package manager.

```julia-repl
(v1.11) pkg> add NetworkDynamics
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
