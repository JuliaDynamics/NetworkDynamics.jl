# NetworkDynamics

The *NetworkDynamics.jl* package simulates the dynamics of complex networks. It provides an interface 
between the [Graphs.jl](https://github.com/JuliaGraphs/Graphs.jl) and the 
[DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl) packages and facilitates the simulation of 
highly efficient dynamic networks by describing the local dynamics of the edges and vertices of the network.

The core idea of this package is to define the **global dynamics** of a complex network in terms of **local dynamics**: each node and each edge exhibits some local dynamics defined as an input-output system.
The graph topology describes, how the local dynamical systems are interconnected. To learn more check out the docs on the [mathematical model](@ref) behind NetworkDynamics.jl.
For basic terminology see the Wikipedia article on [Graph Theory](https://en.wikipedia.org/wiki/Graph_theory).

Main features:
- Clear separation of local dynamics and topology: you can easily change the topology of your system or switch out dynamic components.
- High performance when working with heterogeneous models: you can have different local dynamics in different parts of your network.
- [Symbolic Indexing](@ref) into solutions and states: NetworkDynamics keeps track of the states of each individual subsystem.
- Diverse execution schemes: NetworkDynamics exploits the known interdependencies between components to auto parallelize execution, even on GPUs!
- Equation based models: you can model local dynamics using 
  [ModelingToolkit.jl](https://docs.sciml.ai/ModelingToolkit/dev/) and then combine them into larger networks using `NetworkDynamics.jl`!


## Where to begin?
To learn how to implement your own models and understand the underlying modeling ideas of NetworkDynamics you should 
first read the [Mathematical Model](@ref) section, followed by the [Network Construction](@ref) section.

If you prefer to look at some concrete code first check out the [Getting Started](@ref) tutorial!


## Installation

1. Install Julia:
-   [Julia Installation](https://julialang.org/install/)
-   Find your OS and follow the instructions for the installation

2. Install NetworkDynamics.jl with Julia's package manager:
```julia-repl
(v1.11) pkg> add NetworkDynamics
```

To learn more about how to use Julia you can visit: [Modern Julia Workflows](https://modernjuliaworkflows.org/)


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
Development of this project was in part funded by the *German Federal Ministry for Economic Affairs and Climate Action* 
as part of the *OpPoDyn*-Project (Project ID 01258425/1, 2024-2027).

```@raw html
<img src="assets/bmwk_logo_en.svg" width="300"/>
```
