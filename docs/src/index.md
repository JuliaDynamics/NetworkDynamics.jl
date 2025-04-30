# NetworkDynamics

*NetworkDynamics.jl* is a package ~~for working with~~ *to simulate* dynamical systems ~~on~~ *within* complex networks. ~~NetworkDynamics.jl~~ *It* provides an interface between *the* [Graphs.jl](https://github.com/JuliaGraphs/Graphs.jl) and *the* [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl) *packages* ~~.~~ *and* ~~it allows for high performant  simulation of~~ *faciliates the simulation of highly efficient* dynamic networks by describing the local dynamics on the edges and vertices of the graph.

*Complex network systems are composed by the entities that comprise them (the nodes) and the relationships that connect each entity with one another (the edges). The graphical depictions of such networks are called graphs. The simplest network (which can be seen in Figure 1) is composed of two entities (so two nodes) who are only connected to each other. This connection between the two is the edge of the system. Complex networks are composed of multiple nodes and edges, with most nodes connected to multiple other nodes with multiple edges*** *(@Hans: can you created the graph of such a network and place in here?)*

The behavior of a node or an edge can be described ~~by~~ *through the use of* a) algebraic equations, b) ~~by~~ differential algebraic equation (DAEs) in mass matrix form ~~or~~ c) ~~by~~ ordinary differential equations (ODE). 

The ~~central construction~~ *core of the package* is the function [`Network`](@ref)*.* ~~that~~ *It accepts the* ~~receives~~ functions describing the local dynamics on the edges and nodes of ~~a~~ *the* graph `g` as inputs, and returns a composite function compatible with the DifferentialEquations.jl calling syntax.

```julia
nd = Network(g, vertex_dynamics,  edge_dynamics)
nd(dx, x, p, t)
```

Main features:
- Clear separation of local dynamics and topology: you can easily change *the* topology of your system or switch~~ing~~ out dynamic~~al~~ components.
- High performance when working with heterogeneous models *:* ~~(which means heaving~~ *you can have* different local dynamics in different parts of your network~~)~~.
- [Symbolic Indexing](@ref) into solutions and states: NetworkDynamics keeps track of the states of ~~the~~ *each* individual subsystem~~s~~.
- ~~Different~~ *Diverse* execution schemes: NetworkDynamics exploits the known inter-dependencies between components to auto parallelize execution, even on GPUs!
- Equation based models: ~~use [ModelingToolkit.jl](https://docs.sciml.ai/ModelingToolkit/dev/) to define local dynamics, use `NetworkDynamics.jl` to combine them into large networks!~~ *you can model local dynamics using [ModelingToolkit.jl](https://docs.sciml.ai/ModelingToolkit/dev/) and them combine them into larger networks by using `NetworkDynamics.jl`!*


## Where to begin?
~~Check out the [Mathematical Model](@ref) to understand the underlying modelling ideas of NetworkDynamics followed by the page on [Network Construction](@ref) to learn how to implement you own models.~~

*To learn how to implement your own models and understand the underlying modelling ideas of NetworkDynamics you should first read the [Mathematical Model](@ref) section, followed by section [Network Construction](@ref).*

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
