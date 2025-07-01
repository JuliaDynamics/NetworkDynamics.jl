# NetworkDynamics

*NetworkDynamics.jl* is a package to simulate dynamical systems within complex networks. It provides an interface 
between the [Graphs.jl](https://github.com/JuliaGraphs/Graphs.jl) and the 
[DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl) packages and faciliates the simulation of 
highly efficient dynamic networks by describing the local dynamics on the edges and vertices of the graph.

!!! note
    Complex network systems are composed by the entities that comprise them (the nodes) and the relationships that connect
    each entity with one another (the edges). The graphical depictions of such networks are called graphs. The simplest
    The mathematical structure (used more or less interchangeably with Network) is also called [Graph](https://en.wikipedia.org/wiki/Graph_theory).
    The network (which can be seen in the figure below) is composed of two entities (so two nodes) who are only connected to each other.
    This connection between the two is the edge of the system. Complex networks are composed of multiple nodes and edges,
    with most nodes connected to multiple other nodes with multiple edges *(@Hans: can you created the graph of such a
    network and place in here?)*

```
    @example
    using Graphs, NetworkDynamics, OrdinaryDiffEqTsit5, StableRNGs, GraphMakie, Plots, CairoMakie #hide
    using GraphMakie.NetworkLayout #hide
    CairoMakie.activate!(type="svg") #hide
    g = smallgraph(:bull) #hide
    fig, ax, p = graphplot(g; ilabels=["v$i" for i in 1:nv(g)], #hide
                              elabels=["e$i: $(e.src) ↦ $(e.dst)" for (i, e) in enumerate(edges(g))], #hide
                              layout=Align(Stress()), figure=(;resolution=(400,200))) #hide
    ymin, ymax = extrema(last.(p[:node_pos][])) #hide
    ylims!(ax, (ymin-0.11*(ymax-ymin), ymax+0.11*(ymax-ymin)))#hide
    hidespines!(ax) #hide
    hidedecorations!(ax) #hide
    fig #hide
```

The behavior of a node or an edge can be described through the use of a) algebraic equations, b) differential algebraic 
equation (DAEs) in mass matrix form or c) ordinary differential equations (ODE). 

The core of the package is the function [`Network`](@ref). It accepts the functions describing the local dynamics on the
edges and nodes of the graph `g` as inputs, and returns a composite function compatible with the 
DifferentialEquations.jl calling syntax.

```julia
nd = Network(g, vertex_dynamics,  edge_dynamics)
nd(dx, x, p, t)
```

Main features:
- Clear separation of local dynamics and topology: you can easily change the topology of your system or switch out 
- dynamic components.
- High performance when working with heterogeneous models: you can have different local dynamics in different parts of 
- your network.
- [Symbolic Indexing](@ref) into solutions and states: NetworkDynamics keeps track of the states of each individual 
- subsystem.
- Diverse execution schemes: NetworkDynamics exploits the known inter-dependencies between components to auto 
- parallelize execution, even on GPUs!
- Equation based models: you can model local dynamics using 
- [ModelingToolkit.jl](https://docs.sciml.ai/ModelingToolkit/dev/) and them combine them into larger networks by using 
- `NetworkDynamics.jl`!


## Where to begin?
To learn how to implement your own models and understand the underlying modelling ideas of NetworkDynamics you should 
first read the [Mathematical Model](@ref) section, followed by section [Network Construction](@ref).

If you prefer to look at some concrete code first check out the [Getting Started](@ref) tutorial!


## Installation

1. Install Julia:
-   [Julia Installation](https://julialang.org/install/)
-   Find your OS and follow the instructions for the installation

2. Install NetworkDynamics.jl with Julia's package manager:
```julia-repl
(v1.11) pkg> add NetworkDynamics
```

3. Install the Julia package LiveServer:
```julia-repl
import Pkg; Pkg.add("LiveServer")
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
