# NetworkDynamics

The *NetworkDynamics.jl* package simulates the dynamics of complex networks. It provides an interface 
between the [Graphs.jl](https://github.com/JuliaGraphs/Graphs.jl) and the 
[DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl) packages and facilitates the simulation of 
highly efficient dynamic networks by describing the local dynamics of the edges and vertices of the network.

!!! note
    **Complex Networks in a glance**
    Complex network systems are composed by the entities that comprise them (the nodes) and the relationships that 
    connect each entity with one another (the edges). The mathematical structure is called 
    [Graph](https://en.wikipedia.org/wiki/Graph_theory) and it is used more or less interchangeably with the word 
    Network). The graphical depictions of such networks are also called graphs. You will see both usages in this guide.

In this graph of a simple 5-node network
```@example
using Graphs, GraphMakie, CairoMakie #hide
using GraphMakie.NetworkLayout #hide
CairoMakie.activate!(type="svg") #hide
g = smallgraph(:bull) #hide
fig, ax, p = graphplot(g; ilabels=["v$i" for i in 1:nv(g)], #hide
                          elabels=["e$i: $(e.src) ↦ $(e.dst)" for (i, e) in enumerate(edges(g))], #hide
                          layout=Align(Stress()), figure=(;resolution=(830,450))) #hide
ymin, ymax = extrema(last.(p[:node_pos][])) #hide
ylims!(ax, (ymin-0.11*(ymax-ymin), ymax+0.11*(ymax-ymin)))#hide
hidespines!(ax) #hide
hidedecorations!(ax) #hide
fig #hide
```
we can see that it is composed of nodes (v1 to v5) who are connected to each other. The lines connecting the nodes with 
each other ( e1: 1-->2, e2: 1-->3, e3: 2-->3, e4: 2-->4, e5: 3-->5) are called edges. Complex networks are composed of 
multiple nodes and edges, with most nodes connected to multiple other nodes with multiple edges.

(@Hans after rereading the text I realised that the information about the core of the package and the behaviours of the 
nodes and edges does not belong in the introduction but rather in the mathematical model, so I moved it. If you are ok
with this just delete this comment)

Main features:
- Clear separation of local dynamics and topology: you can easily change the topology of your system or switch out 
  dynamic components)
- High performance when working with heterogeneous models: you can have different local dynamics in different parts of 
  your network)
- [Symbolic Indexing](@ref) into solutions and states: NetworkDynamics keeps track of the states of each individual 
  subsystem.
- Diverse execution schemes: NetworkDynamics exploits the known interdependencies between components to auto 
  parallelize execution, even on GPUs!
- Equation based models: you can model local dynamics using 
  [ModelingToolkit.jl](https://docs.sciml.ai/ModelingToolkit/dev/) and them combine them into larger networks using 
  `NetworkDynamics.jl`!


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
