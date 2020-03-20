# NetworkDynamics

A package for working with dynamical systems on complex networks. NetworkDynamics.jl provides an interface between [LightGraphs.jl](https://github.com/JuliaGraphs/LightGraphs.jl) and [DifferentialEquations.jl](https://github.com/JuliaDiffEq/DifferentialEquations.jl). It allows to define several types of dynamic and static nodes and edges and to link them up in order to create complex network dynamics.

At the moment the behaviour of a node or an edge can be described by algebraic equations or by ordinary differential equations (ODE). Support for stochastic differential equations (SDE) and delay differential equations (DDE) will be added in future releases.

## Installation

Installation is straightforward with Julia's package manager.

```julia-repl
(v1.3) pkg> add NetworkDynamics
```

## Multi-Threading

Since version `0.3.0` multi-threading via the `Threads.@threads` macro is enabled by default. This allows julia to integrate different nodes and edges in different threads, and leads to significant performance gains on parallel architectures.

However, it also causes an overhead in the order of *20 μs* per call to the ODE function which might impair performance on small networks (5 to 20 nodes and 10 to 50 edges) if only a single core is available to julia. To make multiple cores available the environment variable `JULIA_NUM_THREADS` has to be set **before** starting julia.

E.g. to start julia with 4 cores from a bash shell: `$ env JULIA_NUM_THREADS=4 julia`"

If you are using `Juno` for the `Atom` text editor `JULIA_NUM_THREADS` is set to the number of physical cores of your processor by default. This is also the number of threads we recommend to use.

For mor


https://docs.julialang.org/en/v1/manual/environment-variables/index.html#JULIA_NUM_THREADS-1




## PowerDynamics

[PowerDynamics.jl](https://juliaenergy.github.io/PowerDynamics.jl/stable/) is an open-source framework for dynamic power grid modeling and analysis build on top of NetworkDynamics.jl.


# Overview

The key construction is the function [`network_dynamics`](@ref) that takes in
two arrays of functions describing the local dynamics on the edges and nodes of
a graph `g`, and returns a composite function compatible with the
DifferentialEquations.jl calling syntax.

```julia
nd = network_dynamics(vertices!::Array{VertexFunction},  edges!::Array{EdgeFunction}, g)
nd(dx, x, p, t)
```



This page is still in development. :)
