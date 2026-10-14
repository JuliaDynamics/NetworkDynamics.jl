# NetworkDynamics

NetworkDynamics.jl is a Julia package for modeling and simulating dynamical systems on networks.
The core idea is to split the **components** from the **topology**. Each node and each edge is a
small input-output system of differential and/or algebraic equations, and the graph defines how
they are connected. NetworkDynamics.jl assembles these components into one efficient system for
the solvers of [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl).

## Why NetworkDynamics.jl?

Many systems in science and engineering are networks of interacting units: power grids,
oscillator networks, gas and heat networks, mechanical structures or some models in
neuroscience. Realistic models of such systems are often heterogeneous, mix differential and
algebraic equations, and consist of thousands of components. Writing the right-hand side of such
a model by hand quickly becomes impractical. Fully symbolic tools treat the whole network as one
big system of equations, which can get challenging for large systems, because the clear structure
of the network is lost on the way.

NetworkDynamics.jl keeps the components separate from the topology. Instead of throwing that
structure away, it builds on it: components can be initialized one by one, states are indexed by
their position in the graph, and one component model is easily swapped for another. It is aimed
at researchers in complex systems, energy systems and related fields who need to simulate,
initialize and analyze such networks.

The main features are:

- Components can be written as plain Julia functions or as equation-based models in
  [ModelingToolkit.jl](https://docs.sciml.ai/ModelingToolkit/stable/), see [ModelingToolkit Integration](@ref).
- Networks can mix many different component types and are still evaluated fast and without
  allocations.
- States, parameters and observables can be accessed by name, see [Symbolic Indexing](@ref).
- Large networks can be initialized component by component, see [Initialization](@ref initialization-guide).
- Events and callbacks can be defined on the level of single components, see [Callbacks and Events](@ref Callbacks).
- The network can be evaluated on multiple threads or on the GPU.

## Where to begin?
To learn how to implement your own models and understand the underlying modeling ideas of NetworkDynamics you should
first read the [Mathematical Model](@ref) section, followed by the [Network Construction](@ref) section.

If you prefer to look at some concrete code first check out the [Getting Started](@ref) tutorial!


## Installation

1. Install Julia by following the [installation instructions](https://julialang.org/install/)
   for your operating system.
2. Install NetworkDynamics.jl with Julia's package manager:
   ```julia-repl
   pkg> add NetworkDynamics
   ```

If you are new to Julia, the PowerDynamics.jl documentation has a
[three-part introduction](https://juliaenergy.github.io/PowerDynamics.jl/stable/julia_setup/)
covering the setup, package environments and how to structure a research project. It is
written for PowerDynamics.jl, but most of it applies to NetworkDynamics.jl as well.
For more on how to use Julia, see [Modern Julia Workflows](https://modernjuliaworkflows.org/).


## Contributing and support

Questions, bug reports and contributions are welcome. See the
[contributing guide](contributing.md) for how to get help, report issues and contribute code, and
for who maintains the package.


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
as part of the *OpPoDyn*-Project (FKZ 03EI1071A/B, 2024-2027).

```@raw html
<img src="assets/bmwk_logo_en.svg" width="300"/>
```
