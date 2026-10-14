![NetworkDynamics Banner](./logo/banner-animated-dark.svg#gh-dark-mode-only)
![NetworkDynamics Banner](./logo/banner-animated.svg#gh-light-mode-only)

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliadynamics.github.io/NetworkDynamics.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliadynamics.github.io/NetworkDynamics.jl/dev/)
[![Testing](https://github.com/JuliaDynamics/NetworkDynamics.jl/actions/workflows/tests.yml/badge.svg?branch=main)](https://github.com/JuliaDynamics/NetworkDynamics.jl/actions/workflows/tests.yml)
[![codecov](https://codecov.io/gh/JuliaDynamics/NetworkDynamics.jl/graph/badge.svg?branch=main)](https://codecov.io/gh/JuliaDynamics/NetworkDynamics.jl)
[![DOI](https://zenodo.org/badge/169404414.svg)](https://doi.org/10.5281/zenodo.4396192)

NetworkDynamics.jl is a Julia package for modeling and simulating dynamical systems on networks.
The core idea is to split the **components** from the **topology**. Each node and each edge is a
small input-output system of differential and/or algebraic equations, and the graph defines how
they are connected. NetworkDynamics.jl assembles these components into one efficient system for
the solvers of [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl).

## Why NetworkDynamics.jl?

Many systems in science and engineering are networks of interacting units: power grids,
oscillator networks, gas and heat networks or mechanical structures. Realistic models of such
systems are often heterogeneous, mix differential and algebraic equations, and consist of
thousands of components. Writing the right-hand side of such a model by hand quickly becomes
impractical. Fully symbolic tools treat the whole network as one big system of equations, which
can get challenging for large systems, because the clear structure of the network is lost on the
way.

NetworkDynamics.jl keeps the components separate from the topology. Instead of throwing that
structure away, it builds on it: components can be initialized one by one, states are indexed by
their position in the graph, and it is easy to iterate over different dynamical realizations of
your components.

The main features are:

- Components can be written as plain Julia functions or as equation-based models in
  [ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl).
- Networks can mix many different component types and are still evaluated fast and without
  allocations.
- States, parameters and observables can be accessed by name, both in states and in solutions.
- Large networks can be initialized component by component.
- Events and callbacks can be defined on the level of single components.
- The network can be evaluated on multiple threads or on the GPU.

## Getting started

NetworkDynamics.jl is a registered Julia package and can be installed with
```julia-repl
pkg> add NetworkDynamics
```
As a small example, here is a network of 20 Kuramoto oscillators:
```julia
using NetworkDynamics, Graphs, OrdinaryDiffEqTsit5

# edge: coupling between the phases of its two end points
function kuramoto_edge!(e, θ_src, θ_dst, (K,), t)
    e[1] = K * sin(θ_src[1] - θ_dst[1])
end
edge = EdgeModel(; g=AntiSymmetric(kuramoto_edge!), outsym=[:P], psym=[:K=>3.0])

# vertex: phase oscillator driven by the sum of all edge outputs
function kuramoto_vertex!(dθ, θ, esum, (ω,), t)
    dθ[1] = ω + esum[1]
end
vertex = VertexModel(; f=kuramoto_vertex!, g=1, sym=[:θ], psym=[:ω])

nw = Network(watts_strogatz(20, 4, 0.2), vertex, edge)

s0 = NWState(nw)
s0.v[:, :θ] = 2π * rand(20)
s0.p.v[:, :ω] = randn(20)

prob = ODEProblem(nw, s0, (0.0, 10.0))
sol = solve(prob, Tsit5())
sol(10.0; idxs=VIndex(1:20, :θ))
```

The [documentation](https://juliadynamics.github.io/NetworkDynamics.jl/stable/) explains the
[mathematical model](https://juliadynamics.github.io/NetworkDynamics.jl/stable/mathematical_model/)
behind the package and contains several tutorials, starting with
[Getting Started](https://juliadynamics.github.io/NetworkDynamics.jl/stable/generated/getting_started_with_network_dynamics/).

## Contributing and support

Questions, bug reports and contributions are welcome. See [CONTRIBUTING.md](CONTRIBUTING.md) for
how to get help, report issues and contribute code, and for who maintains the package.

## PowerDynamics.jl

[PowerDynamics.jl](https://juliaenergy.github.io/PowerDynamics.jl/stable/) is an open-source framework for dynamic power grid modeling and analysis built on top of NetworkDynamics.jl.


## Citations

If you use NetworkDynamics.jl in your research publications, please cite our [paper](https://aip.scitation.org/doi/10.1063/5.0051387).

```latex
@article{NetworkDynamics.jl-2021,
	author = {Lindner, Michael and Lincoln, Lucas and Drauschke, Fenja and Koulen, Julia M. and Würfel, Hans and Plietzsch, Anton and Hellmann, Frank},
	doi = {10.1063/5.0051387},
	eprint = { https://doi.org/10.1063/5.0051387 },
	journal = {Chaos: An Interdisciplinary Journal of Nonlinear Science},
	number = {6},
	pages = {063133},
	title = {NetworkDynamics.jl—Composing and simulating complex networks in Julia},
	url = { https://doi.org/10.1063/5.0051387 },
	volume = {31},
	year = {2021}
}
```

## Funding
Development of this project was in part funded by the *German Federal Ministry for Economic Affairs and Climate Action* as part of the *OpPoDyn*-Project ([FKZ 03EI1071A/B](https://www.enargus.de/pub/bscw.cgi/?op=enargus.eps2&q=%2201258425/1%22), 2024-2027).

<img src="docs/src/assets/bmwk_logo_en.svg" width="300"/>
