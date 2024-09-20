# NetworkDynamics

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliadynamics.github.io/NetworkDynamics.jl/dev/)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliadynamics.github.io/NetworkDynamics.jl/stable)

A package for working with dynamical systems on complex networks. NetworkDynamics.jl provides an interface between [Graphs.jl](https://github.com/JuliaGraphs/Graphs.jl) and [DifferentialEquations.jl](https://github.com/JuliaDiffEq/DifferentialEquations.jl). It allows to define several types of dynamic and static nodes and edges and to link them up in order to create complex network dynamics.

The behavior of a node or an edge can be described by algebraic equations, by differential algebraic equation (DAEs) in mass matrix form or by ordinary differential equations (ODE). Stochastic ordinary differential equations (SDE) can be implemented as a [two-layer network](https://juliadynamics.github.io/NetworkDynamics.jl/dev/generated/StochasticSystem/). For details see the [docs](https://juliadynamics.github.io/NetworkDynamics.jl/dev/).

## Getting started

Check out the [getting started](https://juliadynamics.github.io/NetworkDynamics.jl/dev/generated/getting_started_with_network_dynamics/) example in our docs.

An [introductory talk](https://www.youtube.com/watch?v=GrmnbDYr6mM) was recorded at JuliaCon2020.

## Benchmarks

In our benchmark on the Kuramoto model NetworkDynamics.jl + DifferentialEquations.jl proved to be an especially performant solution, see https://github.com/PIK-ICoNe/NetworkDynamicsBenchmarks.

## PowerDynamics

[PowerDynamics.jl](https://juliaenergy.github.io/PowerDynamics.jl/stable/) is an open-source framework for dynamic power grid modeling and analysis build on top of NetworkDynamics.jl.

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
Development of this project was in part funded by the *German Federal Ministry for Economic Affairs and Climate Action* as part of the *OpPoDyn*-Project ([Project ID 01258425/1](https://www.enargus.de/pub/bscw.cgi/?op=enargus.eps2&q=%2201258425/1%22), 2024-2027).

<img src="docs/src/assets/bmwk_logo_en.svg" width="300"/>
