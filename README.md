# NetworkDynamics.jl

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliadynamics.github.io/NetworkDynamics.jl/dev/)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliadynamics.github.io/NetworkDynamics.jl/stable)
[![DOI](https://zenodo.org/badge/169404414.svg)](https://doi.org/10.5281/zenodo.4396192)

A package for working with dynamical systems on complex networks. NetworkDynamics.jl provides an interface between [Graphs.jl](https://github.com/JuliaGraphs/Graphs.jl) and [DifferentialEquations.jl](https://github.com/JuliaDiffEq/DifferentialEquations.jl).
It allows modeling of dynamical processes on networks using a *modular* approach: meaning that the overall network dynamics are composed based on dynamical models for the *nodes* and the *edges*.
The dynamical behavior of those components are described by differential algebraic equations with inputs and outputs.

Typical usecases for this modeling appraoch are diffusion processes, oscillator networks or power grids.

## Getting started

Check out our [documentation](https://juliadynamics.github.io/NetworkDynamics.jl/dev/).
A good place to start is the page on the [mathematical model](https://juliadynamics.github.io/NetworkDynamics.jl/dev/mathematical_model/). For a more hands-on approach check out the [getting started example](https://juliadynamics.github.io/NetworkDynamics.jl/dev/generated/getting_started_with_network_dynamics/).

An [introductory talk](https://www.youtube.com/watch?v=GrmnbDYr6mM) for an older version of this package was recorded for JuliaCon2020. Be aware that the API change significantly since then, but it can still give some insights into usecases of the package.

## Benchmarks

In our benchmark on the Kuramoto model NetworkDynamics.jl + DifferentialEquations.jl proved to be an especially performant solution, see https://github.com/PIK-ICoNe/NetworkDynamicsBenchmarks.

## PowerDynamics.jl

> [!IMPORTANT]
> As of 11/2024, the currently available version of PowerDynamics.jl is quite outdated and builds on an old version of NetworkDynamics. However, PowerDynamics will receive an Substantial update as part of an ongoing project. We expect a pre-release of that in the first half of 2025!

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
