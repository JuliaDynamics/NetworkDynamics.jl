# NetworkDynamics

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://pik-icone.github.io/NetworkDynamics.jl/dev)
<!-- [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://pik-icone.github.io/NetworkDynamics.jl/stable) -->

A package for working with dynamical systems on complex networks. NetworkDynamics.jl provides an interface between [LightGraphs.jl](https://github.com/JuliaGraphs/LightGraphs.jl) and [DifferentialEquations.jl](https://github.com/JuliaDiffEq/DifferentialEquations.jl). It allows to define several types of dynamic and static nodes and edges and to link them up in order to create complex network dynamics.

The behavior of a node or an edge can be described by algebraic equations, by differential algebraic equation (DAEs) in mass matrix form, by ordinary differential equations (ODE) or by delay differential equations (DDE). Stochastic ordinary differential equations (SDE) can be implemented as a [two-layer network](https://github.com/pik-icone/NetworkDynamics.jl/blob/master/examples/sde.jl). For details see the [docs](https://pik-icone.github.io/NetworkDynamics.jl/dev).


## Getting started

Check out our step-by-step tutorial as a [jupyter notebook](https://github.com/pik-icone/NetworkDynamics.jl/blob/master/examples/getting_started_with_network_dynamics.ipynb) or [in the docs](https://pik-icone.github.io/NetworkDynamics.jl/dev/getting_started_with_network_dynamics/).

An [introductory talk](https://www.youtube.com/watch?v=GrmnbDYr6mM) was recorded at JuliaCon2020.

## Benchmarks

In our benchmark on the Kuramoto model NetworkDynamics.jl + DifferentialEquations.jl proved to be an especially performant solution, see https://github.com/PIK-ICoNe/NetworkDynamicsBenchmarks.

## PowerDynamics

[PowerDynamics.jl](https://juliaenergy.github.io/PowerDynamics.jl/stable/) is an open-source framework for dynamic power grid modeling and analysis build on top of NetworkDynamics.jl. 

## Citations

If you use NetworkDynamics.jl in your research publications, please cite our [paper](https://arxiv.org/abs/2012.12696).

```latex
@article{lindner2020networkdynamics,
  title={NetworkDynamics. jl--Composing and simulating complex networks in Julia},
  author={Lindner, Michael and Lincoln, Lucas and Drauschke, Fenja and Koulen, Julia Monika and W{\"u}rfel, Hans and Plietzsch, Anton and Hellmann, Frank},
  journal={arXiv preprint arXiv:2012.12696},
  year={2020}
}
```

## Old Documentation

Documentation for relases prior to version 0.5.3 can be found [here](https://pik-icone.github.io/NetworkDynamicsDocumentationHistory/). Current docs are [here](https://pik-icone.github.io/NetworkDynamics.jl/).
