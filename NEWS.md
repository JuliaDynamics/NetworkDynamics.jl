# NetworkDynamics Release Notes

## v0.10.2 Changelog
- [#299](https://github.com/JuliaDynamics/NetworkDynamics.jl/pull/299) enhance metadata system with pattern matching and utility functions:
  - Add String/Regex pattern matching for all metadata functions (`has_metadata`, `get_metadata`, `set_metadata!`, etc.)
  - Add `strip_*!` functions to remove all metadata of a specific type from components
  - Add `free_u()` and `free_p()` functions to identify variables/parameters without default values
  - Support removing metadata by passing `nothing` or `missing` to `set_*!` functions

## v0.10.1 Changelog
- [#294](https://github.com/JuliaDynamics/NetworkDynamics.jl/pull/294) add linear stability analysis functions: `isfixpoint`, `jacobian_eigenvals`, and `is_linear_stable` with support for both ODE and DAE systems
- [#283](https://github.com/JuliaDynamics/NetworkDynamics.jl/pull/283) add automatic sparsity detection using `get_jac_prototype` and `set_jac_prototype!`
- [#285](https://github.com/JuliaDynamics/NetworkDynamics.jl/pull/285) rename `delete_initconstraint!` -> `delete_initconstaints!` and `delete_initformula!` -> `delete_initformulas!`

## v0.10 Changelog
- **BREAKING**: the interface initialization of components has changed: it is now split up in two versions, mutating and non mutating version. Also it errors now if the tolerance bounds are violated. See docs on initialization for more details.

- new `get_graph(::Network)` method to extract graph object from nw
- **improved Initialization System**: Added comprehensive initialization formulas and constraints system:
  - added `@initformula` to add explicit algebraic init equations for specific variables
  - added `@initconstraint` to add additional constraints for the component initialization
- allow access edges via Pairs, i.e. `EIndex(1=>2,:a)` references variable `:a` in edge from vertex 1 to 2. Works also with unique names of vertices like `EIndex(:a=>:b)` [#281](https://github.com/JuliaDynamics/NetworkDynamics.jl/pull/281).

## v0.9 Changelog
### Main changes in this release
NetworkDynamics v0.9 is a complete overhaul of the previous releases of NetworkDynamics.
Users of the package should probably read the new documentation carefully.

The most important changes are:

- Explicit split in `f` and `g` function: There is no split into `ODE` and
  `Static` components anymore, everything is unified in component models with
  internal function `f` and output function `g`.
- Parameters handling: Parameters are allways stored in a flat array. The
  Symbolic Indexing Interfaces helps to set and retrieve parameters.
- Automatic aggregation: vertices no longer receive a list of all connected
  edges. This lead to inhomogeneous call signatures and was a performance
  bottleneck. Now, each `Network` has a `aggregation` function attached to it.
  The backaned will perform a reduction over all connected edges to calculate
  the input for a certain vertex. In practice, for typical flow networks you'll
  allways receive the sum of all flows rather than the individual flows.
- Symbolic Indexing: the order of the states in the state vector changed in
  non-trivial ways. But the old `idx_containing` and `syms_containing` functions
  have been replaced with a much more capable symbolic indexing framework.

### Limitations
- We dropped support for delay differential equations. If you've been using that
  feature please reach out to us.
- Due to the built aggregation, the vertices cannot explicitly handle the inputs
  from edges differently anymore. If you've been relying on those features reach
  out to us.
