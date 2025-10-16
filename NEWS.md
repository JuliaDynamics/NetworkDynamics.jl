# NetworkDynamics Release Notes

## v0.10.11 Changelog
- [#324](https://github.com/JuliaDynamics/NetworkDynamics.jl/pull/324): Add custom ODEProblem constructor which takes a `NWState` object rather than flat arrays.
  Also always generate Network callbacks automatically.
  If you've previously passed `callback=get_callbacks(nw)`, you'll get a deprecation note. For any other usage of the `callback` keyword on ODEProblem Constructor you'll get an error.

## v0.10.10 Changelog
- [#323](https://github.com/JuliaDynamics/NetworkDynamics.jl/pull/323) Add `GuessFormula` system for improving initial guesses in component initialization:
  - New `GuessFormula` type and `@guessformula` macro for defining guess refinement formulas
  - GuessFormulas operate after InitFormulas in the initialization pipeline
  - Unlike InitFormulas (which set defaults), GuessFormulas refine initial guesses for free variables
  - Add `has_guessformula`, `get_guessformulas`, `set_guessformula!`, `add_guessformula!`, `delete_guessformulas!` metadata functions
  - Add `additional_guessformula` keyword to `initialize_component`, `initialize_component!`, and `initialize_componentwise` functions
  - Improved initialization documentation with execution order details

## v0.10.9 Changelog
- [#317](https://github.com/JuliaDynamics/NetworkDynamics.jl/pull/317) Enhanced callback system with negative affect support and runtime callback injection:
  - Add `affect_neg!` parameter to `ContinuousComponentCallback` and `VectorContinuousComponentCallback` for handling downcrossing events
  - Add ability to inject additional callbacks at runtime via `get_callbacks(nw, additional_callbacks)` without storing them in component metadata
  - Add `NetworkDynamics.pretty_f()` debugging utility for pretty-printing MTK-generated functions
  - Improve MTK integration warnings for nested event systems
  - Better initialization error messages with specific variable information when NaN values are detected
- [#314](https://github.com/JuliaDynamics/NetworkDynamics.jl/pull/314) Consolidate deprecated functionality: all deprecated functionality moved to dedicated `src/deprecated.jl` file for easier maintenance

## v0.10.8 Changelog
- [#313](https://github.com/JuliaDynamics/NetworkDynamics.jl/pull/313) Fix FilteringProxy issues and improvements:
  - Fix bug with empty `getindex` operations on FilteringProxy
  - Add pattern highlighting in FilteringProxy display (matches shown in light red)
  - Fix method ambiguity issues in symbolic indexing
  - Improve callback error handling with early validation for wrong signatures
  - Allow name clash reconstruction when it was previously allowed

## v0.10.7 Changelog (PR #312)
- [#312](https://github.com/JuliaDynamics/NetworkDynamics.jl/pull/312) Major rework of symbolic indexing system:
  - New index types: Added `ParamIdx` and `StateIdx` for explicit parameter vs state numeric indexing
  - Enhanced proxy system: Replaced internal `VProxy`/`EProxy` with new `FilteringProxy` system that provides more powerful interactive filtering and inspection capabilities
  - Improved display system: Major enhancements to `show` methods with compact printing and matched name highlighting
  - New exports: `ParamIdx`, `StateIdx`, `generate_indices`, `FilteringProxy`
  - All user-facing API remains backward compatible - `.v` and `.e` properties work as before but are now more powerful

## v0.10.6 Changelog
- [#309](https://github.com/JuliaDynamics/NetworkDynamics.jl/pull/309) improve error handling in initialization system:
  - Add custom exception types: `NetworkInitError` and `ComponentInitError` with detailed error messages
  - Enhanced error detection for NaN values, time-dependent systems, and RHS evaluation failures
  - Improved `find_fixpoint` function with better input validation and time handling support (fixes [#308](https://github.com/JuliaDynamics/NetworkDynamics.jl/pull/308))
  - Add equality (`==`) and approximate equality (`isapprox`) methods for `NWState` and `NWParameter`
  - Minor documentation fixes and spelling corrections
  - fix [#310](https://github.com/JuliaDynamics/NetworkDynamics.jl/pull/310) (allow `VPIndex(i)` to index into network objects)
  - implement [#307](https://github.com/JuliaDynamics/NetworkDynamics.jl/pull/307), allow nw[VIndex(:)], nw[[VIndex(1),VIndex(2)]]

## v0.10.4 Changelog
- [#303](https://github.com/JuliaDynamics/NetworkDynamics.jl/pull/303) update for ModelingToolkit.jl v10 compatibility:
  - rename all `ODESystem` -> `System` (follows MTK v10 API)
  - MTK extension now uses `mtkcompile` instead of `structural_simplify` internally
  - Add new `implicit_output` function to handle fully implicit output variables in MTK models
  - Add documentation for handling fully implicit outputs in MTK integration
  - Update minimum ModelingToolkit.jl requirement from v9.67 to v10

## v0.10.3 Changelog
- [#301](https://github.com/JuliaDynamics/NetworkDynamics.jl/pull/301) improve callback system performance and flexibility:
  - Add callback batching for better DiscreteComponentCallback performance
  - Allow `EIndex(1=>2)` as standalone edge index with relaxed type constraints
  - Optimize CallbackSet construction to prevent performance bottlenecks
  - Add important documentation warning about parameter array copying in callbacks
  - **Fixed spelling**: `ContinousComponentCallback` → `ContinuousComponentCallback` and `VectorContinousComponentCallback` → `VectorContinuousComponentCallback` (old names maintained as deprecated aliases for backward compatibility)

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
