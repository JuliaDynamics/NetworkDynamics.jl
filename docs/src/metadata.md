# Metadata
Component models such as [`VertexModel`](@ref) and [`EdgeModel`](@ref) can store metadata. We distinguish between two kinds of metadata: component metadata and symbol metadata.

## Component Metadata
Component metadata is a `Dict{Symbol,Any}` attached to each component to store various information. Use [`metadata`](@ref) to retrieve the full dict.

To access the data, you can use the methods `has_metadata`, `get_metadata`, `set_metadata!` and `delete_metadata!` (see [Component Metadata API](@ref)).

Special metadata:

- `:graphelement`: optional field to specialize the graphelement for each
  component (`vidx`) for vertices, `(;src,dst)` named tuple of either vertex
  names or vertex indices for edges. Has special accessors `has_/get_/set_graphelement`.
- `:callback`: optional field to define callback functions on the component level. See [Callbacks](@ref) and [Callbacks API](@ref) for more information.
- `:position`: Store a tuple `(x, y)` with position of the node for plotting. Has special accessors `has_/get_/set_position`.
- `:marker`: Store a `Symbol` for the graph plot. Possible values could be `:circle`, `:rect`, `:utriangle`, `:cross`, `:diamond`, `:dtriangle`, `:pentagon`, `:xcross` or anything which works as a `marker` keyword argument in Makie.
- `:initconstraint`: Store additional initialization constrains. Has special `has_/get_/set_/delete_initconstraint` accessors. See [Initialization](@ref initialization-guide) for further
details.
- `:initformula`: Similar to initconstraint, but is a straight forward explicit mapping to initialize some variables.


## Symbol Metadata
Each component stores symbol metadata. The symbol metadata is a `Dict{Symbol, Dict{Symbol, Any}}` which stores a metadata dict per symbol. Symbols are everything that appears in [`sym`](@ref), [`psym`](@ref), [`obssym`](@ref) and [`insym`](@ref).

To access the data, you can use the methods `has_metadata`, `get_metadata`, `set_metadata!` and `delete_metadata!` (see [Per Symbol Metadata API](@ref)). These functions also support pattern matching using String or Regex patterns to match symbol names, making it easier to work with symbols containing special characters or when you only know part of the symbol name.

Special cases for symbol metadata are:

- `default`: Stores default values for states/parameters. In initialization, those are considered fixed.
- `guess`: Stores a guess for a state/parameter which needs to be solved during initialization ("free" variables).
- `bounds`: Stores bounds for variables/parameters
- `init`: Stores the solution of the "free" variables, this is rarely set manually but instead when calling [`initialize_component!`](@ref).
- `scope`: Marks the [`ParameterScope`](@ref) of a (parameter) variable as `:local`, `:component` or `:global`. See [Parameter Scope](@ref) below.

For those, there are special functions `has_*`, `get_*`, `set_*!`, `delete_*!` and `strip_*!`. The `strip_*!` functions remove all metadata of a specific type from all symbols in a component. See [Per Symbol Metadata API](@ref).

These are closely aligned with the [metadata use in ModelingToolkit](@extref ModelingToolkit symbolic_metadata). They are automatically copied from the `System` if you use MTK models to create NetworkDynamics models.

## Parameter Scope
The [`ParameterScope`](@ref) metadata (`scope`) declares on which level a parameter is expected to be consistent:

- `:local` (default): the parameter is purely local to the component and is not checked.
- `:component`: the parameter must be consistent **within a single** [`VertexModel`](@ref)/[`EdgeModel`](@ref) (e.g. across its subcomponents).
- `:global`: the parameter must be consistent **across the whole network**, i.e. all parameters sharing the same trailing symbol name (matching `r"NAME$"`) must hold the same value.

In `@mtkmodel`/`@component` style definitions the scope can be attached as variable metadata and is automatically copied to the component:

```julia
@parameters begin
    Sbase = 100, [scope = :global]
end
```

On the component level it can also be set manually using [`set_scope!`](@ref) (and queried via [`get_scope`](@ref)).

If the same trailing name is declared with **different** scopes (for example `:global` on one parameter and `:component` or `:local` on another sharing the name), each scope is checked independently — values are never compared across scopes — and a warning is emitted, since this is almost always a mistake.

Consistency of scoped parameters can be checked with [`chk_global_parameters`](@ref), which accepts a [`Network`](@ref), [`NWState`](@ref) or [`NWParameter`](@ref). For a `Network` the metadata **defaults** are compared, for `NWState`/`NWParameter` the **current values**. This check also runs automatically on `ODEProblem` construction and can be toggled via `NetworkDynamics.CHECK_GLOBAL_PARAMETERS[]`.

## Parameter Defaults From Another Parameter
The [`ParameterDefaultFrom`](@ref) metadata (`default_from`) lets a parameter take its **default value** from another parameter, so you only have to specify a value once. A common use is a base quantity like `Vbase` that several parts of a model should share. The source is named by its full parameter name (namespaced with `₊`), in one of two forms:

- a `Symbol` for a parameter elsewhere in the **same** [`VertexModel`](@ref)/[`EdgeModel`](@ref), e.g. `[default_from = :busbar₊Vbase]`, or
- a `Tuple` `(:src, name)` / `(:dst, name)` for a parameter of the **source/destination vertex** of an edge, e.g. `[default_from = (:src, :busbar₊Vbase)]`.

```julia
@parameters begin
    Vbase, [default_from = :busbar₊Vbase]         # same component
    Vbase, [default_from = (:src, :busbar₊Vbase)] # from the src vertex of an edge
end
```

It is resolved automatically when components and the [`Network`](@ref) are built (you can also trigger it manually with [`resolve_default_from!`](@ref)). Chains (`a` feeds `b` feeds `c`) propagate in a single pass. Same-component sources are resolved already in the [`VertexModel`](@ref)/[`EdgeModel`](@ref) constructor — for MTK-based components this is the point where the MTK system is compiled into a NetworkDynamics component (see [ModelingToolkit Integration](@ref)). Cross-component (`:src`/`:dst`) sources need the graph and are resolved in the `Network` constructor.

The behavior is:

- **A default you set manually is never overwritten.** It only fills in values it provided itself; a conflict with a manual default is silently kept (set `verbose=true` to get a warning).
- **A missing source always warns**, regardless of `verbose`, and the parameter is left untouched.
- **Misuse is an error:** `default_from` on something that is not a parameter, or an `:src`/`:dst` source on a vertex parameter.

!!! note "It copies a value, it does not link the model"
    This is **not** the same as writing `Vbase = busbar.Vbase` in the model equations. `default_from` simply copies a default value from a parameter you name; it does not create an equation or otherwise couple the two parameters. You assert that the source (e.g. `busbar₊Vbase`) exists — if it does not, you get a warning.

!!! warning "It only touches defaults"
    `default_from` sets parameter **defaults** during construction. If you later change a value on an [`NWState`](@ref)/[`NWParameter`](@ref), the dependent parameters are **not** updated to match — you have to keep them in sync yourself. [`chk_global_parameters`](@ref) can help catch resulting inconsistencies.

## Metadata Utils
Accessing metadata (especially defaults) of states and parameters is a very
common task. We provide several helper methods to do so. Please check out their docstrings for further explanation:

- [`dump_state`](@ref)
- [`dump_initial_state`](@ref)
- [`get_initial_state`](@ref)
- [`free_u`](@ref) - find variables without default values
- [`free_p`](@ref) - find parameters without default values
- [`describe_vertices`](@ref) (needs `DataFrames.jl` loaded)
- [`describe_edges`](@ref) (needs `DataFrames.jl` loaded)
