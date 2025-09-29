# Metadata
Component models such as [`VertexModel`](@ref) and [`EdgeModel`](@ref) can store metadata. We distinguish between two 
kinds of metadata: component metadata and symbol metadata.

## Component Metadata
Component metadata is a `Dict{Symbol,Any}` attached to each component to store various information. Use 
[`metadata`](@ref) to retrieve the full dict.

To access the data, you can use the methods `has_metadata`, `get_metadata`, `set_metadata!` and `delete_metadata!` 
(see [Component Metadata API](@ref)).

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

For those, there are special functions `has_*`, `get_*`, `set_*!`, `delete_*!` and `strip_*!`. The `strip_*!` functions remove all metadata of a specific type from all symbols in a component. See [Per Symbol Metadata API](@ref).

These are closely aligned with the [metadata use in ModelingToolkit](@extref ModelingToolkit symbolic_metadata). They are automatically copied from the `System` if you use MTK models to create NetworkDynamics models.

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
