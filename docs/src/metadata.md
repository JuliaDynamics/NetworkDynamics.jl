# Metadata

Component function such as `ODEVertex`, `StaticEdge`, ... can store metadta. We distinguish between two kinds of metadata: component metadata and symbol metadata.

## Component Metadata
Component metadata is a `Dict{Symbol,Any}` attached to each component to store various information. Use [`metadata`](@ref) to retrieve the full dict.

To access the data, you can use the methods `has_metadata`, `get_metadata` and `set_metadata!` (see [Component Metadata API](@ref)).

Special uses: after [component wise initialization](@ref), the field `:init_residual` stores the residual vector of the nonlinear problem.

## Symbol Metadata
Each component stores symbol metadata. The symbol metadata is a `Dict{Symbol, Dict{Symbol, Any}}` which stores a metadate dict per symbol. Symbols are everything that appears in [`sym`](@ref), [`psym`](@ref), [`obssym`](@ref) and [`inputsym`](@ref).

To access the data, you can use the methods `has_metadata`, `get_metadata` and `set_metadata!` (see [Per Symbol Metadata API](@ref)).

Special cases for symbol metadata are:

- `default`: Stores default values for states/parameters. In initialization, those are considered fixed.
- `guess`: Stores a guess for a state/parameter which needs to solved during initialization ("free" variables).
- `bounds`: Store bounds for variables/parameters
- `init`: Stores the solution of the "free" variables during initialization.

Fore those, there are special functions `has_*`, `get_*` and `set_*!`. See [Per Symbol Metadata API](@ref).


Those are closely aligned to the [metadata use in ModelingToolkit](https://docs.sciml.ai/ModelingToolkit/stable/basics/Variable_metadata/). They are automaticially copied from the `ODESystem` if you use MTK models to create NetworkDynamic models.