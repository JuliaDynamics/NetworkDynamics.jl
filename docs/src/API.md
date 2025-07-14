# API

The following functions are designed for public use.

## Network Construction API
```@docs
Network
get_graph
dim(::Network)
pdim(::Network)
```

## Component Models
```@docs
VertexModel()
EdgeModel()
```

## Component Models with MTK
```@docs
VertexModel(::ModelingToolkit.ODESystem, ::Any, ::Any)
EdgeModel(::ModelingToolkit.ODESystem, ::Any, ::Any, ::Any, ::Any)
EdgeModel(::ModelingToolkit.ODESystem, ::Any, ::Any, ::Any)
```

### Output Function Helpers/Wrappers 
```@docs
StateMask
Symmetric
AntiSymmetric
Directed
Fiducial
```

### Accessors for Component Properties
```@docs
fftype
dim(::NetworkDynamics.ComponentModel)
sym
outdim
outsym
pdim(::NetworkDynamics.ComponentModel)
psym
obssym
hasinsym
insym
hasindim
indim
```

### `FeedForwardType`-Traits
```@docs
FeedForwardType
PureFeedForward
FeedForward
NoFeedForward
PureStateMap
```

## Symbolic Indexing API
### Network Parameter Object
```@docs
NWParameter
NWParameter(::Any)
NWParameter(::NWParameter)
NWParameter(::SciMLBase.DEIntegrator)
pflat
parameter_symbols
```

### Network State Object
```@docs
NWState
NWState(::Any)
NWState(::NWState)
NWState(::NWParameter)
NWState(::SciMLBase.DEIntegrator)
uflat
variable_symbols
```

### Symbolic Indices
```@docs
VIndex
EIndex
VPIndex
EPIndex
@obsex
```

### Index generators
```@docs
vidxs
eidxs
vpidxs
epidxs
```

## Metadata API
### Component Metadata API
```@docs
metadata
has_metadata(::NetworkDynamics.ComponentModel, ::Symbol)
get_metadata(::NetworkDynamics.ComponentModel, ::Symbol)
set_metadata!(::NetworkDynamics.ComponentModel, ::Symbol, ::Any)
delete_metadata!(::NetworkDynamics.ComponentModel, ::Symbol)
has_graphelement
get_graphelement
set_graphelement!
has_position
get_position
set_position!
has_marker
get_marker
set_marker!
```
### Per-Symbol Metadata API
```@docs
symmetadata
get_metadata(::NetworkDynamics.ComponentModel, ::Symbol, ::Symbol)
has_metadata(::NetworkDynamics.ComponentModel, ::Symbol, ::Symbol)
set_metadata!(::NetworkDynamics.ComponentModel, ::Symbol, ::Symbol, ::Any)
delete_metadata!(::NetworkDynamics.ComponentModel, ::Symbol, ::Symbol)
has_default
get_default
set_default!
delete_default!
has_guess
get_guess
set_guess!
delete_guess!
has_init
get_init
set_init!
delete_init!
has_bounds
get_bounds
set_bounds!
delete_bounds!
set_defaults!
set_interface_defaults!
get_defaults_dict
get_guesses_dict
get_bounds_dict
get_inits_dict
```

### Metadata and Inspection Utils
```@docs
dump_state
dump_initial_state
get_initial_state
describe_vertices
describe_edges
```

## Initialization
```@docs
find_fixpoint
initialize_componentwise
initialize_componentwise!
initialize_component
initialize_component!
init_residual
InitConstraint
@initconstraint
set_initconstraint!
delete_initconstraint!
has_initconstraint
get_initconstraints
add_initconstraint!
InitFormula
@initformula
has_initformula
get_initformulas
set_initformula!
add_initformula!
delete_initformula!
interface_values
```

## Callbacks API
### Define Callbacks
```@docs
NetworkDynamics.ComponentCallback
ContinousComponentCallback
VectorContinousComponentCallback
DiscreteComponentCallback
PresetTimeComponentCallback
ComponentCondition
ComponentAffect
SymbolicView
get_callbacks(::NetworkDynamics.Network)
```
### Attach Callbacks to Edge/VertexModels
```@docs
has_callback
get_callbacks(::NetworkDynamics.ComponentModel)
set_callback!
add_callback!
```

## Sparsity Detection
```@docs
get_jac_prototype
set_jac_prototype!
```

## Execution Types
```@docs
ExecutionStyle
SequentialExecution
PolyesterExecution
ThreadedExecution
KAExecution
```

## Aggregators
```@docs
Aggregator
SequentialAggregator
SparseAggregator
ThreadedAggregator
PolyesterAggregator
KAAggregator
```

## Utils
```@docs
save_parameters!
ff_to_constraint
Base.copy(::NetworkDynamics.ComponentModel)
extract_nw
```

## NetworkDynamicsInspector API
```@docs
inspect
dump_app_state
set_sol!
set_state!
set_graphplot!
set_timeseries!
define_timeseries!
```
