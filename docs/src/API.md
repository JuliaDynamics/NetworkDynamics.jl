# API

The following functions are designed for public use.

## Network Construction API
```@docs
Network
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
```

### Network State Object
```@docs
NWState
NWState(::Any)
NWState(::NWState)
NWState(::NWParameter)
NWState(::SciMLBase.DEIntegrator)
uflat
pflat
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
has_graphelement
get_graphelement
set_graphelement!
```
### Per-Symbol Metadata API
```@docs
symmetadata
get_metadata(::NetworkDynamics.ComponentModel, ::Symbol, ::Symbol)
has_metadata(::NetworkDynamics.ComponentModel, ::Symbol, ::Symbol)
set_metadata!(::NetworkDynamics.ComponentModel, ::Symbol, ::Symbol, ::Any)
has_default
get_default
set_default!
has_guess
get_guess
set_guess!
has_init
get_init
set_init!
has_bounds
get_bounds
set_bounds!
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

## Initialization
```@docs
find_fixpoint
initialize_component!
init_residual
get_initial_state
dump_initial_state
set_defaults!
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
```
