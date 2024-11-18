# API

The following functions are designed for public use.

## Network Construction API
```@docs
Network
dim(::Network)
pdim(::Network)
```

## Component Functions
```@docs
VertexFunction
EdgeFunction
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
dim(::NetworkDynamics.ComponentFunction)
sym
outdim
outsym
pdim(::NetworkDynamics.ComponentFunction)
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
has_metadata(::NetworkDynamics.ComponentFunction, ::Symbol)
get_metadata(::NetworkDynamics.ComponentFunction, ::Symbol)
set_metadata!(::NetworkDynamics.ComponentFunction, ::Symbol, ::Any)
has_graphelement
get_graphelement
set_graphelement!
```
### Per-Symbol Metadata API
```@docs
symmetadata
get_metadata(::NetworkDynamics.ComponentFunction, ::Symbol, ::Symbol)
has_metadata(::NetworkDynamics.ComponentFunction, ::Symbol, ::Symbol)
set_metadata!(::NetworkDynamics.ComponentFunction, ::Symbol, ::Symbol, ::Any)
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

## Initialization
```@docs
find_fixpoint
initialize_component!
init_residual
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
Base.copy(::NetworkDynamics.ComponentFunction)
```
