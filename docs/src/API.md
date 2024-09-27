# API

The following functions are designed for public use

## Network Construction
```@docs
Network
dim(::Network)
pdim(::Network)
```

## Vertex Functions
```@docs
StaticVertex
ODEVertex
```

## Edge Functions
```@docs
StaticEdge
ODEEdge
```

### Coupling types
```@docs
Symmetric
AntiSymmetric
Directed
Fiducial
```

## Symbolic Indexing
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

### Index generators
```@docs
vidxs
eidxs
vpidxs
epidxs
```

### Symbolic Indices
```@docs
VIndex
EIndex
VPIndex
EPIndex
```

## Metadata API
### Component Accessors
```@docs
dim(::NetworkDynamics.ComponentFunction)
sym
pdim(::NetworkDynamics.ComponentFunction)
psym
obssym
depth
hasinputsym
inputsym
coupling
```
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

## Utils
### Callbacks
```@docs
save_parameters!
```
