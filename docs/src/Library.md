# API

The following functions are designed for public use

```@index
```

## Network Construction
```@docs
Network
dim
pdim
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

## Utils
### Callbacks
```@docs
save_parameters!
```
