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
### Network State Object
```@docs
NWState
NWState(::Any)
NWState(::NWState)
NWState(::NWParameter)
NWState(::SciMLBase.DEIntegrator)
```

### Network Parameter Object
```@docs
NWParameter
NWParameter(::Any)
NWParameter(::NWParameter)
NWParameter(::SciMLBase.DEIntegrator)
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
