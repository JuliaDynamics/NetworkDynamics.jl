"""
Together with Constructors this module forms the backbone of the core API.
It provide the basic types to construct Arrays of VertexFunction and
EdgeFunction which can be handled by network_dynamics.
"""
module NDFunctions

using LinearAlgebra
# using SparseArrays

import Base.convert
import Base.promote_rule


export VertexFunction
export EdgeFunction
export StaticVertex
export StaticEdge
export ODEVertex
export ODEEdge
#export DDEVertex
#export DDEEdge
"""
Abstract supertype for all vertex functions.
"""
abstract type VertexFunction end
"""
Abstract supertype for all edge functions.
"""
abstract type EdgeFunction end

"""
    StaticVertex(f!, dim, sym)

Wrapper that ensures compatibility of a **mutating** function **`f!`** with
the key constructor `network_dynamics`.

**`f!`**  describes the local behaviour at a static node and has to respect
the following calling syntax

```julia
f!(v, e_s, e_t, p, t) -> nothing
```

Here  `v`, `p` and `t` are the usual arguments, while
`e_s` and `e_d` are arrays containing the edges for which the
described vertex is the source or the destination respectively.

**`dim`** is the number of independent variables in the vertex equations and
**`sym`** is an array of symbols for these variables.

For more details see the documentation.
"""
@Base.kwdef struct StaticVertex{T} <: VertexFunction
    f!::T # (v, e_s, e_t, p, t) -> nothing
    dim::Int # number of dimensions of x
    sym=[:v for i in 1:dim] # Symbols for the dimensions
end

"""
    StaticEdge(f!, dim, sym)

Wrapper that ensures compatibility of a **mutating** function **`f!`** with
the key constructor `network_dynamics`.

**`f!`**  describes the local behaviour at a static edge and has to respect
the following calling syntax

```julia
f!(e, v_s, v_t, p, t) -> nothing
```

Here  `e`, `p` and `t` are the usual arguments, while
`v_s` and `v_d` are arrays containing the vertices which are
the source and destination of the described edge.

**`dim`** is the number of independent variables in the edge equations and
**`sym`** is an array of symbols for these variables.

For more details see the documentation.
"""
@Base.kwdef struct StaticEdge{T} <: EdgeFunction
    f!::T # (e, v_s, v_t, p, t) -> nothing
    dim::Int # number of dimensions of x
    sym=[:e for i in 1:dim] # Symbols for the dimensions
end

"""
    ODEVertex(f!, dim, mass_matrix, sym)

Wrapper that ensures compatibility of a **mutating** function **`f!`** with
the key constructor `network_dynamics`.

**`f!`**  describes the local behaviour at a dynamic node and has to respect
the following calling syntax

```julia
f!(dv, v, e_s, e_t, p, t) -> nothing
```

Here `dv`, `v`, `p` and `t` are the usual ODE arguments, while
`e_s` and `e_d` are arrays containing the edges for which the
described vertex is the source or the destination respectively.

**`dim`** is the number of independent variables in the vertex equations and
**`sym`** is an array of symbols for these variables.
**`mass_matrix`** is an optional argument that defaults to the identity
matrix `I`. If a mass matrix M is given the system `M * dv = f!` will be
solved.

For more details see the documentation.
"""
@Base.kwdef struct ODEVertex{T} <: VertexFunction
    f!::T # The function with signature (dx, x, e_s, e_t, p, t) -> nothing
    dim::Int # number of dimensions of x
    mass_matrix=I # Mass matrix for the equation
    sym=[:v for i in 1:dim] # Symbols for the dimensions
end

"""
    ODEEdge(f!, dim, mass_matrix, sym)

Wrapper that ensures compatibility of a **mutating** function **`f!`** with
the key constructor `network_dynamics`.

**`f!`**  describes the local behaviour at a dynamic edge and has to respect
the following calling syntax

```julia
f!(de, e, v_s, v_t, p, t) -> nothing
```

Here  `de`, `e`, `p` and `t` are the usual arguments, while
`v_s` and `v_d` are arrays containing the vertices which are
the source and destination of the described edge.

**`dim`** is the number of independent variables in the edge equations and
**`sym`** is an array of symbols for these variables. For more details see
the documentation.
**`mass_matrix`** is an optional argument that defaults to the identity
matrix `I`. If a mass matrix M is given the system `M * de = f!` will be
solved.

For more details see the documentation.
"""
@Base.kwdef struct ODEEdge{T} <: EdgeFunction
    f!::T # The function with signature (dx, x, e_s, e_t, p, t) -> nothing
    dim::Int # number of dimensions of x
    mass_matrix=I # Mass matrix for the equation
    sym=[:e for i in 1:dim] # Symbols for the dimensions
end



convert(::Type{ODEVertex}, x::StaticVertex) = ODEVertex(x)
promote_rule(::Type{ODEVertex}, ::Type{StaticVertex}) = ODEVertex

# Not sure if the next line does something?
promote_rule(::Type{ODEVertex{T}}, ::Type{ODEVertex{U}}) where {T, U} = ODEVertex

convert(::Type{ODEEdge}, x::StaticEdge) = ODEEdge(x)
promote_rule(::Type{ODEEdge}, ::Type{StaticEdge}) = ODEEdge
"""
Promotes a StaticVertex to an ODEVertex.
"""
struct ODE_from_Static{T}
    f!::T
end
function (ofs::ODE_from_Static)(dx,x,args...)
    # If mass matrix = 0 the differential equation sets dx = 0.
    # To set x to the value calculated by f! we first write the value calculated
    # by f! into dx, then subtract x. This leads to the  constraint
    # 0 = - x + f(...)
    # where f(...) denotes the value that f!(a, ...) writes into a.
    ofs.f!(dx,args...)
    dx .-= x
    nothing
end

function ODEVertex(sv::StaticVertex)
    ODEVertex(ODE_from_Static(sv.f!), sv.dim, 0., sv.sym)
end

function ODEEdge(se::StaticEdge)
    ODEEdge(ODE_from_Static(se.f!), se.dim, 0., se.sym)
end


end
