"""
Together with Constructors this module forms the backbone of the core API.
It provide the basic types to construct Arrays of VertexFunction and
EdgeFunction which can be handled by network_dynamics.
"""
module ComponentFunctions

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
export DDEVertex
export StaticDelayEdge
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

`dim` is the number of independent variables in the vertex equations and
`sym` is an array of symbols for these variables.

For more details see the documentation.
"""
@Base.kwdef struct StaticVertex{T} <: VertexFunction
    f!::T # (v, e_s, e_t, p, t) -> nothing
    dim::Int # number of dimensions of x
    sym=[:v for i in 1:dim] # Symbols for the dimensions
end

"""
    StaticEdge(f!, dim, sym)

Wrapper that ensures compatibility of a **mutating** function `f!` with
the key constructor `network_dynamics`.

`f!`  describes the local behaviour at a static edge and has to respect
the following calling syntax

```julia
f!(e, v_s, v_t, p, t) -> nothing
```

Here  `e`, `p` and `t` are the usual arguments, while
`v_s` and `v_d` are arrays containing the vertices which are
the source and destination of the described edge.

- `dim` is the number of independent variables in the edge equations and
- `sym` is an array of symbols for these variables.
- `coupling` is a Symbol describing if the EdgeFunction is intended for a directed graph (`:directed`) or for an undirected graph (`{:undirected, :symmetric, :antisymmetric, :fiducial}`). `:directed` is intended for directed graphs. `:undirected` is the default option and is only compatible with SimpleGraph. in this case f! should specify the coupling from a source vertex to a destination vertex. `:symmetric` and `:antisymmetric` trigger performance optimizations, if `f!` has that symmetry property. `:fiducial` lets the user specify both the coupling from src to dst, as well as the coupling from dst to src and is intended for advanced users.

For more details see the documentation.
"""
@Base.kwdef struct StaticEdge{T} <: EdgeFunction
    f!::T # (e, v_s, v_t, p, t) -> nothing
    dim::Int # number of dimensions of e
    coupling = :unspecified # :directed, :symmetric, :antisymmetric, :undirected
    sym=[:e for i in 1:dim] # Symbols for the dimensions

    function StaticEdge{T}(user_f!::T,
                           dim::Int,
                           coupling::Symbol,
                           sym::Vector{Symbol}) where T
        
        coupling_types = [:unspecified, :directed, :fiducial, :undirected, :symmetric,
                          :antisymmetric]

        coupling ∈ coupling_types ? nothing : error("Coupling type not recognized. Choose
                                                     from $coupling_types.")

        dim > 0 ? nothing : error("dim has to be a positive number.")

        dim == length(sym) ? nothing : error("Please specify a symbol for every dimension.")

        if coupling == :unspecified || :directed
            return new(user_f!, dim, coupling, sym)
        end

        dim % 2 == 0 ? nothing : error("Undirected edges are required to have even dim.")

        if coupling == :fiducial
            return new(user_f!, dim, coupling, sym)

        elseif coupling == :undirected
            # This might cause unexpected behaviour if source and destination vertex don't
            # have the same internal arguments.
            # Make sure to explicitly define the edge is :fiducial in that case.
            f! = @inline (e, v_s, v_d, p, t) -> begin
                @inbounds user_f!(view(e,1:dim), v_s, v_d, p, t)
                @inbounds user_f!(view(e,dim+1:2dim), v_d, v_s, p, t)
            end
        elseif coupling == :antisymmetric
            f! = @inline (e, v_s, v_d, p, t) -> begin
                @inbounds user_f!(view(e,1:dim), v_s, v_d, p, t)
                @inbounds for i in 1:dim
                    e[dim + i] = -1.0 * e[i]
                end
            end
        elseif coupling == :symmetric
            f! = @inline (e, v_s, v_d, p, t) -> begin
                @inbounds user_f!(view(e,1:dim), v_s, v_d, p, t)
                @inbounds for i in 1:dim
                    e[dim + i] = e[i]
                end
            end
        end
        # For edges with mass matrix this will be a little more complicated
        return new(f!, 2dim, coupling, repeat(sym, 2))
    end
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


"""
    DDEVertex(f!, dim, mass_matrix, sym)

Wrapper that ensures compatibility of a **mutating** function **`f!`** with
the key constructor `network_dynamics`.

**`f!`**  describes the local behaviour at a dynamic node and has to respect
the following calling syntax

```julia
f!(dv, v, e_s, e_t, h, p, t) -> nothing
```

Here `dv`, `v`, `p` and `t` are the usual ODE arguments, while
`e_s` and `e_d` are arrays containing the edges for which the
described vertex is the source or the destination respectively. `h` is the history array for `v`.

**`dim`** is the number of independent variables in the vertex equations and
**`sym`** is an array of symbols for these variables.
**`mass_matrix`** is an optional argument that defaults to the identity
matrix `I`. If a mass matrix M is given the system `M * dv = f!` will be
solved.

For more details see the documentation.
"""
@Base.kwdef struct DDEVertex{T} <: VertexFunction
    f!::T # The function with signature (dx, x, e_s, e_t, h, p, t) -> nothing
    dim::Int # number of dimensions of x
    mass_matrix=I # Mass matrix for the equation
    sym=[:v for i in 1:dim] # Symbols for the dimensions
end

function DDEVertex(ov::ODEVertex)
    f! = (dv, v, e_s, e_d, h_v, p, t) -> ov.f!(dv, v, e_s, e_d, p, t)
    DDEVertex(f!, ov.dim, ov.mass_matrix, ov.sym)
end

# function DDEVertex(sv::StaticVertex)
#     f! = (dv, v, e_s, e_d, h_v, p, t) -> sv.f!(v, e_s, e_d, p, t)
#     DDEVertex(f! = f!, dim = sv.dim, sym= sv.sym)
# end
function DDEVertex(sv::StaticVertex)
    f! = (dv, v, e_s, e_d, h_v, p, t) -> ODE_from_Static(sv.f!)(dv, v, e_s, e_d, p, t)
    DDEVertex(f!, sv.dim, 0., sv.sym)
end

# Promotion rules [eventually there might be too many types to hardcode everyhting]

convert(::Type{DDEVertex}, x::StaticVertex) = DDEVertex(x)
promote_rule(::Type{DDEVertex}, ::Type{StaticVertex}) = DDEVertex


convert(::Type{DDEVertex}, x::ODEVertex) = DDEVertex(x)
promote_rule(::Type{DDEVertex}, ::Type{ODEVertex}) = DDEVertex


"""
Like a static edge but with extra arguments for the history of the source and destination vertices. This is NOT a DDEEdge.
"""
 @Base.kwdef struct StaticDelayEdge{T} <: EdgeFunction
    f!::T # (e, v_s, v_t, p, t) -> nothing
    dim::Int # number of dimensions of x
    sym=[:e for i in 1:dim] # Symbols for the dimensions
end

function StaticDelayEdge(se::StaticEdge)
    f! = (e, v_s, v_d, h_v_s, h_v_d, p, t) -> se.f!(e, v_s, v_d, p, t)
    StaticDelayEdge(f!, se.dim, se.sym)
end

# Promotion rules

convert(::Type{StaticDelayEdge}, x::StaticEdge) = StaticDelayEdge(x)
promote_rule(::Type{StaticDelayEdge}, ::Type{StaticEdge}) = StaticDelayEdge


#= not used at the moment
"""
Like a ODEEdge but with extra arguments for the history of the source and destination vertices. This is NOT a DDEEdge.
"""
@Base.kwdef struct ODEDelayEdge{T} <: EdgeFunction
   f!::T # (e, v_s, v_t, p, t) -> nothing
   dim::Int # number of dimensions of x
   mass_matrix=I # Mass matrix for the equation
   sym=[:e for i in 1:dim] # Symbols for the dimensions
end

function ODEDelayEdge(se::ODEEdge)
   f! = (de, e, v_s, v_d, h_v_s, h_v_d, p, t) -> se.f!(de, e, v_s, v_d, p, t)
   ODEDelayEdge(f!, se.dim, se.sym)
end

# promotions rules are more complicated for this case

=#


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
