module NDFunctions

using LinearAlgebra
using SparseArrays

import Base.convert
import Base.promote_rule

#=
This module needs to be documented. Together with Constructors
it forms the backbone of the core API, users should provide us
with Arrays of VertexFunction and EdgeFunction as well as a graph and that's it.
=#

export StaticVertex
export StaticEdge
export ODEVertex
export ODEEdge
export VertexFunction
export EdgeFunction
# export DDEVertex
# export DDEEdge

@Base.kwdef struct StaticVertex{T}
    f!::T # (v, e_s, e_t, p, t) -> nothing
    dim::Int # number of dimensions of x
    sym=[:v for i in 1:dim] # Symbols for the dimensions
end


@Base.kwdef struct StaticEdge{T}
    f!::T # (e, v_s, v_t, p, t) -> nothing
    dim::Int # number of dimensions of x
    sym=[:e for i in 1:dim] # Symbols for the dimensions
end


@Base.kwdef struct ODEVertex{T}
    f!::T # The function with signature (dx, x, e_s, e_t, p, t) -> nothing
    dim::Int # number of dimensions of x
    mass_matrix=I # Mass matrix for the equation
    sym=[:v for i in 1:dim] # Symbols for the dimensions
end


@Base.kwdef struct ODEEdge{T}
    f!::T # The function with signature (dx, x, e_s, e_t, p, t) -> nothing
    dim::Int # number of dimensions of x
    mass_matrix=I # Mass matrix for the equation
    sym=[:e for i in 1:dim] # Symbols for the dimensions
end


const VertexFunction = Union{ODEVertex, StaticVertex}
const EdgeFunction = Union{ODEEdge, StaticEdge}

convert(::Type{ODEVertex}, x::StaticVertex) = ODEVertex(x)
promote_rule(::Type{ODEVertex}, ::Type{StaticVertex}) = ODEVertex

# Not sure if the next line does something?
promote_rule(::Type{ODEVertex{T}}, ::Type{ODEVertex{U}}) where {T, U} = ODEVertex

convert(::Type{ODEEdge}, x::StaticEdge) = ODEEdge(x)
promote_rule(::Type{ODEEdge}, ::Type{StaticEdge}) = ODEEdge

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
