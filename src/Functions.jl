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
export DDEVertex
export DDEEdge

#=All the structs contain the field sym, this has not been used yet in the main
constructor network_dynamics. Shouldn't be too hard though. More of a design question. StaticVertex
not implemented yet. =#

@Base.kwdef struct StaticVertex
    f! # (v, e_s, e_t, p, t) -> nothing
    dim # number of dimensions of x
    sym=[:v for i in 1:dim] # Symbols for the dimensions
end


@Base.kwdef struct StaticEdge
    f! # (e, v_s, v_t, p, t) -> nothing
    dim # number of dimensions of x
    sym=[:e for i in 1:dim] # Symbols for the dimensions
end


@Base.kwdef struct ODEVertex
    f! # The function with signature (dx, x, e_s, e_t, p, t) -> nothing
    dim # number of dimensions of x
    mass_matrix=I # Mass matrix for the equation
    sym=[:v for i in 1:dim] # Symbols for the dimensions
end


@Base.kwdef struct ODEEdge
    f! # The function with signature (dx, x, e_s, e_t, p, t) -> nothing
    dim # number of dimensions of x
    mass_matrix=I # Mass matrix for the equation
    sym=[:e for i in 1:dim] # Symbols for the dimensions
end


@Base.kwdef struct DDEVertex
    f! # The function with signature (dv, v, h_v, e_s, e_t, h_e_s, h_e_d, p, t) -> nothing where h is the history function
    dim # number of dimensions
    mass_matrix=I # Mass matrix for the equation
    sym=[:v for i in 1:dim] # Symbols for the dimensions
    tau_s=[] # Array of Delays for the incoming currents of different variables
    tau_d=[] # Array of Delays for the outgoing currents of different variables
end


@Base.kwdef struct DDEEdge
    f! # The function with signature (de, e, h_e, v_s, v_d, h_s, h_d, p, t)
    dim # number of variables
    mass_matrix=I # Mass matrix for the equation
    sym=[:v for i in 1:dim] # Symbols for the dimensions
end

const VertexFunction = Union{ODEVertex, StaticVertex, DDEVertex}
const EdgeFunction = Union{ODEEdge, StaticEdge, DDEEdge}

convert(::Type{ODEVertex}, x::StaticVertex) = ODEVertex(x)
promote_rule(::Type{ODEVertex}, ::Type{StaticVertex}) = ODEVertex

struct ODE_from_Static
    f!
end
function (ofs::ODE_from_Static)(dx,x,args...)
    # If mass matrix = 0 the differential equation sets dx = 0.
    # To set x to the value calculated by f! we first write the value calculated
    # by f! into dx, then subtract x. This leads to the  constraint
    # 0 = - x + f(...)
    # where f(...) denotes the value that f!(a, ...) writes into a.
    f!(dx,x,args...)
    dx .-= x
    nothing
end

function ODEVertex(sv::StaticVertex)
    ODEVertex(ODE_from_Static(sv.f!), sv.dim, 0., sv.sym)
end

function ODEEdge(se::StaticEdge)
    ODEEdge(ODE_from_Static(se.f!), se.dim, 0., se.sym)
end

# investigate whether convert(::Type(ODEEdge), se::StaticEdge) = ODEEdge(se)
# is useful, check out PowerDynamics for what it does...



end
