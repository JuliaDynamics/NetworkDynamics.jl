module NDFunctions

using LinearAlgebra
using SparseArrays

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

struct StaticVertex #???
    f! # ToDo
    dim # number of dimensions of x
    sym # Symbols for the dimensions
end

struct StaticEdge
    f! # (l, v_s, v_t, p, t) -> nothing
    dim # number of dimensions of x
end

struct ODEVertex
    f! # The function with signature (dx, x, e_s, e_t, p, t) -> nothing
    dim # number of dimensions of x
    massmatrix # Mass matrix for the equation
    sym # Symbols for the dimensions
end

function ODEVertex(f!, dim)
    ODEVertex(f!, dim, sparse(1.0I, dim, dim), nothing)
end

struct ODEEdge
    f! # The function with signature (dx, x, e_s, e_t, p, t) -> nothing
    dim # number of dimensions of x
    massmatrix # Mass matrix for the equation
    sym # Symbols for the dimensions
end

function ODEEdge(f!, dim)
    ODEEdge(f!, dim, sparse(1.0I, dim, dim), [:e for i in 1:dim])
end

struct DDEVertex
    f! # The function with signature (dv, v, h_v, e_s, e_t, h_e_s, h_e_d, p, t) -> nothing where h is the history function
    dim # number of dimensions
    massmatrix # Mass matrix for the equation
    sym # Symbols for the dimensions
    tau_s # Array of Delays for the incoming currents of different variables
    tau_d # Array of Delays for the outgoing currents of different variables
end

function DDEVertex(f!, dim, tau_s, tau_d)
    DDEVertex(f!, dim, sparse(1.0I, dim, dim), [:v for i in 1:dim], tau_s, tau_d)
end

struct DDEEdge
    f! # The function with signature (de, e, h_e, v_s, v_d, h_s, h_d, p, t)
    dim # number of variables
    massmatrix # Mass matrix for the equation
    sym # Symbols for the variables
end

function DDEEdge(f!,dim)
    DDEEdge(f!,dim,sparse(1.0I, dim, dim), [:e for i in 1:dim])
end

const VertexFunction = Union{ODEVertex, StaticVertex, DDEVertex}
const EdgeFunction = Union{ODEEdge, StaticEdge, DDEEdge}

function ODE_f(f!, de, e, v_s, v_t, p, t)
    # If mass matrix = 0 the differential equation sets de = 0.
    # To set e to the value calculated by f! we first write the value calculated
    # by f! into de, the subtract e. This leads to the  constraint
    # 0 = - e + f(v_s, v_p, p, t)
    # where f(...) denotes the value that f!(a, ...) writes into a.
    f!(de, v_s, v_t, p, t)
    de .-= e
    nothing
end

#DDE_f is used to transform StaticEdge and ODEEdge to an DDEEdge by changing its arguments accordingly.

function ODEEdge(se::StaticEdge)
    ODEEdge(
    ODE_f!(se.f!, de, e, v_s, v_t, p, t),
    massmatrix = sparse(0.0I,se.dim,se.dim),
    se.dim,
    [:e for i in 1:se.dim])
end

# investigate whether convert(::Type(ODEEdge), se::StaticEdge) = ODEEdge(se)
# is useful, check out PowerDynamics for what it does...


end
