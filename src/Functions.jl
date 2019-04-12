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
    massmatrix=nothing # Mass matrix for the equation
    sym=[:v for i in 1:dim] # Symbols for the dimensions
end


@Base.kwdef struct ODEEdge
    f! # The function with signature (dx, x, e_s, e_t, p, t) -> nothing
    dim # number of dimensions of x
    massmatrix=nothing # Mass matrix for the equation
    sym=[:e for i in 1:dim] # Symbols for the dimensions
end


@Base.kwdef struct DDEVertex
    f! # The function with signature (dv, v, h_v, e_s, e_t, h_e_s, h_e_d, p, t) -> nothing where h is the history function
    dim # number of dimensions
    massmatrix=nothing # Mass matrix for the equation
    sym=[:v for i in 1:dim] # Symbols for the dimensions
    tau_s=[] # Array of Delays for the incoming currents of different variables
    tau_d=[] # Array of Delays for the outgoing currents of different variables
end


@Base.kwdef struct DDEEdge
    f! # The function with signature (de, e, h_e, v_s, v_d, h_s, h_d, p, t)
    dim # number of variables
    massmatrix=nothing # Mass matrix for the equation
    sym=[:v for i in 1:dim] # Symbols for the dimensions
end

const VertexFunction = Union{ODEVertex, StaticVertex, DDEVertex}
const EdgeFunction = Union{ODEEdge, StaticEdge, DDEEdge}

#Experimental. Needs work.
function ODE_f!(f!)
    # If mass matrix = 0 the differential equation sets de = 0.
    # To set e to the value calculated by f! we first write the value calculated
    # by f! into de, the subtract e. This leads to the  constraint
    # 0 = - e + f(v_s, v_p, p, t)
    # where f(...) denotes the value that f!(a, ...) writes into a.
    func = (de,e,v_s,v_t,p,t) -> de .= f!(e,v_s,v_t,p,t) - e
    func
end
#DDE_f is used to transform StaticEdge and ODEEdge to a DDEEdge by changing its arguments accordingly.

function DDE_vertex_f!(f!)
    func = (dv,v,h_v,e_s,e_d,h_s,h_d,p,t) -> dv .= f!(dv,v,e_s,e_d,p,t)
    func
end

function DDEVertex(ov::ODEVertex)
    DDEVertex(DDE_vertex_f!(ov.f!),ov.dim,ov.massmatrix,ov.sym,[],[])
end

function DDEVertex(dv::DDEVertex)
    dv
end

function DDE_edge_f!(f!)
    func = (de,e,h_e,v_s,v_d,h_s,h_d,p,t) -> de .= f!(de,e,v_s,v_d,p,t)
    func
end

function DDEEdge(oe::ODEEdge)
    DDEEdge(DDE_edge_f!(oe.f!),oe.dim,oe.massmatrix,oe.sym)
end

function DDEEdge(de::DDEEdge)
    de
end

#This rightnow gives a Singular Exception Error because of the zero mass matrix. Easy workaround has yet to be found...
function ODEEdge(se::StaticEdge)
    ODEEdge(ODE_f!(se.f!), se.dim, sparse(0.0I,se.dim,se.dim), [:e for i in 1:se.dim])
end

# investigate whether convert(::Type(ODEEdge), se::StaticEdge) = ODEEdge(se)
# is useful, check out PowerDynamics for what it does...


end
