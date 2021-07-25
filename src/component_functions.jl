export ODEVertex, StaticEdge, ODEEdge

"""
Abstract supertype for all vertex functions.
"""
abstract type VertexFunction end
"""
Abstract supertype for all edge functions.
"""
abstract type EdgeFunction end

@Base.kwdef struct ODEVertex{T} <: VertexFunction
    f::T
    dim::Int
    pdim::Int
end

@Base.kwdef struct StaticEdge{T} <: EdgeFunction
    f::T # (e, v_s, v_t, p, t) -> nothing
    dim::Int # number of dimensions of e
    pdim::Int
end

@Base.kwdef struct ODEEdge{T} <: EdgeFunction
    f::T # (de, e, v_s, v_t, p, t) -> nothing
    dim::Int # number of dimensions of e
    pdim::Int
end
