export ODEVertex, StaticEdge, ODEEdge
export Symmetric, AntiSymmetric, Directed, Fiducial

"""
Abstract supertype for all vertex functions.
"""
abstract type VertexFunction end

abstract type Coupling end

"""
Abstract supertype for all edge functions.
"""
abstract type EdgeFunction{C<:Coupling} end

struct AntiSymmetric <: Coupling end
struct Symmetric <: Coupling end
struct Directed <: Coupling end
struct Fiducial <: Coupling end

coupling(::EdgeFunction{C}) where C = C()
coupling(::Type{<:EdgeFunction{C}}) where C = C()

@Base.kwdef struct ODEVertex{F} <: VertexFunction
    f::F
    dim::Int
    pdim::Int
end

struct StaticEdge{F,C} <: EdgeFunction{C}
    f::F # (e, v_s, v_t, p, F) -> nothing
    dim::Int # number of dimensions of e
    pdim::Int
end
function StaticEdge(; f, dim, pdim, coupling)
    # todo probably a method check?
    StaticEdge{typeof(f), typeof(coupling)}(f, dim, pdim)
end

struct ODEEdge{F,C} <: EdgeFunction{C}
    f::F # (de, e, v_s, v_t, p, t) -> nothing
    dim::Int # number of dimensions of e
    pdim::Int
end
function ODEEdge(; f, dim, pdim, coupling)
    # todo probably a method check?
    StaticEdge{typeof(f), typeof(coupling)}(f, dim, pdim)
end
