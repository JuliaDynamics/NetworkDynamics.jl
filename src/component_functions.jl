export ODEVertex, StaticEdge, ODEEdge
export Symmetric, AntiSymmetric, Directed, Fiducial

abstract type Coupling end
struct AntiSymmetric <: Coupling end
struct Symmetric <: Coupling end
struct Directed <: Coupling end
struct Fiducial <: Coupling end
const CouplingUnion = Union{AntiSymmetric,Symmetric,Directed,Fiducial}

abstract type ComponentFunction end
"""
Abstract supertype for all vertex functions.
"""
abstract type VertexFunction <: ComponentFunction end

"""
Abstract supertype for all edge functions.
"""
# abstract type EdgeFunction{C<:Coupling} <: ComponentFunction end
abstract type EdgeFunction{C} <: ComponentFunction end


coupling(::EdgeFunction{C}) where {C} = C()
coupling(::Type{<:EdgeFunction{C}}) where {C} = C()

Base.@kwdef struct ODEVertex{F} <: VertexFunction
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
    StaticEdge{typeof(f),typeof(coupling)}(f, dim, pdim)
end

struct ODEEdge{F,C} <: EdgeFunction{C}
    f::F # (de, e, v_s, v_t, p, t) -> nothing
    dim::Int # number of dimensions of e
    pdim::Int
end
function ODEEdge(; f, dim, pdim, coupling)
    # todo probably a method check?
    ODEEdge{typeof(f),typeof(coupling)}(f, dim, pdim)
end

dim(e::Union{EdgeFunction,VertexFunction}) = e.dim

pdim(e::Union{EdgeFunction,VertexFunction}) = e.pdim
aggrdepth(e::EdgeFunction) = e.dim
aggrdepth(e::EdgeFunction{Fiducial}) = floor(Int, e.dim / 2)

statetype(::T) where {T<:ComponentFunction} = statetype(T)
statetype(::Type{<:ODEVertex}) = Dynamic()
statetype(::Type{<:StaticEdge}) = Static()
statetype(::Type{<:ODEEdge}) = Dynamic()
