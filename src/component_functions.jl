export ODEVertex, StaticEdge, ODEEdge
export Symmetric, AntiSymmetric, Directed, Fiducial

abstract type Coupling end
struct AntiSymmetric <: Coupling end
struct Symmetric <: Coupling end
struct Directed <: Coupling end
struct Fiducial <: Coupling end
const CouplingUnion = Union{AntiSymmetric,Symmetric,Directed,Fiducial}

abstract type ComponentFunction end

Mixers.@premix struct Component{F, OF}
    """Component function"""
    f::F
    "bax"
    dim::Int
    "fbaxobar"
    sym::Vector{Symbol} = [Symbol("s$i") for i in 1:dim]
    def::Vector{Union{Nothing,Float64}} = [nothing for _ in 1:dim]
    pdim::Int;
    psym::Vector{Symbol} = [Symbol("p$i") for i in 1:pdim];
    pdef::Vector{Union{Nothing,Float64}} =[nothing for _ in 1:dim]
    obsf::OF = nothing;
    obssym::Vector{Symbol} = Symbol[];
    @assert dim == length(sym) "dim ≠ length(sym)"
    @assert pdim == length(psym) "pdim ≠ length(psym)"
end

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

"""
$(TYPEDEF)

# Fields
$(FIELDS)
"""
@Component @with_kw_noshow struct ODEVertex{MM} <: VertexFunction
    name::Symbol = :ODEVertex
    mass_matrix::MM = LinearAlgebra.I
    # dfdp dfdv dfde
end

@Component @with_kw_noshow struct StaticEdge{C} <: EdgeFunction{C}
    name::Symbol = :StaticEdge
    coupling::C
end

@Component @with_kw_noshow struct ODEEdge{C,MM} <: EdgeFunction{C}
    name::Symbol = :ODEEdge
    coupling::C
    mass_matrix::MM = LinearAlgebra.I
end

compf(c::ComponentFunction) = c.f
dim(c::ComponentFunction) = c.dim
pdim(c::ComponentFunction) = c.pdim


aggrdepth(e::EdgeFunction) = e.dim
aggrdepth(e::EdgeFunction{Fiducial}) = floor(Int, e.dim / 2)

statetype(::T) where {T<:ComponentFunction} = statetype(T)
statetype(::Type{<:ODEVertex}) = Dynamic()
statetype(::Type{<:StaticEdge}) = Static()
statetype(::Type{<:ODEEdge}) = Dynamic()


batchequal(a, b) = false
function batchequal(a::EdgeFunction, b::EdgeFunction)
    for f in (compf, dim, pdim, coupling)
        f(a) == f(b) || return false
    end
    return true
end
function batchequal(a::VertexFunction, b::VertexFunction)
    for f in (compf, dim, pdim)
        f(a) == f(b) || return false
    end
    return true
end
