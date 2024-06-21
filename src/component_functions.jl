abstract type Coupling end
"""
    struct AntiSymmetric <: Coupling end

AntiSymmetric coupling type. The edge function f is evaluated once:

 - the dst vertex receives the first `d` values of the edge state,
 - the src vertex receives (-1) of that.

Here, `d` is the edge depth of the Network.
"""
struct AntiSymmetric <: Coupling end
"""
    struct Symmetric <: Coupling end

Symmetric coupling type. The edge function f is evaluated once:

 - the dst vertex receives the first `d` values of the edge state,
 - the src vertex receives the same.

Here, `d` is the edge depth of the Network.
"""
struct Symmetric <: Coupling end
"""
    struct Directed <: Coupling end

Directed coupling type. The edge function f is evaluated once:

 - the dst vertex receives the first `d` values of the edge state,
 - the src vertex receives nothing.

Here, `d` is the edge depth of the Network.
"""
struct Directed <: Coupling end
"""
    struct Fiducial <: Coupling end

Fiducial coupling type. The edge function f is evaluated once:

 - the dst vertex receives the `1:d` values of the edge state,
 - the src vertex receives the `d+1:2d` values of the edge state.

Here, `d` is the edge depth of the Network.
"""
struct Fiducial <: Coupling end
const CouplingUnion = Union{AntiSymmetric,Symmetric,Directed,Fiducial}

abstract type ComponentFunction end

Mixers.@pour CommonFields begin
    name::Symbol
    f::F
    dim::Int
    sym::Vector{Symbol}
    def::Vector{Union{Nothing,Float64}}
    depth::Int
    pdim::Int
    psym::Vector{Symbol}
    pdef::Vector{Union{Nothing,Float64}}
    obsf::OF
    obssym::Vector{Symbol}
end
compf(c::ComponentFunction) = c.f
dim(c::ComponentFunction)::Int = c.dim
sym(c::ComponentFunction)::Vector{Symbol} = c.sym
def(c::ComponentFunction)::Vector{Union{Nothing,Float64}} = c.def
pdim(c::ComponentFunction)::Int = c.pdim
psym(c::ComponentFunction)::Vector{Symbol} = c.psym
pdef(c::ComponentFunction)::Vector{Union{Nothing,Float64}} = c.pdef
obsf(c::ComponentFunction) = c.obsf
obssym(c::ComponentFunction)::Vector{Symbol} = c.obssym
depth(c::ComponentFunction)::Int = c.depth

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
struct ODEVertex{F,OF,MM} <: VertexFunction
    @CommonFields
    mass_matrix::MM
    # dfdp dfdv dfde
end
ODEVertex(; kwargs...) = _construct_comp(ODEVertex, kwargs)
ODEVertex(f; kwargs...) = ODEVertex(;f, kwargs...)
ODEVertex(f, dim; kwargs...) = ODEVertex(;f, _dimsym(dim)..., kwargs...)
ODEVertex(f, dim, pdim; kwargs...) = ODEVertex(;f, _dimsym(dim, pdim)..., kwargs...)

struct StaticVertex{F,OF} <: VertexFunction
    @CommonFields
end
StaticVertex(; kwargs...) = _construct_comp(StaticVertex, kwargs)
StaticVertex(f; kwargs...) = StaticVertex(;f, kwargs...)
StaticVertex(f, dim; kwargs...) = StaticVertex(;f, _dimsym(dim)..., kwargs...)
StaticVertex(f, dim, pdim; kwargs...) = StaticVertex(;f, _dimsym(dim, pdim)..., kwargs...)
function ODEVertex(sv::StaticVertex)
    d = Dict{Symbol,Any}()
    for prop in propertynames(sv)
        d[prop] = getproperty(sv, prop)
    end
    d[:f]  = let _f = sv.f
        (dx, x, esum, p, t) -> begin
            _f(dx, esum, p, t)
            @inbounds for i in eachindex(dx)
                dx[i] = dx[i] - x[i]
            end
            return nothing
        end
    end
    d[:mass_matrix] = 0.0
    ODEVertex(; d...)
end

struct StaticEdge{C,F,OF} <: EdgeFunction{C}
    @CommonFields
    coupling::C
end
StaticEdge(; kwargs...) = _construct_comp(StaticEdge, kwargs)
StaticEdge(f; kwargs...) = StaticEdge(;f, kwargs...)
StaticEdge(f, dim, coupling; kwargs...) = StaticEdge(;f, _dimsym(dim)..., coupling, kwargs...)
StaticEdge(f, dim, pdim, coupling; kwargs...) = StaticEdge(;f, _dimsym(dim, pdim)..., coupling, kwargs...)

struct ODEEdge{C,F,OF,MM} <: EdgeFunction{C}
    @CommonFields
    coupling::C
    mass_matrix::MM
end
ODEEdge(; kwargs...) = _construct_comp(ODEEdge, kwargs)
ODEEdge(f; kwargs...) = ODEEdge(;f, kwargs...)
ODEEdge(f, dim, coupling; kwargs...) = ODEEdge(;f, _dimsym(dim)..., coupling, kwargs...)
ODEEdge(f, dim, pdim, coupling; kwargs...) = ODEEdge(;f, _dimsym(dim, pdim)..., coupling, kwargs...)

statetype(::T) where {T<:ComponentFunction} = statetype(T)
statetype(::Type{<:ODEVertex}) = Dynamic()
statetype(::Type{<:StaticVertex}) = Static()
statetype(::Type{<:StaticEdge}) = Static()
statetype(::Type{<:ODEEdge}) = Dynamic()
isdynamic(::T) where {T<:ComponentFunction} = isdynamic(T)
isdynamic(x::Type{<:ComponentFunction}) = statetype(x) == Dynamic()

"""
    comptT(<:ComponentFunction) :: Type{<:ComponentFunction}

Returns the dispatch type of the component. Does not include unecessary type parameters.
"""
compT(::T) where {T<:ComponentFunction} = compT(T)
compT(::Type{<:ODEVertex}) = ODEVertex
compT(T::Type{<:StaticEdge}) = StaticEdge{typeof(coupling(T))}
compT(T::Type{<:ODEEdge}) = ODEEdge{typeof(coupling(T))}

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

# helper functions to dispatch on correct dim/sym keywords based on type
const _sym_T = Union{Vector, Pair, Symbol}
_dimsym(dim::Number) = (; dim)
_dimsym(sym::_sym_T) = (; sym)
_dimsym(dim::Number, pdim::Number) = (; dim, pdim)
_dimsym(dim::Number, psym::_sym_T) = (; dim, psym)
_dimsym(sym::_sym_T, pdim::Number) = (; sym, pdim)
_dimsym(sym::_sym_T, psym::_sym_T) = (; sym, psym)

"""
    _construct_comp(::Type{T}, kwargs) where {T}

Internal function to construct a component function from keyword arguments.
Fills up kw arguments with default values and performs sanity checks.
"""
function _construct_comp(::Type{T}, kwargs) where {T}
    dict = _fill_defaults(T, kwargs)

    # pop check keyword
    check = pop!(dict, :check, true)

    if !all(in(keys(dict)), fieldnames(T))
        throw(ArgumentError("Cannot construct $T: arguments $(setdiff(fieldnames(T), keys(dict))) missing."))
    end
    if !all(in(fieldnames(T)), keys(dict))
        throw(ArgumentError("Cannot construct $T: got additional arguments $(setdiff(keys(dict), fieldnames(T)))."))
    end

    args = map(fieldtypes(T), fieldnames(T)) do FT, name
        convert(FT, dict[name])
    end

    c = T(args...)
    check && chk_component(c)
    return c
end

"""
    _fill_defaults(T, kwargs)

Fill up keyword arguments `kwargs` for type T with default values.
Also perfoms sanity check some properties like mass matrix, depth, ...
"""
function _fill_defaults(T, kwargs)
    dict = Dict{Symbol, Any}(kwargs)

    # syms might be provided as single pairs or symbols, wrap in vector
    _maybewrap!(dict, :sym, Union{Symbol, Pair})
    _maybewrap!(dict, :psym, Union{Symbol, Pair})
    _maybewrap!(dict, :obssym, Symbol)

    # sym & dim
    if !haskey(dict, :dim)
        if haskey(dict, :sym)
            dict[:dim] = length(dict[:sym])
        else
            throw(ArgumentError("Either `dim` or `sym` must be provided to construct $T."))
        end
    end
    if haskey(dict, :sym)
        if _has_defaults(dict[:sym])
            if haskey(dict, :def)
                throw(ArgumentError("Provide defaults either as pairs in syms or as values in def!"))
            end
            dict[:sym], dict[:def] = _split_defaults(dict[:sym])
        end
    else
        if haskey(dict, :dim)
            dim = dict[:dim]
            if T <: VertexFunction
                dict[:sym] = [dim>1 ? Symbol("v", subscript(i)) : :s for i in 1:dict[:dim]]
            else
                dict[:sym] = [dim>1 ? Symbol("e", subscript(i)) : :e for i in 1:dict[:dim]]
            end
        else
            throw(ArgumentError("Either `dim` or `sym` must be provided to construct $T."))
        end
    end
    if !haskey(dict,:def)
        dict[:def] = Union{Float64,Nothing}[nothing for _ in 1:dict[:dim]]
    end
    @argcheck length(dict[:sym]) == length(dict[:def]) == dict[:dim] "Length of sym & def must match dim."

    # psym & pdim
    if !haskey(dict, :pdim)
        if haskey(dict, :psym)
            dict[:pdim] = length(dict[:psym])
        else
            dict[:pdim] = 0
        end
    end
    if haskey(dict, :psym)
        if _has_defaults(dict[:psym])
            if haskey(dict, :pdef)
                throw(ArgumentError("Provide defaults either as pairs in psyms or as values in pdef!"))
            end
            dict[:psym], dict[:pdef] = _split_defaults(dict[:psym])
        end
    else
        pdim = dict[:pdim]
        dict[:psym] = [pdim>1 ? Symbol("p", subscript(i)) : :p for i in 1:dict[:pdim]]
    end
    if !haskey(dict,:pdef)
        dict[:pdef] = Union{Float64,Nothing}[nothing for _ in 1:dict[:pdim]]
    end
    @argcheck length(dict[:psym]) == length(dict[:pdef]) == dict[:pdim] "Length of psym & pdef must match pdim."

    # obsf & obssym
    if !haskey(dict, :obsf)
        dict[:obsf] = nothing
    end
    if !haskey(dict, :obssym)
        if dict[:obsf] != nothing
            throw(ArgumentError("If `obsf` is provided, `obssym` must be provided as well."))
        else
            dict[:obssym] = Symbol[]
        end
    end
    if isnothing(dict[:obsf]) && !isempty(dict[:obssym])
        throw(ArgumentError("Cannot provide nonempty obssym without obsf."))
    end

    # name
    if !haskey(dict, :name)
        dict[:name] = _default_name(T)
    end

    # mass_matrix
    if isdynamic(T)
        if !haskey(dict, :mass_matrix)
            dict[:mass_matrix] = LinearAlgebra.I
        else
            mm = dict[:mass_matrix]
            if mm isa UniformScaling
            elseif mm isa Vector # convert to diagonal
                if length(mm) == dict[:dim]
                    dict[:mass_matrix] = LinearAlgebra.Diagonal(mm)
                else
                    throw(ArgumentError("If given as a vector, mass matrix must have length equal to dimension of component."))
                end
            elseif mm isa Number # convert to uniform scaling
                dict[:mass_matrix] = LinearAlgebra.UniformScaling(mm)
            elseif mm isa AbstractMatrix
                @argcheck size(mm) == (dict[:dim], dict[:dim]) "Size of mass matrix must match dimension of component."
            else
                throw(ArgumentError("Mass matrix must be a vector, square matrix,\
                                     a uniform scaling, or scalar. Got $(mm)."))
            end
        end
    end

    # coupling
    if T<:EdgeFunction && !haskey(dict, :coupling)
        throw(ArgumentError("Coupling type must be provided to construct $T."))
    end

    # depth
    if !haskey(dict, :depth)
        if T<:VertexFunction
            dict[:depth] = dict[:dim]
        elseif T<:EdgeFunction
            coupling = dict[:coupling]
            dim = dict[:dim]
            dict[:depth] = coupling==Fiducial() ? floor(Int, dim/2) : dim
        else
            throw(ArgumentError("Cannot construct $T: default depth not known."))
        end
    end
    if haskey(dict, :coupling) && dict[:coupling]==Fiducial() && dict[:depth] > floor(dict[:dim]/2)
        throw(ArgumentError("Depth cannot exceed half the dimension for Fiducial coupling."))
    elseif dict[:depth] > dict[:dim]
        throw(ArgumentError("Depth cannot exceed half the dimension."))
    end
    return dict
end

_default_name(::Type{StaticVertex}) = :StaticVertex
_default_name(::Type{ODEVertex}) = :ODEVertex
_default_name(::Type{StaticEdge}) = :StaticEdge
_default_name(::Type{ODEEdge}) = :ODEEdge

_has_defaults(vec::AbstractVector{<:Symbol}) = false
_has_defaults(vec::AbstractVector{<:Pair}) = true
_has_defaults(vec::AbstractVector) = any(el -> el isa Pair, vec)
function _split_defaults(input)
    Base.require_one_based_indexing(input)
    syms = Vector{Symbol}(undef, length(input))
    defs = Vector{Union{Nothing,Float64}}(undef, length(input))
    for i in eachindex(input)
        if input[i] isa Pair
            syms[i] = first(input[i])
            defs[i] = last(input[i])
        else
            syms[i] = input[i]
            defs[i] = nothing
        end
    end
    syms, defs
end

"If index `s` in `d` exists and isa `T` wrap in vector."
function _maybewrap!(d, s, T)
    if haskey(d, s)
        v = d[s]
        if v isa T
            d[s] = [v]
        end
    end
end
