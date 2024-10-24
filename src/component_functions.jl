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

"""
Abstract supertype for all vertex functions.
"""
abstract type VertexFunction <: ComponentFunction end

"""
Abstract supertype for all edge functions.
"""
# abstract type EdgeFunction{C<:Coupling} <: ComponentFunction end
abstract type EdgeFunction{C} <: ComponentFunction end

Mixers.@pour CommonFields begin
    name::Symbol
    f::F
    sym::Vector{Symbol}
    depth::Int
    psym::Vector{Symbol}
    inputsym:: Union{Nothing, Vector{Symbol}, @NamedTuple{src::Vector{Symbol},dst::Vector{Symbol}}}
    obsf::OF
    obssym::Vector{Symbol}
    symmetadata::Dict{Symbol,Dict{Symbol, Any}}
    metadata::Dict{Symbol,Any}
end

"""
    compf(c::ComponentFunction)

Retrieve the actual dynamical function the component.
"""
compf(c::ComponentFunction) = c.f

"""
    dim(c::ComponentFunction)::Int

Retrieve the dimension of the component.
"""
dim(c::ComponentFunction)::Int = length(sym(c))

"""
    sym(c::ComponentFunction)::Vector{Symbol}

Retrieve the symbols of the component.
"""
sym(c::ComponentFunction)::Vector{Symbol} = c.sym

"""
    pdim(c::ComponentFunction)::Int

Retrieve the parameter dimension of the component.
"""
pdim(c::ComponentFunction)::Int = length(psym(c))

"""
    psym(c::ComponentFunction)::Vector{Symbol}

Retrieve the parameter symbols of the component.
"""
psym(c::ComponentFunction)::Vector{Symbol} = c.psym

"""
    obsf(c::ComponentFunction)

Retrieve the observation function of the component.
"""
obsf(c::ComponentFunction) = c.obsf

"""
    obssym(c::ComponentFunction)::Vector{Symbol}

Retrieve the observation symbols of the component.
"""
obssym(c::ComponentFunction)::Vector{Symbol} = c.obssym

"""
    depth(c::ComponentFunction)::Int

Retrieve the depth of the component. The depth is the number of "outputs". For
example, in a vertex function, a depth `N` means, that the connected edges receive
the states `1:N`.
"""
depth(c::ComponentFunction)::Int = c.depth

"""
    symmetadata(c::ComponentFunction)::Dict{Symbol,Dict{Symbol,Any}}

Retrieve the metadata dictionary for the symbols. Keys are the names of the
symbols as they appear in [`sym`](@ref), [`psym`](@ref), [`obssym`](@ref) and [`inputsym`](@ref).

See also [`symmetadata`](@ref)
"""
symmetadata(c::ComponentFunction)::Dict{Symbol,Dict{Symbol,Any}} = c.symmetadata

"""
    metadata(c::ComponentFunction)

Retrieve metadata object for the component.

See also [`metadata`](@ref)
"""
metadata(c::ComponentFunction)::Dict{Symbol,Any} = c.metadata

"""
    hasinputsym(c::ComponentFunction)

Checks if the optioan field `inputsym` is present in the component function.
"""
hasinputsym(c::ComponentFunction) = !isnothing(c.inputsym)

"""
    hasinputsym(c::VertexFunction)::Vector{Symbol}
    hasinputsym(c::EdgeFunction)::NamedTuple with :src and :dst keys

Musst be called *after* [`hasinputsym`](@ref) returned true.
Gives the `inputsym` vector(s). For vertex functions just a single vector, for
edges it returns a named tuple `(; src, dst)` with two symbol vectors.
"""
inputsym(c::VertexFunction)::Vector{Symbol} = c.inputsym
inputsym(c::EdgeFunction)::@NamedTuple{src::Vector{Symbol},dst::Vector{Symbol}} = c.inputsym

"""
    coupling(::EdgeFunction)
    coupling(::Type{<:EdgeFunction})

Returns the coupling of the given `EdgeFunction`.
"""
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
ODEVertex(v::ODEVertex; kwargs...) = _reconstruct_comp(ODEVertex, v, kwargs)

struct StaticVertex{F,OF} <: VertexFunction
    @CommonFields
end
StaticVertex(; kwargs...) = _construct_comp(StaticVertex, kwargs)
StaticVertex(f; kwargs...) = StaticVertex(;f, kwargs...)
StaticVertex(f, dim; kwargs...) = StaticVertex(;f, _dimsym(dim)..., kwargs...)
StaticVertex(f, dim, pdim; kwargs...) = StaticVertex(;f, _dimsym(dim, pdim)..., kwargs...)
StaticVertex(v::StaticVertex; kwargs...) = _reconstruct_comp(StaticVertex, v, kwargs)
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
StaticEdge(e::StaticEdge; kwargs...) = _reconstruct_comp(StaticEdge, e, kwargs)

struct ODEEdge{C,F,OF,MM} <: EdgeFunction{C}
    @CommonFields
    coupling::C
    mass_matrix::MM
end
ODEEdge(; kwargs...) = _construct_comp(ODEEdge, kwargs)
ODEEdge(f; kwargs...) = ODEEdge(;f, kwargs...)
ODEEdge(f, dim, coupling; kwargs...) = ODEEdge(;f, _dimsym(dim)..., coupling, kwargs...)
ODEEdge(f, dim, pdim, coupling; kwargs...) = ODEEdge(;f, _dimsym(dim, pdim)..., coupling, kwargs...)
ODEEdge(e::ODEEdge; kwargs...) = _reconstruct_comp(ODEEdge, e, kwargs)

statetype(::T) where {T<:ComponentFunction} = statetype(T)
statetype(::Type{<:ODEVertex}) = Dynamic()
statetype(::Type{<:StaticVertex}) = Static()
statetype(::Type{<:StaticEdge}) = Static()
statetype(::Type{<:ODEEdge}) = Dynamic()

isdynamic(x) = statetype(x) == Dynamic()
isstatic(x)  = statetype(x) == Static()

"""
    dispatchT(<:ComponentFunction) :: Type{<:ComponentFunction}

Returns the type "essence" of the component used for dispatch.
Fills up type parameters with `nothing` to ensure `Core.compiler.isconstType`
for GPU compatibility.
"""
dispatchT(::T) where {T<:ComponentFunction} = dispatchT(T)
dispatchT(::Type{<:StaticVertex}) = StaticVertex{nothing,nothing}
dispatchT(::Type{<:ODEVertex}) = ODEVertex{nothing,nothing,nothing}
dispatchT(T::Type{<:StaticEdge}) = StaticEdge{typeof(coupling(T)),nothing,nothing}
dispatchT(T::Type{<:ODEEdge}) = ODEEdge{typeof(coupling(T)),nothing,nothing,nothing}

batchequal(a, b) = false
function batchequal(a::EdgeFunction, b::EdgeFunction)
    compf(a) == compf(b) || return false
    dim(a) == dim(b) || return false
    pdim(a) == pdim(b) || return false
    coupling(a) == coupling(b) || return false
    return true
end
function batchequal(a::VertexFunction, b::VertexFunction)
    compf(a) == compf(b) || return false
    dim(a) == dim(b) || return false
    pdim(a) == pdim(b) || return false
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

    # check signature of f
    # if !_valid_signature(T, dict[:f])
    #     throw(ArgumentError("Function f does not take the correct number of arguments."))
    # end

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

function _reconstruct_comp(::Type{T}, cf::ComponentFunction, kwargs) where {T}
    fields = fieldnames(T)
    dict = Dict{Symbol, Any}()
    for f in fields
        dict[f] = getproperty(cf, f)
    end
    for (k, v) in kwargs
        dict[k] = v
    end
    _construct_comp(T, dict)
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

    symmetadata = get!(dict, :symmetadata, Dict{Symbol,Dict{Symbol,Any}}())

    metadata = try
        convert(Dict{Symbol,Any}, get!(dict, :metadata, Dict{Symbol,Any}()))
    catch e
        throw(ArgumentError("Provided metadata keyword musst be a Dict{Symbol,Any}. Got $(repr(dict[:metadata]))."))
    end

    if haskey(dict, :graphelement)
        ge = pop!(dict, :graphelement)
        metadata[:graphelement] = ge
    end
    if haskey(dict, :vidx) && T <: VertexFunction
        vidx = pop!(dict, :vidx)
        metadata[:graphelement] = vidx
    end
    if haskey(dict, :src) && haskey(dict, :dst) && T <: EdgeFunction
        src = pop!(dict, :src)
        dst = pop!(dict, :dst)
        metadata[:graphelement] = (; src, dst)
    end

    # sym & dim
    haskey(dict, :dim) || haskey(dict, :sym) || throw(ArgumentError("Either `dim` or `sym` must be provided to construct $T."))
    if haskey(dict, :sym)
        if haskey(dict, :dim)
            if dict[:dim] != length(dict[:sym])
                throw(ArgumentError("Length of sym and dim must match."))
            end
            # @warn "Unnecessary kw dim, can be infered from sym."
            delete!(dict, :dim)
        end
        if _has_metadata(dict[:sym])
            dict[:sym], _metadata = _split_metadata(dict[:sym])
            mergewith!(merge!, symmetadata, _metadata)
        end
    else
        _dim = pop!(dict, :dim)
        if T <: VertexFunction
            dict[:sym] = [_dim>1 ? Symbol("v", subscript(i)) : :s for i in 1:_dim]
        else
            dict[:sym] = [_dim>1 ? Symbol("e", subscript(i)) : :e for i in 1:_dim]
        end
    end
    dim = length(dict[:sym])
    if haskey(dict,:def)
        _def = pop!(dict, :def)
        @argcheck length(_def)==dim  "Length of sym & def must match dim."
        for (sym, def) in zip(dict[:sym], _def)
            if isnothing(def)
                continue
            end
            if haskey(symmetadata, sym) && haskey(symmetadata[sym], :default)
                throw(ArgumentError("Default value for $sym is already provided in metadata."))
            else
                mt = get!(symmetadata, sym, Dict{Symbol,Any}())
                mt[:default] = def
            end
        end
    end

    # psym & pdim
    if !haskey(dict, :pdim) && !haskey(dict, :psym)
        dict[:pdim] = 0
    end
    if haskey(dict, :psym)
        if haskey(dict, :pdim)
            if dict[:pdim] != length(dict[:psym])
                throw(ArgumentError("Length of sym and dim must match."))
            end
            # @warn "Unnecessary kw pdim, can be infered from psym."
            delete!(dict, :pdim)
        end
        if _has_metadata(dict[:psym])
            dict[:psym], _metadata = _split_metadata(dict[:psym])
            mergewith!(merge!, symmetadata, _metadata)
        end
    else
        _pdim = pop!(dict, :pdim)
        dict[:psym] = [_pdim>1 ? Symbol("p", subscript(i)) : :p for i in 1:_pdim]
    end
    if haskey(dict,:pdef)
        _pdef = pop!(dict, :pdef)
        @argcheck length(_pdef) == length(dict[:psym]) "Length of sym & def must match dim."
        for (sym, def) in zip(dict[:psym], _pdef)
            if isnothing(def)
                continue
            end
            if haskey(symmetadata, sym) && haskey(symmetadata[sym], :default)
                throw(ArgumentError("Default value for $sym is already provided in metadata."))
            else
                mt = get!(symmetadata, sym, Dict{Symbol,Any}())
                mt[:default] = def
            end
        end
    end

    # obsf & obssym
    if haskey(dict, :obsf) || haskey(dict, :obssym)
        if !(haskey(dict, :obsf) && haskey(dict, :obssym))
            throw(ArgumentError("If `obsf` is provided, `obssym` must be provided as well."))
        end
        if _has_metadata(dict[:obssym])
            dict[:obssym], _metadata = _split_metadata(dict[:obssym])
            mergewith!(merge!, symmetadata, _metadata)
        end
    else
        dict[:obsf] = nothing
        dict[:obssym] = Symbol[]
    end

    # inputsym
    if T <: VertexFunction
        if haskey(dict, :inputsym_src) || haskey(dict, :inputsym_dst)
            throw(ArgumentError("Keywords `inputsym_src` and `inputsym_dst` are not valid for $T."))
        end
        if haskey(dict, :inputsym) && !isnothing(dict[:inputsym])
            if _has_metadata(dict[:inputsym])
                dict[:inputsym], _metadata = _split_metadata(dict[:inputsym])
                mergewith!(merge!, symmetadata, _metadata)
            end
        else
            dict[:inputsym] = nothing
        end
    elseif T <: EdgeFunction
        haskey(dict, :inputsym) && throw(ArgumentError("Keywords `inputsym` is not valid for $T."))
        if haskey(dict, :inputsym_src) || haskey(dict, :inputsym_dst)
            if !(haskey(dict, :inputsym_src) && haskey(dict, :inputsym_dst))
                throw(ArgumentError("Both `inputsym_src` and `inputsym_dst` must be provided."))
            end
            src, _metadata = _split_metadata(pop!(dict, :inputsym_src))
            mergewith!(merge!, symmetadata, _metadata)
            dst, _metadata = _split_metadata(pop!(dict, :inputsym_dst))
            mergewith!(merge!, symmetadata, _metadata)
            dict[:inputsym] = (; src, dst)
        else
            dict[:inputsym] = nothing
        end
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
                if length(mm) == dim
                    dict[:mass_matrix] = LinearAlgebra.Diagonal(mm)
                else
                    throw(ArgumentError("If given as a vector, mass matrix must have length equal to dimension of component."))
                end
            elseif mm isa Number # convert to uniform scaling
                dict[:mass_matrix] = LinearAlgebra.UniformScaling(mm)
            elseif mm isa AbstractMatrix
                @argcheck size(mm) == (dim, dim) "Size of mass matrix must match dimension of component."
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
            dict[:depth] = dim
        elseif T<:EdgeFunction
            coupling = dict[:coupling]
            dict[:depth] = coupling==Fiducial() ? floor(Int, dim/2) : dim
        else
            throw(ArgumentError("Cannot construct $T: default depth not known."))
        end
    end
    if haskey(dict, :coupling) && dict[:coupling]==Fiducial() && dict[:depth] > floor(dim/2)
        throw(ArgumentError("Depth cannot exceed half the dimension for Fiducial coupling."))
    elseif dict[:depth] > dim
        throw(ArgumentError("Depth cannot exceed half the dimension."))
    end

    # check for name clashes (at the end because only now sym, psym, obssym are initialized)
    _s  = dict[:sym]
    _ps = dict[:psym]
    _os = dict[:obssym]
    __is = dict[:inputsym]
    _is = if isnothing(__is)
        Symbol[]
    elseif __is isa NamedTuple
        vcat(__is.src, __is.dst)
    else
        __is
    end
    if !allunique(vcat(_s, _ps, _os, _is))
        throw(ArgumentError("Symbol names must be unique. There are clashes in sym, psym, obssym and inputsym."))
    end

    return dict
end

_default_name(::Type{StaticVertex}) = :StaticVertex
_default_name(::Type{ODEVertex}) = :ODEVertex
_default_name(::Type{StaticEdge}) = :StaticEdge
_default_name(::Type{ODEEdge}) = :ODEEdge

_has_metadata(vec::AbstractVector{<:Symbol}) = false
_has_metadata(vec::AbstractVector{<:Pair}) = true
_has_metadata(vec::AbstractVector) = any(el -> el isa Pair, vec)
function _split_metadata(input)
    Base.require_one_based_indexing(input)
    syms = Vector{Symbol}(undef, length(input))
    metadata = Dict{Symbol,Dict{Symbol,Any}}()
    for i in eachindex(input)
        if input[i] isa Pair
            sym = input[i].first
            dat  = input[i].second
            syms[i] = sym
            metadata[sym] = if dat isa Number
                Dict(:default => dat)
            elseif dat isa NamedTuple
                Dict(zip(keys(dat), values(dat)))
            else
                dat
            end
        else
            syms[i] = input[i]
        end
    end
    syms, metadata
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

_valid_signature(::Type{<:StaticVertex}, f) = _takes_n_vectors(f, 3) #(u, edges, p, t)
_valid_signature(::Type{<:ODEVertex}, f) = _takes_n_vectors(f, 4) #(du, u, edges, p, t)
_valid_signature(::Type{<:StaticEdge}, f) = _takes_n_vectors(f, 4) #(u, src, dst, p, t)
_valid_signature(::Type{<:ODEEdge}, f) = _takes_n_vectors(f, 5) #(du, u, src, dst, p, t)

_takes_n_vectors(f, n) = hasmethod(f, (Tuple(Vector{Float64} for i in 1:n)..., Float64))

"""
    copy(c::NetworkDynamics.ComponentFunction)

Shallow copy of the component function. Creates a deepcopy of `metadata` and `symmetadata`
but references the same objects everywhere else.
"""
@generated function Base.copy(c::ComponentFunction)
    fields = fieldnames(c)
    # fields to copy
    cfields = (:metadata, :symmetadata)
    # normal fields
    nfields = setdiff(fields, cfields)
    assign = Expr(:block,
        (:($(field) = c.$field) for field in nfields)...,
        (:($(field) = deepcopy(c.$field)) for field in cfields)...)
    construct = Expr(:call, c, [:($field) for field in fields]...)

    quote
       $assign
       $construct
    end
end

Base.hash(cf::ComponentFunction, h::UInt) = hash_fields(cf, h)
function Base.:(==)(cf1::ComponentFunction, cf2::ComponentFunction)
    typeof(cf1) == typeof(cf2) && equal_fields(cf1, cf2)
end
