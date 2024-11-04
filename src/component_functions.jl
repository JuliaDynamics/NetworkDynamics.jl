abstract type FeedForwardType end
struct PureFeedForward <: FeedForwardType end
struct FeedForward <: FeedForwardType end
struct NoFeedForward <: FeedForwardType end
struct PureStateMap <: FeedForwardType end
hasfftype(::Any) = false
hasff(x) = fftype(x) isa Union{PureFeedForward, FeedForward}

_hassymbolmapping(::Any) = false


struct StateMask{N,I}
    idxs::I
    StateMask(i) = new{length(i), typeof(i)}(i)
end
hasfftype(::StateMask) = true
fftype(::StateMask) = PureStateMap()

@inline function (s::StateMask)(out, u)
    @inbounds for i in eachindex(s.idxs)
        out[i] = u[s.idxs[i]]
    end
    nothing
end


abstract type Coupling{FF} end
hasfftype(::Coupling) = true
fftype(::Coupling{FF}) where {FF} = FF()

struct AntiSymmetric{FF,G} <: Coupling{FF}
    g::G
    AntiSymmetric(g; ff=_infer_ss_fftype(g)) = new{typeof(ff), typeof(g)}(g)
end
@inline function (c::AntiSymmetric)(osrc, odst, args...)
    @inline c.g(odst, args...)
    @inbounds for i in 1:length(osrc)
        osrc[i] = -odst[i]
    end
    nothing
end

struct Symmetric{FF,G} <: Coupling{FF}
    g::G
    Symmetric(g; ff=_infer_ss_fftype(g)) = new{typeof(ff), typeof(g)}(g)
end
@inline function (c::Symmetric)(osrc, odst, args...)
    @inline c.g(odst, args...)
    @inbounds for i in 1:length(osrc)
        osrc[i] = odst[i]
    end
    nothing
end

struct Directed{FF} <: Coupling{FF} end
# TODO directed

struct Fiducial{FF,GS,GD} <: Coupling{FF}
    src::GS
    dst::GD
    function Fiducial(src, dst; ff=nothing)
        if isnothing(ff)
            ffsrc = _infer_ss_fftype(src)
            ffdst = _infer_ss_fftype(dst)
            if ffsrc == ffdst
                ff = ffsrc
            else
                error("Both src and dst coupling functions have different ff types.")
            end
        end
        new{typeof(ff), typeof(src), typeof(dst)}(src, dst)
    end
end
Fiducial(src::UnitRange, dst::UnitRange; ff=nothing) = Fiducial(StateMask(src), StateMask(dst); ff)
Fiducial(;src, dst, ff=nothing) = Fiducial(src, dst; ff)

@inline function (c::Fiducial)(osrc, odst, args...)
    @inline c.src(osrc, args...)
    @inline c.dst(odst, args...)
    nothing
end


abstract type ComponentFunction end

"""
Abstract supertype for all vertex functions.
"""
abstract type VertexFunction <: ComponentFunction end

"""
Abstract supertype for all edge functions.
"""
# abstract type EdgeFunction{C<:Coupling} <: ComponentFunction end
abstract type EdgeFunction <: ComponentFunction end

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

outsym(c::ComponentFunction) = c.outsym
outdim(c::VertexFunction) = length(outsym(c))
outdim(c::EdgeFunction) = (; dst=outdim_dst(c), src=outdim_src(c))
outdim_src(c::EdgeFunction) = length(outsym(c).src)
outdim_dst(c::EdgeFunction) = length(outsym(c).dst)

struct UnifiedVertex{F,G,FFT,OF,MM} <: VertexFunction
    name::Symbol
    # main function
    f::F
    sym::Vector{Symbol}
    mass_matrix::MM
    # outputs
    g::G
    outsym::Vector{Symbol}
    ff::FFT
    # parameters and option input sym
    psym::Vector{Symbol}
    inputsym::Union{Nothing, Vector{Symbol}}
    # observed
    obsf::OF
    obssym::Vector{Symbol}
    # optional inputsyms
    symmetadata::Dict{Symbol,Dict{Symbol, Any}}
    metadata::Dict{Symbol,Any}
end
UnifiedVertex(; kwargs...) = _construct_comp(UnifiedVertex, kwargs)

struct UnifiedEdge{F,G,FFT,OF,MM} <: EdgeFunction
    name::Symbol
    # main function
    f::F
    sym::Vector{Symbol}
    mass_matrix::MM
    # outputs
    g::G
    outsym::@NamedTuple{dst::Vector{Symbol},src::Vector{Symbol}}
    ff::FFT
    # parameters and option input sym
    psym::Vector{Symbol}
    inputsym::Union{Nothing, @NamedTuple{src::Vector{Symbol},dst::Vector{Symbol}}}
    # observed
    obsf::OF
    obssym::Vector{Symbol}
    # metadata
    symmetadata::Dict{Symbol,Dict{Symbol, Any}}
    metadata::Dict{Symbol,Any}
end
UnifiedEdge(; kwargs...) = _construct_comp(UnifiedEdge, kwargs)

fftype(c::Union{<:UnifiedVertex, <:UnifiedEdge}) = c.ff

function _infer_fftype(::Type{<:VertexFunction}, g, dim)
    pureff = _takes_n_vectors_and_t(g, 3) && iszero(dim)  # (out, ein, p, t)
    ff     = _takes_n_vectors_and_t(g, 4)                 # (out, u, ein, p, t)
    noff   = _takes_n_vectors_and_t(g, 3) && !iszero(dim) # (out, u, p, t)
    pureu  = _takes_n_vectors_no_t(g, 2)                  # (out, states)

    ff+noff+pureff+pureu > 1 && error("Could not determinen output map type from g signature. Provide :ff keyword explicitly!")
    ff+noff+pureff+pureu == 0 && error("Method signature of `g` seems invalid!")

    pureff && return PureStateMap()
    ff && return FeedForward()
    noff && return NoFeedForward()
    pureu && return PureStateMap()
end

_infer_ss_fftype(g) = _infer_fftype(EdgeFunction, g, -1; singlesided=true)
function _infer_fftype(::Type{<:EdgeFunction}, g, dim; singlesided=false)
    # if singlesided we allways have one less
    pureff = _takes_n_vectors_and_t(g, 5-singlesided) # ([osrc], odst, src, dst, p, t)
    ff     = _takes_n_vectors_and_t(g, 6-singlesided) # ([osrc], odst, u, src, dst, p, t)
    noff   = _takes_n_vectors_and_t(g, 4-singlesided) # ([osrc], odst, u, p, t)
    pureu  = _takes_n_vectors_no_t(g, 3-singlesided)  # ([osrc], odst, u)

    ff+noff+pureff+pureu > 1 && error("Could not determinen output map type from g signature. Provide :ff keyword explicitly!")
    ff+noff+pureff+pureu == 0 && error("Method signature of `g` seems invalid!")

    pureff && return PureFeedForward()
    ff && return FeedForward()
    noff && return NoFeedForward()
    pureu && return PureStateMap()
end

compg(c::ComponentFunction) = c.g

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

struct StaticEdge{C,F,OF} <: EdgeFunction
    @CommonFields
    coupling::C
end
StaticEdge(; kwargs...) = _construct_comp(StaticEdge, kwargs)
StaticEdge(f; kwargs...) = StaticEdge(;f, kwargs...)
StaticEdge(f, dim, coupling; kwargs...) = StaticEdge(;f, _dimsym(dim)..., coupling, kwargs...)
StaticEdge(f, dim, pdim, coupling; kwargs...) = StaticEdge(;f, _dimsym(dim, pdim)..., coupling, kwargs...)
StaticEdge(e::StaticEdge; kwargs...) = _reconstruct_comp(StaticEdge, e, kwargs)

struct ODEEdge{C,F,OF,MM} <: EdgeFunction
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
statetype(::Type{<:UnifiedEdge}) = Dynamic()
statetype(::Type{<:UnifiedVertex}) = Dynamic()

isdynamic(x) = statetype(x) == Dynamic()
isstatic(x)  = statetype(x) == Static()

"""
    dispatchT(<:ComponentFunction) :: Type{<:ComponentFunction}

Returns the type "essence" of the component used for dispatch.
Fills up type parameters with `nothing` to ensure `Core.Compiler.isconstType`
for GPU compatibility.
"""
dispatchT(::T) where {T<:ComponentFunction} = dispatchT(T)
dispatchT(::Type{<:StaticVertex}) = StaticVertex{nothing,nothing}
dispatchT(::Type{<:ODEVertex}) = ODEVertex{nothing,nothing,nothing}
dispatchT(T::Type{<:StaticEdge}) = StaticEdge{typeof(coupling(T)),nothing,nothing}
dispatchT(T::Type{<:ODEEdge}) = ODEEdge{typeof(coupling(T)),nothing,nothing,nothing}
dispatchT(T::Type{<:UnifiedVertex}) = UnifiedVertex{nothing,nothing,nothing,nothing,nothing}
dispatchT(T::Type{<:UnifiedEdge}) = UnifiedEdge{nothing,nothing,nothing,nothing,nothing}

batchequal(a, b) = false
function batchequal(a::ComponentFunction, b::ComponentFunction)
    for f in (compf, compg, dim, pdim, odim)
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
    _maybewrap!(dict, :outsym, Symbol)
    _maybewrap!(dict, :obssym, Symbol)

    if !haskey(dict, :g)
        throw(ArgumentError("Output function g musst be provided!"))
    end

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

    # set f to nothing if not present
    if !haskey(dict, :f)
        dict[:f] = nothing
        if haskey(dict, :dim) && dict[:dim] != 0
            throw(ArgumentError("Cannot provide dim != 0 without f."))
        end
        if haskey(dict, :sym) && !isempty(dict[:sym])
            throw(ArgumentError("Cannot provide sym without f."))
        end
        dict[:sym] = Symbol[]
        dict[:dim] = 0
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

    # infer fftype
    if !haskey(dict, :ff)
        g = dict[:g]
        dict[:ff] = hasfftype(g) ? fftype(g) : _infer_fftype(T, g, dim)
    end

    # outsym & outdim
    # _hassymbolmapping(dict[:g]) || haskey(dict, :outdim) || haskey(dict, :outsym) || throw(ArgumentError("Either `outdim` or `outsym` must be provided to construct $T with arbitray g."))
    if haskey(dict, :outsym)
        if haskey(dict, :outdim)
            if dict[:outdim] != length(dict[:outsym])
                throw(ArgumentError("Length of outsym and outdim must match."))
            end
            # @warn "Unnecessary kw outdim, can be infered from outsym."
            delete!(dict, :outdim)
        end
        outsym = dict[:outsym]

        # extract and merge metadata from outsym
        if T <: VertexFunction
            if has_metadata(outsym)
                dict[:outsym], _metadata = _split_metadata(outsym)
                mergewith!(merge!, symmetadata, _metadata)
            end
        elseif T <: EdgeFunction
            (; src, dst) = outsym
            src = if has_metadata(src)
                newsrc, _metadata = _split_metadata(src)
                mergewith!(merge!, symmetadata, _metadata)
                newsrc
            end
            dst = if has_metadata(dst)
                newdst, _metadata = _split_metadata(dst)
                mergewith!(merge!, symmetadata, _metadata)
                newdst
            end
            dict[:outsym] = (; dst, src)
        end
    elseif _hassymbolmapping(dict[:g])
        dict[:outsym] = _mapsymbols(dict[:g], dict[:sym])
        if haskey(dict, :outdim)
            delete!(dict, :outdim)
            lnegth
        end
    elseif haskey(dict, :outdim)
        _dim = pop!(dict, :outdim)
        if T <: VertexFunction
            dict[:outsym] = [_dim>1 ? Symbol("vout", subscript(i)) : :sout for i in 1:_dim]
        elseif T <: EdgeFunction
            base_syms = [_dim>1 ? Symbol("out", subscript(i)) : :out for i in 1:_dim]
            dst = map(x -> Symbol("dst₊", x), base_syms)
            src = map(x -> Symbol("src₊", x), base_syms)
            dict[:outsym] = (; dst, src)
        else
            error()
        end
    else
        throw(ArgumentError("Either `outdim` or `outsym` must be provided to construct $T \
            with arbitray g. Outsyms can be inferred for `StateMask` output functions."))
    end

    outdim = length(dict[:outsym])

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

# define the symbolmapping to infer output symbols from state symbols
_hassymbolmapping(::StateMask) = true
_hassymbolmapping(::AntiSymmetric{<:Any, <:StateMask}) = true
_hassymbolmapping(::Symmetric{<:Any, <:StateMask}) = true
_hassymbolmapping(::Fiducial{<:Any, <:StateMask, <:StateMask}) = true

_mapsymbols(g::StateMask, s) = s[g.idxs]
function _mapsymbols(g::AntiSymmetric{<:Any, <:StateMask}, s)
    dst = _mapsymbols(g.g, s)
    src = map(x->Symbol("₋", x), dst)
    (;dst, src)
end
function _mapsymbols(g::Symmetric{<:Any, <:StateMask}, s)
    dst = _mapsymbols(g.g, s)
    src = dst
    (;dst, src)
end
function _mapsymbols(g::Fiducial{<:Any, <:StateMask, <:StateMask}, s)
    dst = _mapsymbols(g.dst, s)
    src = _mapsymbols(g.src, s)
    (;dst, src)
end


_default_name(::Type{StaticVertex}) = :StaticVertex
_default_name(::Type{ODEVertex}) = :ODEVertex
_default_name(::Type{UnifiedVertex}) = :UnifiedVertex
_default_name(::Type{StaticEdge}) = :StaticEdge
_default_name(::Type{ODEEdge}) = :ODEEdge
_default_name(::Type{UnifiedEdge}) = :UnifiedEdge

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

_valid_signature(::Type{<:StaticVertex}, f) = _takes_n_vectors_and_t(f, 3) #(u, edges, p, t)
_valid_signature(::Type{<:ODEVertex}, f) = _takes_n_vectors_and_t(f, 4) #(du, u, edges, p, t)
_valid_signature(::Type{<:StaticEdge}, f) = _takes_n_vectors_and_t(f, 4) #(u, src, dst, p, t)
_valid_signature(::Type{<:ODEEdge}, f) = _takes_n_vectors_and_t(f, 5) #(du, u, src, dst, p, t)

_takes_n_vectors_and_t(f, n) = hasmethod(f, (Tuple(Vector{Float64} for i in 1:n)..., Float64))
_takes_n_vectors_no_t(f, n) = hasmethod(f, Tuple(Vector{Float64} for i in 1:n))

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
