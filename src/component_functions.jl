abstract type FeedForwardType end
struct PureFeedForward <: FeedForwardType end
struct FeedForward <: FeedForwardType end
struct NoFeedForward <: FeedForwardType end
struct PureStateMap <: FeedForwardType end
hasfftype(::Any) = false
hasff(x) = fftype(x) isa Union{PureFeedForward, FeedForward}

_has_sym_to_outsym_mapping(::Any) = false

struct StateMask{N,I}
    idxs::I
    function StateMask(i::AbstractArray)
        if !isbitstype(typeof(i))
           i = SVector{length(i)}(i)
        end
        new{length(i), typeof(i)}(i)
    end
end
StateMask(i::Number) = StateMask(SVector{1}(i))
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
    function AntiSymmetric(g; ff=nothing)
        ff = isnothing(ff) ? _infer_ss_fftype(g) : ff
        new{typeof(ff), typeof(g)}(g)
    end
end
AntiSymmetric(g::Union{AbstractVector,Number}; ff=nothing) = AntiSymmetric(StateMask(g); ff)
@inline function (c::AntiSymmetric)(osrc, odst, args...)
    @inline c.g(odst, args...)
    @inbounds for i in 1:length(osrc)
        osrc[i] = -odst[i]
    end
    nothing
end

struct Symmetric{FF,G} <: Coupling{FF}
    g::G
    function Symmetric(g; ff=nothing)
        ff = isnothing(ff) ? _infer_ss_fftype(g) : ff
        new{typeof(ff), typeof(g)}(g)
    end
end
Symmetric(g::Union{AbstractVector,Number}; ff=nothing) = Symmetric(StateMask(g); ff)
@inline function (c::Symmetric)(osrc, odst, args...)
    @inline c.g(odst, args...)
    @inbounds for i in 1:length(osrc)
        osrc[i] = odst[i]
    end
    nothing
end

struct Directed{FF,G} <: Coupling{FF}
    g::G
    function Directed(g; ff=nothing)
        ff = isnothing(ff) ? _infer_ss_fftype(g) : ff
        new{typeof(ff), typeof(g)}(g)
    end
end
Directed(g::Union{AbstractVector,Number}; ff=nothing) = Directed(StateMask(g); ff)
@inline function (c::Directed)(osrc, odst, args...)
    @inline c.g(odst, args...)
    nothing
end

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
function Fiducial(src::Union{AbstractVector,Number}, dst::Union{AbstractVector,Number}; ff=nothing)
    Fiducial(StateMask(src), StateMask(dst); ff)
end
Fiducial(;src, dst, ff=nothing) = Fiducial(src, dst; ff)

@inline function (c::Fiducial)(osrc, odst, args...)
    @inline c.src(osrc, args...)
    @inline c.dst(odst, args...)
    nothing
end


abstract type ComponentFunction end

struct VertexFunction{F,G,FFT,OF,MM} <: ComponentFunction
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
    insym::Union{Nothing, Vector{Symbol}}
    # observed
    obsf::OF
    obssym::Vector{Symbol}
    # optional insyms
    symmetadata::Dict{Symbol,Dict{Symbol, Any}}
    metadata::Dict{Symbol,Any}
    # cached symbol collections
    _outsym_flat::Vector{Symbol} # outsyms as they appear in outbuf
    _obssym_all::Vector{Symbol}  # collection of true "observed" (flat_out\statesym) ∪ obssym
end
VertexFunction(; kwargs...) = _construct_comp(VertexFunction, Base.inferencebarrier(kwargs))
VertexFunction(v::VertexFunction; kwargs...) = _reconstruct_comp(VertexFunction, v, Base.inferencebarrier(kwargs))

struct EdgeFunction{F,G,FFT,OF,MM} <: ComponentFunction
    name::Symbol
    # main function
    f::F
    sym::Vector{Symbol}
    mass_matrix::MM
    # outputs
    g::G
    outsym::@NamedTuple{src::Vector{Symbol},dst::Vector{Symbol}}
    ff::FFT
    # parameters and option input sym
    psym::Vector{Symbol}
    insym::Union{Nothing, @NamedTuple{src::Vector{Symbol},dst::Vector{Symbol}}}
    # observed
    obsf::OF
    obssym::Vector{Symbol}
    # metadata
    symmetadata::Dict{Symbol,Dict{Symbol, Any}}
    metadata::Dict{Symbol,Any}
    # cached symbol collections
    _outsym_flat::Vector{Symbol} # outsyms as they appear in outbuf
    _obssym_all::Vector{Symbol}  # collection of true "observed" (flat_out\statesym) ∪ obssym
end
EdgeFunction(; kwargs...) = _construct_comp(EdgeFunction, Base.inferencebarrier(kwargs))
EdgeFunction(v::EdgeFunction; kwargs...) = _reconstruct_comp(EdgeFunction, v, Base.inferencebarrier(kwargs))

"""
    compf(c::ComponentFunction)

Retrieve internal function `f` of the component.
"""
compf(c::ComponentFunction) = c.f

"""
    compg(c::ComponentFunction)

Retrieve output function `g` of the component.
"""
compg(c::ComponentFunction) = c.g

"""
    fftype(c::ComponentFunction)

Retrieve the feed forward type of the component.
"""
fftype(c::ComponentFunction) = c.ff

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
    outdim(c::EdgeFunction)::@NamedTuple(src::Int, dst::Int)
    outdim(c::VertexFunction)::Int

Retrieve the output dimension of the component
"""
outdim(c::VertexFunction) = length(outsym(c))
outdim(c::EdgeFunction) = (;src=outdim_src(c), dst=outdim_dst(c))

outdim_src(c::EdgeFunction) = length(outsym(c).src)
outdim_dst(c::EdgeFunction) = length(outsym(c).dst)

"""
   outsym(c::VertexFunction)::Vector{Symbol}
   outsym(c::EdgeFunction)::@NamedTuple{src::Vector{Symbol}, dst::Vector{Symbol}}

Retrieve the output symbols of the component.
"""
outsym(c::ComponentFunction) = c.outsym

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
    symmetadata(c::ComponentFunction)::Dict{Symbol,Dict{Symbol,Any}}

Retrieve the metadata dictionary for the symbols. Keys are the names of the
symbols as they appear in [`sym`](@ref), [`psym`](@ref), [`obssym`](@ref) and [`insym`](@ref).

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
    hasinsym(c::ComponentFunction)

Checks if the optioan field `insym` is present in the component function.
"""
hasinsym(c::ComponentFunction) = !isnothing(c.insym)
hasindim(c::ComponentFunction) = hasinsym(c)

"""
    insym(c::VertexFunction)::Vector{Symbol}
    insym(c::EdgeFunction)::NamedTuple with :src and :dst keys

Musst be called *after* [`hasinsym`](@ref)/[`hasindim`](@ref) returned true.
Gives the `insym` vector(s). For vertex functions just a single vector, for
edges it returns a named tuple `(; src, dst)` with two symbol vectors.
"""
insym(c::VertexFunction)::Vector{Symbol} = c.insym
insym(c::EdgeFunction)::@NamedTuple{src::Vector{Symbol},dst::Vector{Symbol}} = c.insym

"""
    indim(c::VertexFunction)::Int
    indim(c::EdgeFunction)::@NamedTuple{src::Int,dst::Int}

Musst be called *after* [`hasinsym`](@ref)/[`hasindim`](@ref) returned true.
Gives the input dimension(s).
"""
indim(c::VertexFunction)::Int = length(insym(c))
indim(c::EdgeFunction)::@NamedTuple{src::Int,dst::Int} = (; src=length(insym(c).src), dst=length(insym(c).dst))

# return both "observed" outputs (those that do not shadow states) and true observed
outsym_flat(c::ComponentFunction) = c._outsym_flat
obssym_all(c::ComponentFunction) = c._obssym_all

_infer_ss_fftype(g) = _infer_fftype(g, 1, 2, nothing)
_infer_fftype(::Type{<:VertexFunction}, g, dim) = _infer_fftype(g, 1, 1, dim)
_infer_fftype(::Type{<:EdgeFunction}, g, dim) = _infer_fftype(g, 2, 2, dim)

function _infer_fftype(g, nout, nin, dim)
    pureff = _takes_n_vecs_and_t(g, nout + nin + 1)     # (outs..., ins..., p, t)
    ff     = _takes_n_vecs_and_t(g, nout + 1 + nin + 1) # (outs..., u, ins..., p, t)
    noff   = _takes_n_vecs_and_t(g, nout + 1 + 1)       # (outs..., u, p, t)
    pureu  = _takes_n_vecs_no_t(g, nout + 1)            # (outs..., u)

    # if nin = 1 we cannot distiguish between pureff and noff based on arguments,
    # resolve using dimension
    if nin==1 && pureff && noff
        pureff = iszero(dim)
        noff = !iszero(dim)
    end

    ff+noff+pureff+pureu > 1 && error("Could not determinen output map type from g signature. Provide :ff keyword explicitly!")
    ff+noff+pureff+pureu == 0 && error("Method signature of `g` seems invalid!")

    pureff && return PureFeedForward()
    ff && return FeedForward()
    noff && return NoFeedForward()
    pureu && return PureStateMap()
end

_takes_n_vecs_and_t(f, n) = hasmethod(f, (Tuple(Vector{Float64} for i in 1:n)..., Float64))
_takes_n_vecs_no_t(f, n) = hasmethod(f, Tuple(Vector{Float64} for i in 1:n))

"""
    dispatchT(<:ComponentFunction) :: Type{<:ComponentFunction}

Returns the type "essence" of the component used for dispatch.
Fills up type parameters with `nothing` to ensure `Core.Compiler.isconstType`
for GPU compatibility.
"""
dispatchT(::T) where {T<:ComponentFunction} = dispatchT(T)
dispatchT(T::Type{<:VertexFunction}) = VertexFunction{nothing,nothing,nothing,nothing,nothing}
dispatchT(T::Type{<:EdgeFunction}) = EdgeFunction{nothing,nothing,nothing,nothing,nothing}

# TODO: introduce batchequal hash for faster batching of component functions
batchequal(a, b) = false
function batchequal(a::ComponentFunction, b::ComponentFunction)
    compf(a) == compf(b) || return false
    compg(a) == compg(b) || return false
    dim(a)   == dim(b)   || return false
    pdim(a)  == pdim(b)  || return false
    outdim(a)  == outdim(b)  || return false
    return true
end

"""
    _construct_comp(::Type{T}, kwargs) where {T}

Internal function to construct a component function from keyword arguments.
Fills up kw arguments with default values and performs sanity checks.
"""
function _construct_comp(::Type{T}, @nospecialize(kwargs)) where {T}
    dict = _fill_defaults(T, Base.inferencebarrier(kwargs))

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
function _fill_defaults(T, @nospecialize(kwargs))
    dict = Dict{Symbol, Any}(kwargs)

    # syms might be provided as single pairs or symbols, wrap in vector
    _maybewrap!(dict, :sym, Union{Symbol, Pair})
    _maybewrap!(dict, :psym, Union{Symbol, Pair})
    _maybewrap!(dict, :insym, Union{Symbol,Pair}; allownt=T <: EdgeFunction)
    _maybewrap!(dict, :outsym, Union{Symbol,Pair}; allownt=T <: EdgeFunction)
    _maybewrap!(dict, :obssym, Symbol)

    symmetadata = get!(dict, :symmetadata, Dict{Symbol,Dict{Symbol,Any}}())

    metadata = try
        convert(Dict{Symbol,Any}, get!(dict, :metadata, Dict{Symbol,Any}()))
    catch e
        throw(ArgumentError("Provided metadata keyword musst be a Dict{Symbol,Any}. Got $(repr(dict[:metadata]))."))
    end

    ####
    #### graphelement
    ####
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

    ####
    #### f and g
    ####
    if !haskey(dict, :g)
        throw(ArgumentError("Output function g musst be provided!"))
    end
    g = dict[:g]
    if g isa Union{AbstractVector, Number}
        g = dict[:g] = StateMask(g)
    end
    if T <: EdgeFunction && g isa StateMask
        throw(ArgumentError("StateMask cannot be used directly for EdgeFunction. Wrap in Fiducial, Symmetric, AntiSymmetric or Directed instead."))
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
    end


    ####
    #### variable sym
    ####
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
    sym = dict[:sym]

    # infer fftype (needs dim)
    if !haskey(dict, :ff)
        dict[:ff] = hasfftype(g) ? fftype(g) : _infer_fftype(T, g, dim)
    end

    ####
    #### parameter sym
    ####
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
    psym = dict[:psym]

    ####
    #### output sym
    ####
    if !haskey(dict, :outsym)
        if _has_sym_to_outsym_mapping(g)
            dict[:outsym] = _sym_to_outsym(g, sym)
        elseif haskey(dict, :outdim)
            _outdim = pop!(dict, :outdim)
            dict[:outsym] = [(_outdim>1 ? Symbol("o", subscript(i)) : :o) for i in 1:_outdim]
        else
            throw(ArgumentError("Either `outdim` or `outsym` must be provided to construct $T \
                with arbitray g. Outsyms can be inferred for `StateMask` output functions."))
        end
    end
    outsym = dict[:outsym]
    if T <: EdgeFunction && outsym isa AbstractVector
        if _has_metadata(outsym)
            throw(ArgumentError("Metadata for outsyms can only be provided when using full (;src, dst) form"))
        end
        outsym = dict[:outsym] = _symvec_to_sym_tup(g, outsym)
    end

    if haskey(dict, :outdim)
        _outdim = pop!(dict, :outdim)
        if T <: EdgeFunction
            if _outdim isa NamedTuple && (_outdim.src != length(outsym.src) || _outdim.dst != length(outsym.dst))
                throw(ArgumentError("Length of outsym and outdim must match."))
            elseif _outdim != length(outsym.dst) || !(g isa Directed) && _outdime != length(outsym.src)
                throw(ArgumentError("Length of outsym and outdim must match."))
            end
        elseif T <: VertexFunction
            if _outdim != length(outsym)
                throw(ArgumentError("Length of outsym and outdim must match."))
            end
        end
    end

    # extract and merge metadata from outsym
    if T <: VertexFunction
        if _has_metadata(outsym)
            dict[:outsym], _metadata = _split_metadata(outsym)
            mergewith!(merge!, symmetadata, _metadata)
        end
    elseif T <: EdgeFunction
        (; src, dst) = outsym
        if _has_metadata(src)
            src, _metadata = _split_metadata(src)
            mergewith!(merge!, symmetadata, _metadata)
        end
        if _has_metadata(dst)
            dst, _metadata = _split_metadata(dst)
            mergewith!(merge!, symmetadata, _metadata)
        end
        dict[:outsym] = (; src, dst)
    end

    ####
    #### insym
    ####
    if !haskey(dict, :insym)
        if haskey(dict, :indim)
            _indim = pop!(dict, :indim)
            dict[:insym] = [(_indim>1 ? Symbol("i", subscript(i)) : :i) for i in 1:_indim]
        else
            dict[:insym] = nothing
        end
    end
    insym = dict[:insym]
    if T <: EdgeFunction && insym isa AbstractVector
        if _has_metadata(insym)
            throw(ArgumentError("Metadata for insyms can only be provided when using full (;src, dst) form"))
        end
        insym = dict[:insym] = _symvec_to_sym_tup(g, insym)
    end
    if haskey(dict, :indim)
        _indim = pop!(dict, :indim)
        if T <: EdgeFunction
            if _indim isa NamedTuple && (_indim.src != length(insym.src) || _indim.dst != length(insym.dst))
                throw(ArgumentError("Length of insym and indim must match."))
            elseif _indim != length(insym.src) || _indime != length(insym.dst)
                throw(ArgumentError("Length of insym and indim must match."))
            end
        elseif T <: VertexFunction
            if _indim != length(outsym)
                throw(ArgumentError("Length of insym and indim must match."))
            end
        end
    end
    # extract and merge metadata from insym
    if !isnothing(insym)
        if T <: VertexFunction
            if _has_metadata(insym)
                dict[:insym], _metadata = _split_metadata(insym)
                mergewith!(merge!, symmetadata, _metadata)
            end
        elseif T <: EdgeFunction
            (; src, dst) = insym
            if _has_metadata(src)
                src, _metadata = _split_metadata(src)
                mergewith!(merge!, symmetadata, _metadata)
            end
            if _has_metadata(dst)
                dst, _metadata = _split_metadata(dst)
                mergewith!(merge!, symmetadata, _metadata)
            end
            dict[:insym] = (; src, dst)
        end
    end

    ####
    #### obsf & obssym
    ####
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
    obssym = dict[:obssym]

    ####
    #### name
    ####
    if !haskey(dict, :name)
        dict[:name] = _default_name(T)
    end

    ####
    #### mass matrix
    ####
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

    ####
    #### Cached outsymflat/outsymall
    ####
    _outsym_flat = if T <: VertexFunction
        outsym
    elseif T <: EdgeFunction
        vcat(outsym.src, outsym.dst)
    else
        error()
    end
    dict[:_outsym_flat] = _outsym_flat
    dict[:_obssym_all] = setdiff(_outsym_flat, sym) ∪ obssym


    # check for name clashes (at the end because only now sym, psym, obssym are initialized)
    _s  = sym
    _ps = psym
    _obss = obssym
    _os = outsym
    _os = if _has_sym_to_outsym_mapping(g)
       setdiff(outsym, sym) # allow name clashes for statemask derivatives
    else
        outsym
    end
    __is = insym

    _is = if isnothing(__is)
        Symbol[]
    elseif __is isa NamedTuple
        vcat(__is.src, __is.dst)
    else
        __is
    end
    if !allunique(vcat(_s, _ps, _obss, _is, _os))
        throw(ArgumentError("Symbol names must be unique. There are clashes in sym, psym, outsym, obssym and insym."))
    end

    return dict
end

# define the symbolmapping to infer output symbols from state symbols
_has_sym_to_outsym_mapping(::StateMask) = true
_has_sym_to_outsym_mapping(::Directed{<:Any, <:StateMask}) = true
_has_sym_to_outsym_mapping(::AntiSymmetric{<:Any, <:StateMask}) = true
_has_sym_to_outsym_mapping(::Symmetric{<:Any, <:StateMask}) = true
_has_sym_to_outsym_mapping(::Fiducial{<:Any, <:StateMask, <:StateMask}) = true

_sym_to_outsym(g::StateMask, s::AbstractVector{Symbol}) = s[g.idxs]
function _sym_to_outsym(g::AntiSymmetric{<:Any, <:StateMask}, s::AbstractVector{Symbol})
    s = _sym_to_outsym(g.g, s)
    _symvec_to_sym_tup(g, s)
end
function _sym_to_outsym(g::Symmetric{<:Any, <:StateMask}, s::AbstractVector{Symbol})
    s = _sym_to_outsym(g.g, s)
    _symvec_to_sym_tup(g, s)
end
function _sym_to_outsym(g::Directed{<:Any, <:StateMask}, s::AbstractVector{Symbol})
    s = _sym_to_outsym(g.g, s)
    _symvec_to_sym_tup(g, s)
end
function _sym_to_outsym(g::Fiducial{<:Any, <:StateMask, <:StateMask}, s::AbstractVector{Symbol})
    dst = _sym_to_outsym(g.dst, s)
    src = _sym_to_outsym(g.src, s)
    (; src, dst)
end
function _symvec_to_sym_tup(g, s::AbstractVector{Symbol})
    src = map(x->Symbol("src₊", x), s)
    dst = map(x->Symbol("dst₊", x), s)
    (;src, dst)
end
function _symvec_to_sym_tup(g::Symmetric, s::AbstractVector{Symbol})
    (;src=s, dst=s)
end
function _symvec_to_sym_tup(g::AntiSymmetric, s::AbstractVector{Symbol})
    (; src=map(x->Symbol("₋", x), s), dst=s)
end
function _symvec_to_sym_tup(g::Directed, dst::AbstractVector{Symbol})
    (;src=Symbol[], dst)
end


_default_name(::Type{VertexFunction}) = :VF
_default_name(::Type{EdgeFunction}) = :EF

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
function _maybewrap!(d, s, T; allownt=false)
    if haskey(d, s)
        v = d[s]
        if v isa T
            d[s] = [v]
        elseif allownt && v isa NamedTuple
            d[s] = map(Iterators.filter(_->true, pairs(v))) do (ntk, ntv)
                if ntv isa T
                    ntk => [ntv]
                else
                    ntk => ntv
                end
            end |> NamedTuple
        end
    end
end

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

# has functions/equality for use in dictionarys
Base.hash(cf::ComponentFunction, h::UInt) = hash_fields(cf, h)
function Base.:(==)(cf1::ComponentFunction, cf2::ComponentFunction)
    typeof(cf1) == typeof(cf2) && equal_fields(cf1, cf2)
end
