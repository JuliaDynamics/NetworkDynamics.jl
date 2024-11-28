"""
    abstract type FeedForwardType end

Abstract supertype for the FeedForwardType traits.
"""
abstract type FeedForwardType end
"""
    PureFeedForward <: FeedForwardType

Trait for component output functions `g` that have pure feed forward behavior
(do not depend on x):

    g!(outs..., ins..., p, t)

See also [`FeedForward`](@ref), [`NoFeedForward`](@ref) and [`PureStateMap`](@ref).
"""
struct PureFeedForward <: FeedForwardType end
"""
    FeedForward <: FeedForwardType

Trait for component output functions `g` that have feed forward behavior. May
depend on everything:

    g!(outs..., x, ins..., p, t)

See also [`PureFeedForward`](@ref), [`NoFeedForward`](@ref) and [`PureStateMap`](@ref).
"""
struct FeedForward <: FeedForwardType end
"""
    NoFeedForward <: FeedForwardType

Trait for component output functions `g` that have no feed forward behavior (do
not depend on inputs):

    g!(outs..., x, p, t)

See also [`PureFeedForward`](@ref), [`FeedForward`](@ref) and [`PureStateMap`](@ref).
"""
struct NoFeedForward <: FeedForwardType end
"""
    PureStateMap <: FeedForwardType

Trait for component output functions `g` that only depends on state:

    g!(outs..., x)

See also [`PureFeedForward`](@ref), [`FeedForward`](@ref) and [`NoFeedForward`](@ref).
"""
struct PureStateMap <: FeedForwardType end

"""
    fftype(x)

Retrieve the feed forward trait of `x`.
"""
function fftype end

"""
    hasfftype(x)

Defaults to false. Musst be overloaded for objects for which `fftype(x)` is defined.
"""
hasfftype(::Any) = false
hasff(x) = fftype(x) isa Union{PureFeedForward, FeedForward}

"""
    StateMask(i::AbstractArray)
    StateMaks(i::Number)

A `StateMask` is a predefined output function. It can be used to define
the output of a component model by picking from the internal state.

I.e. `g=StateMask(2:3)` in a vertex function will output the internal states 2 and 3.
In many contexts, `StateMask`s can be constructed implicitly by just providing
the indices, e.g. `g=1:2`.

For [`EdgeModel`](@ref) this needs to be combined with a [`Directed`](@ref),
[`Symmetric`](@ref), [`AntiSymmetric`](@ref) or [`Fiducial`](@ref) coupling, e.g.
`g=Fiducial(1:2, 3:4)` forwards states 1:2 to dst and states 3:4 to src.
"""
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


abstract type SingleSidedOutputWrapper end

"""
    AntiSymmetric(g_dst)

Wraps a single-sided output function `g_dst` turns it into a double sided
output function which applies

    y_dst = g_dst(...)
    y_src = -y_dst

`g_dst` can be a `Number`/`AbstractArray` to impicitly wrap the corresponding [`StateMask`](@ref).

See also [`Symmetric`](@ref), [`Directed`](@ref), [`Fiducial`](@ref) and [`StateMask`](@ref).
"""
struct AntiSymmetric{G} <: SingleSidedOutputWrapper
    g::G
end
AntiSymmetric(g::Union{AbstractVector,Number}) = AntiSymmetric(StateMask(g))
@inline function (c::AntiSymmetric)(osrc, odst, args...)
    @inline c.g(odst, args...)
    @inbounds for i in 1:length(osrc)
        osrc[i] = -odst[i]
    end
    nothing
end

"""
    Symmetric(g)

Wraps a single-sided output function `g` turns it into a double sided
output function which applies

    y_dst = g(...)
    y_src = y_dst

`g` can be a `Number`/`AbstractArray` to impicitly wrap the corresponding [`StateMask`](@ref).

See also [`AntiSymmetric`](@ref), [`Directed`](@ref), [`Fiducial`](@ref) and [`StateMask`](@ref).
"""
struct Symmetric{G} <: SingleSidedOutputWrapper
    g::G
end
Symmetric(g::Union{AbstractVector,Number}) = Symmetric(StateMask(g))
@inline function (c::Symmetric)(osrc, odst, args...)
    @inline c.g(odst, args...)
    @inbounds for i in 1:length(osrc)
        osrc[i] = odst[i]
    end
    nothing
end

"""
    Directed(g_dst)

Wraps a single-sided output function `g_dst` turns it into a double sided
output function which applies

    y_dst = g_dst(...)

With `Directed` there is no output for the `src` side.
`g_dst` can be a `Number`/`AbstractArray` to impicitly wrap the corresponding [`StateMask`](@ref).

See also [`AntiSymmetric`](@ref), [`Symmetric`](@ref), [`Fiducial`](@ref) and [`StateMask`](@ref).
"""
struct Directed{G} <: SingleSidedOutputWrapper
    g::G
end
Directed(g::Union{AbstractVector,Number}) = Directed(StateMask(g))
@inline function (c::Directed)(osrc, odst, args...)
    @inline c.g(odst, args...)
    nothing
end

"""
    Fiducial(g_src, g_dst)

Wraps two single-sided output function `g_src` and `g_dst` and turns them
into a double sided output function which applies

    y_dst = g_src(...)
    y_src = g_dst(...)

`g` can be a `Number`/`AbstractArray` to impicitly wrap the corresponding [`StateMask`](@ref).

See also [`AntiSymmetric`](@ref), [`Directed`](@ref), [`Fiducial`](@ref) and [`StateMask`](@ref).
"""
struct Fiducial{GS,GD} <: SingleSidedOutputWrapper
    src::GS
    dst::GD
end
function Fiducial(src::Union{AbstractVector,Number}, dst::Union{AbstractVector,Number})
    Fiducial(StateMask(src), StateMask(dst))
end
Fiducial(;src, dst) = Fiducial(src, dst)

@inline function (c::Fiducial)(osrc, odst, args...)
    @inline c.src(osrc, args...)
    @inline c.dst(odst, args...)
    nothing
end

"""
    AnnotatedSym{F}

Wrapper to annotate a vector of symbols as AntiSymmetric, Symmetric or Directed.
"""
struct AnnotatedSym{F}
    s::Vector{Symbol}
end
function AnnotatedSym(wrapper, s::AbstractVector{<:Symbol})
    if wrapper ∉ (AntiSymmetric, Symmetric, Directed)
        throw(ArgumentError("Only AntiSymmetric, Symmetric and Directed are allowed."))
    end
    AnnotatedSym{wrapper}(s)
end
sym(s::AnnotatedSym) = s.s
wrapper(s::AnnotatedSym{W}) where {W} = W
Base.show(io::IO, ::MIME"text/plain", s::AnnotatedSym) = print(io, wrapper(s), "(", sym(s), ")")

"""
    AntiSymmetric(s::AbstractVector{<:Symbol})

Annotate a vector of output-symbols as `AntiSymmetric`, used when creating `EdgeModel`s from
single-sided MTK models.
"""
AntiSymmetric(s::Symbol) = AntiSymmetric([s])
AntiSymmetric(s::AbstractVector{<:Symbol}) = AnnotatedSym(AntiSymmetric, s)
"""
    Symmetric(s::AbstractVector{<:Symbol})

Annotate a vector of output-symbols as `Symmetric`, used when creating `EdgeModel`s from
single-sided MTK models.
"""
Symmetric(s::Symbol) = Symmetric([s])
Symmetric(s::AbstractVector{<:Symbol}) = AnnotatedSym(Symmetric, s)
"""
    Directed(s::AbstractVector{<:Symbol})

Annotate a vector of output-symbols as `Directed`, used when creating `EdgeModel`s from
single-sided MTK models.
"""
Directed(s::Symbol) = Directed([s])
Directed(s::AbstractVector{<:Symbol}) = AnnotatedSym(Directed, s)


abstract type ComponentModel end

struct VertexModel{F,G,FFT,OF,MM} <: ComponentModel
    name::Symbol
    # main function
    f::F
    sym::Vector{Symbol}
    mass_matrix::MM
    # outputs
    g::G
    outsym::Vector{Symbol}
    ff::FFT
    # parameters, optional input sym and optional external inputs
    psym::Vector{Symbol}
    insym::Union{Nothing, Vector{Symbol}}
    extsym::Vector{SymbolicIndex}
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
"""
    VertexModel(; kwargs...)

Build a `VertexModel` according to the keyword arguments.

Main Arguments:
- `f=nothing`: Dynamic function of the component. Can be nothing if `dim` is 0.
- `g`: Output function of the component. Usefull helpers: [`StateMask`](@ref)
- `sym`/`dim`: Symbolic names of the states. If `dim` is provided, `sym` is set automaticially.
- `outsym`/`outdim`:
   Symbolic names of the outputs. If `outdim` is provided, `outsym` is set automaticially.
   Can be infered automaticially if `g` isa `StateMask`.
- psym`/`pdim=0`: Symbolic names of the parameters. If `pdim` is provided, `psym` is set automaticially.
- `mass_matrix=I`: Mass matrix of component. Can be a vector `v` and is then interpreted as `Diagonal(v)`.
- `name=dim>0 ? :VertexM : :StaticVertexM`: Name of the component.

Optional Arguments:
- `insym`/`indim`: Symbolic names of the inputs. If `indim` is provided, `insym` is set automaticially.
- `vidx`: Index of the vertex in the graph, enables graphless constructor.
- `ff`: `FeedForwardType` of component. Will be typically infered from `g` automaticially.
- `obssym`/`obsf`: Define additional "observable" states.
- `symmetadata`/`metadata`: Provide prefilled metadata dictionaries.

All Symbol arguments can be used to set default values, i.e. `psym=[:K=>1, :p]`.
"""
VertexModel(; kwargs...) = _construct_comp(VertexModel, Base.inferencebarrier(kwargs))
VertexModel(v::VertexModel; kwargs...) = _reconstruct_comp(VertexModel, v, Base.inferencebarrier(kwargs))

struct EdgeModel{F,G,FFT,OF,MM} <: ComponentModel
    name::Symbol
    # main function
    f::F
    sym::Vector{Symbol}
    mass_matrix::MM
    # outputs
    g::G
    outsym::@NamedTuple{src::Vector{Symbol},dst::Vector{Symbol}}
    ff::FFT
    # parameters, optional input sym and optional external inputs
    psym::Vector{Symbol}
    insym::Union{Nothing, @NamedTuple{src::Vector{Symbol},dst::Vector{Symbol}}}
    extsym::Vector{SymbolicIndex}
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
"""
    EdgeModel(; kwargs...)

Build a `EdgeModel` according to the keyword arguments.

Main Arguments:
- `f=nothing`: Dynamic function of the component. Can be nothing if `dim` is 0.
- `g`: Output function of the component. Usefull helpers: [`AntiSymmetric`](@ref), [`Symmetric`](@ref), [`Fiducial`](@ref), [`Directed`](@ref) and [`StateMask`](@ref).
- `sym`/`dim`: Symbolic names of the states. If `dim` is provided, `sym` is set automaticially.
- `outsym`/`outdim`:
   Symbolic names of the outputs. If `outdim` is provided, `outsym` is set automaticially.
   In general, outsym for edges isa named tuple `(; src, dst)`. However, depending on the `g` function,
   it might be enough to provide a single vector or even nothing (e.g. `AntiSymmetric(StateMask(1:2))`).
   See [Building `EdgeModel`s](@ref) for examples.
- psym`/`pdim=0`: Symbolic names of the parameters. If `pdim` is provided, `psym` is set automaticially.
- `mass_matrix=I`: Mass matrix of component. Can be a vector `v` and is then interpreted as `Diagonal(v)`.
- `name=dim>0 ? :EdgeM : :StaticEdgeM`: Name of the component.

Optional Arguments:
- `insym`/`indim`: Symbolic names of the inputs. If `indim` is provided, `insym` is set automaticially.
   For edges, `insym` is a named tuple `(; src, dst)`. If give as vector tuple is created automaticially.
- `src`/`dst`: Index or name of the vertices at src and dst end. Enables graphless constructor.
- `ff`: `FeedForwardType` of component. Will be typically infered from `g` automaticially.
- `obssym`/`obsf`: Define additional "observable" states.
- `symmetadata`/`metadata`: Provide prefilled metadata dictionaries.

All Symbol arguments can be used to set default values, i.e. `psym=[:K=>1, :p]`.
"""
EdgeModel(; kwargs...) = _construct_comp(EdgeModel, Base.inferencebarrier(kwargs))
EdgeModel(v::EdgeModel; kwargs...) = _reconstruct_comp(EdgeModel, v, Base.inferencebarrier(kwargs))

"""
    compf(c::ComponentModel)

Retrieve internal function `f` of the component.
"""
compf(c::ComponentModel) = c.f

"""
    compg(c::ComponentModel)

Retrieve output function `g` of the component.
"""
compg(c::ComponentModel) = c.g

"""
    fftype(c::ComponentModel)

Retrieve the feed forward type of the component.
"""
fftype(c::ComponentModel) = c.ff

"""
    dim(c::ComponentModel)::Int

Retrieve the dimension of the component.
"""
dim(c::ComponentModel)::Int = length(sym(c))

"""
    sym(c::ComponentModel)::Vector{Symbol}

Retrieve the symbols of the component.
"""
sym(c::ComponentModel)::Vector{Symbol} = c.sym


"""
    outdim(c::VertexModel)::Int
    outdim(c::EdgeModel)::@NamedTuple(src::Int, dst::Int)

Retrieve the output dimension of the component
"""
outdim(c::VertexModel) = length(outsym(c))
outdim(c::EdgeModel) = (;src=outdim_src(c), dst=outdim_dst(c))

outdim_src(c::EdgeModel) = length(outsym(c).src)
outdim_dst(c::EdgeModel) = length(outsym(c).dst)

"""
   outsym(c::VertexModel)::Vector{Symbol}
   outsym(c::EdgeModel)::@NamedTuple{src::Vector{Symbol}, dst::Vector{Symbol}}

Retrieve the output symbols of the component.
"""
outsym(c::ComponentModel) = c.outsym

"""
    pdim(c::ComponentModel)::Int

Retrieve the parameter dimension of the component.
"""
pdim(c::ComponentModel)::Int = length(psym(c))

"""
    psym(c::ComponentModel)::Vector{Symbol}

Retrieve the parameter symbols of the component.
"""
psym(c::ComponentModel)::Vector{Symbol} = c.psym

"""
    obsf(c::ComponentModel)

Retrieve the observation function of the component.
"""
obsf(c::ComponentModel) = c.obsf

"""
    obssym(c::ComponentModel)::Vector{Symbol}

Retrieve the observation symbols of the component.
"""
obssym(c::ComponentModel)::Vector{Symbol} = c.obssym

"""
    symmetadata(c::ComponentModel)::Dict{Symbol,Dict{Symbol,Any}}

Retrieve the metadata dictionary for the symbols. Keys are the names of the
symbols as they appear in [`sym`](@ref), [`psym`](@ref), [`obssym`](@ref) and [`insym`](@ref).

See also [`symmetadata`](@ref)
"""
symmetadata(c::ComponentModel)::Dict{Symbol,Dict{Symbol,Any}} = c.symmetadata

"""
    metadata(c::ComponentModel)

Retrieve metadata object for the component.

See also [`metadata`](@ref)
"""
metadata(c::ComponentModel)::Dict{Symbol,Any} = c.metadata

"""
    hasinsym(c::ComponentModel)

Checks if the optioan field `insym` is present in the component model.
"""
hasinsym(c::ComponentModel) = !isnothing(c.insym)
"""
    hasindim(c::ComponentModel)

Checks if the optioan field `insym` is present in the component model.
"""
hasindim(c::ComponentModel) = hasinsym(c)

"""
    insym(c::VertexModel)::Vector{Symbol}
    insym(c::EdgeModel)::@NamedTuple{src::Vector{Symbol}, dst::Vector{Symbol}}

Musst be called *after* [`hasinsym`](@ref)/[`hasindim`](@ref) returned true.
Gives the `insym` vector(s). For vertex model just a single vector, for
edges it returns a named tuple `(; src, dst)` with two symbol vectors.
"""
insym(c::VertexModel)::Vector{Symbol} = c.insym
insym(c::EdgeModel)::@NamedTuple{src::Vector{Symbol},dst::Vector{Symbol}} = c.insym

"""
    indim(c::VertexModel)::Int
    indim(c::EdgeModel)::@NamedTuple{src::Int,dst::Int}

Musst be called *after* [`hasinsym`](@ref)/[`hasindim`](@ref) returned true.
Gives the input dimension(s).
"""
indim(c::VertexModel)::Int = length(insym(c))
indim(c::EdgeModel)::@NamedTuple{src::Int,dst::Int} = (; src=length(insym(c).src), dst=length(insym(c).dst))

"""
    extsym(c::ComponentModel)::Vector{Symbol}

Retrieve the external input symbols of the component.
"""
extsym(c::ComponentModel) = c.extsym

"""
    extdim(c::ComponentModel)::Int

Retrieve the external input dimension of the component.
"""
extdim(c::ComponentModel) = length(extsym(c))

# return both "observed" outputs (those that do not shadow states) and true observed
outsym_flat(c::ComponentModel) = c._outsym_flat
obssym_all(c::ComponentModel) = c._obssym_all

# normalized means, that we'll always return a tuple of values, thus we can generalize better over Edges/Vertices
insym_normalized(c::EdgeModel) = values(insym(c))
insym_normalized(c::VertexModel) = (insym(c),)
indim_normalized(c::ComponentModel) = map(length, insym_normalized(c))
outsym_normalized(c::EdgeModel) = values(outsym(c))
outsym_normalized(c::VertexModel) = (outsym(c),)
outdim_normalized(c::ComponentModel) = map(length, outsym_normalized(c))

infer_fftype(::Type{<:VertexModel}, g, dim, hasext) = _infer_fftype(g, 1, 1+hasext, dim)
infer_fftype(::Type{<:EdgeModel}, g, dim, hasext) = _infer_fftype(g, 2, 2+hasext, dim)
# special cases for wrapped output functions
function infer_fftype(::Type{<:EdgeModel}, g::Union{Symmetric, AntiSymmetric, Directed}, dim, hasext)
    _infer_fftype(g.g, 1, 2+hasext, dim)
end
function infer_fftype(::Type{<:EdgeModel}, g::Fiducial, dim, hasext)
    ffsrc = _infer_fftype(g.src, 1, 2+hasext, dim)
    ffdst = _infer_fftype(g.dst, 1, 2+hasext, dim)
    if ffsrc == ffdst
        return ffsrc
    else
        error("Both src and dst coupling functions have different ff types.")
    end
end
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
    dispatchT(<:ComponentModel) :: Type{<:ComponentModel}

Returns the type "essence" of the component used for dispatch.
Fills up type parameters with `nothing` to ensure `Core.Compiler.isconstType`
for GPU compatibility.
"""
dispatchT(::T) where {T<:ComponentModel} = dispatchT(T)
dispatchT(T::Type{<:VertexModel}) = VertexModel{nothing,nothing,nothing,nothing,nothing}
dispatchT(T::Type{<:EdgeModel}) = EdgeModel{nothing,nothing,nothing,nothing,nothing}

# TODO: introduce batchequal hash for faster batching of component models
batchequal(a, b) = false
function batchequal(a::ComponentModel, b::ComponentModel)
    compf(a)  == compf(b)  || return false
    compg(a)  == compg(b)  || return false
    fftype(a) == fftype(b) || return false
    dim(a)    == dim(b)    || return false
    outdim(a) == outdim(b) || return false
    pdim(a)   == pdim(b)   || return false
    extdim(a) == extdim(b) || return false
    return true
end

"""
    _construct_comp(::Type{T}, kwargs) where {T}

Internal function to construct a component model from keyword arguments.
Fills up kw arguments with default values and performs sanity checks.
"""
function _construct_comp(::Type{T}, @nospecialize(kwargs)) where {T}
    dict = _fill_defaults(T, Base.inferencebarrier(kwargs))

    # check signature of f
    # if !_valid_signature(T, dict[:f])
    #     throw(ArgumentError("Function f does not take the correct number of arguments."))
    # end

    # pop check keyword
    check = pop!(dict, :check, CHECK_COMPONENT[])

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

function _reconstruct_comp(::Type{T}, cf::ComponentModel, kwargs) where {T}
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
    allow_output_sym_clash = pop!(dict, :allow_output_sym_clash, false)

    # syms might be provided as single pairs or symbols, wrap in vector
    _maybewrap!(dict, :sym, Union{Symbol, Pair})
    _maybewrap!(dict, :psym, Union{Symbol, Pair})
    _maybewrap!(dict, :insym, Union{Symbol,Pair}; allownt=T <: EdgeModel)
    _maybewrap!(dict, :outsym, Union{Symbol,Pair}; allownt=T <: EdgeModel)
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
    if haskey(dict, :vidx) && T <: VertexModel
        vidx = pop!(dict, :vidx)
        metadata[:graphelement] = vidx
    end
    if haskey(dict, :src) && haskey(dict, :dst) && T <: EdgeModel
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
    if T <: EdgeModel && g isa StateMask
        throw(ArgumentError("StateMask cannot be used directly for EdgeModel. Wrap in Fiducial, Symmetric, AntiSymmetric or Directed instead."))
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
        if T <: VertexModel
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
    if T <: EdgeModel && outsym isa AbstractVector
        if _has_metadata(outsym)
            outsym, _metadata = _split_metadata(outsym)
            mergewith!(merge!, symmetadata, _metadata)
        end
        outsym = dict[:outsym] = _symvec_to_sym_tup(g, outsym)
    end

    if haskey(dict, :outdim)
        _outdim = pop!(dict, :outdim)
        if T <: EdgeModel
            if _outdim isa NamedTuple && (_outdim.src != length(outsym.src) || _outdim.dst != length(outsym.dst))
                throw(ArgumentError("Length of outsym and outdim must match."))
            elseif _outdim != length(outsym.dst) || !(g isa Directed) && _outdime != length(outsym.src)
                throw(ArgumentError("Length of outsym and outdim must match."))
            end
        elseif T <: VertexModel
            if _outdim != length(outsym)
                throw(ArgumentError("Length of outsym and outdim must match."))
            end
        end
    end

    # extract and merge metadata from outsym
    if T <: VertexModel
        if _has_metadata(outsym)
            dict[:outsym], _metadata = _split_metadata(outsym)
            mergewith!(merge!, symmetadata, _metadata)
        end
    elseif T <: EdgeModel
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
    outsym = dict[:outsym]

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
    if T <: EdgeModel && insym isa AbstractVector
        if _has_metadata(insym)
            throw(ArgumentError("Metadata for insyms can only be provided when using full (;src, dst) form"))
        end
        insym = dict[:insym] = _symvec_to_sym_tup(nothing, insym)
    end
    if haskey(dict, :indim)
        _indim = pop!(dict, :indim)
        if T <: EdgeModel
            if _indim isa NamedTuple && (_indim.src != length(insym.src) || _indim.dst != length(insym.dst))
                throw(ArgumentError("Length of insym and indim must match."))
            elseif _indim != length(insym.src) || _indime != length(insym.dst)
                throw(ArgumentError("Length of insym and indim must match."))
            end
        elseif T <: VertexModel
            if _indim != length(outsym)
                throw(ArgumentError("Length of insym and indim must match."))
            end
        end
    end
    # extract and merge metadata from insym
    if !isnothing(insym)
        if T <: VertexModel
            if _has_metadata(insym)
                dict[:insym], _metadata = _split_metadata(insym)
                mergewith!(merge!, symmetadata, _metadata)
            end
        elseif T <: EdgeModel
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
        dict[:name] = if T <: VertexModel
            dim > 0 ? :VertexM : :StaticVertexM
        elseif T <: EdgeModel
            dim > 0 ? :EdgeM : :StaticEdgeM
        else
            error()
        end
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
    _outsym_flat = if T <: VertexModel
        outsym
    elseif T <: EdgeModel
        vcat(outsym.src, outsym.dst)
    else
        error()
    end
    dict[:_outsym_flat] = _outsym_flat
    dict[:_obssym_all] = setdiff(_outsym_flat, sym) ∪ obssym

    ####
    #### Extsym
    ####
    if haskey(dict, :extsym)
        @assert dict[:extsym] isa Vector{<:SymbolicIndex}
    else
        dict[:extsym] = SymbolicIndex[]
    end

    # infer fftype (needs dim and extdim)
    if !haskey(dict, :ff)
        dict[:ff] = if hasfftype(g)
            fftype(g)
        else
            infer_fftype(T, g, dim, !isempty(dict[:extsym]))
        end
    end

    # check for name clashes (at the end because only now sym, psym, obssym are initialized)
    _s  = sym
    _ps = psym
    _obss = obssym
    _os = if allow_output_sym_clash || _has_sym_to_outsym_mapping(g)
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
_has_sym_to_outsym_mapping(::Any) = false
_has_sym_to_outsym_mapping(::StateMask) = true
_has_sym_to_outsym_mapping(::Directed{<:StateMask}) = true
_has_sym_to_outsym_mapping(::AntiSymmetric{<:StateMask}) = true
_has_sym_to_outsym_mapping(::Symmetric{<:StateMask}) = true
_has_sym_to_outsym_mapping(::Fiducial{<:StateMask, <:StateMask}) = true

_sym_to_outsym(g::StateMask, s::AbstractVector{Symbol}) = s[g.idxs]
function _sym_to_outsym(g::AntiSymmetric{<:StateMask}, s::AbstractVector{Symbol})
    s = _sym_to_outsym(g.g, s)
    _symvec_to_sym_tup(g, s)
end
function _sym_to_outsym(g::Symmetric{<:StateMask}, s::AbstractVector{Symbol})
    s = _sym_to_outsym(g.g, s)
    _symvec_to_sym_tup(g, s)
end
function _sym_to_outsym(g::Directed{<:StateMask}, s::AbstractVector{Symbol})
    s = _sym_to_outsym(g.g, s)
    _symvec_to_sym_tup(g, s)
end
function _sym_to_outsym(g::Fiducial{<:StateMask, <:StateMask}, s::AbstractVector{Symbol})
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
    copy(c::NetworkDynamics.ComponentModel)

Shallow copy of the component model. Creates a deepcopy of `metadata` and `symmetadata`
but references the same objects everywhere else.
"""
@generated function Base.copy(c::ComponentModel)
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
Base.hash(cf::ComponentModel, h::UInt) = hash_fields(cf, h)
function Base.:(==)(cf1::ComponentModel, cf2::ComponentModel)
    typeof(cf1) == typeof(cf2) && equal_fields(cf1, cf2)
end

function compfg(c)
    (outs, du, u, ins, p, t) -> begin
        f = compf(c)
        isnothing(f) || f(du, u, ins..., p, t)
        compg(c)(_gargs(fftype(c), outs, du, u, ins, p, t)...)
        nothing
    end
end
_gargs(::PureFeedForward, outs, du, u, ins, p, t) = (outs..., ins..., p, t)
_gargs(::FeedForward, outs, du, u, ins, p, t) = (outs..., u, ins..., p, t)
_gargs(::NoFeedForward, outs, du, u, ins, p, t) = (outs..., u, p, t)
_gargs(::PureStateMap, outs, du, u, ins, p, t) = (outs..., u)

"""
    ff_to_constraint(v::VertexModel)

Takes `VertexModel` `v` with feed forward and turns all algebraic output
states into internal states by defining algebraic constraints
contraints `0 = out - g(...)`. The new output function is
just a [`StateMask`](@ref) into the extended internal state vector.

Returns the transformed `VertexModel`.
"""
function ff_to_constraint(v::VertexModel)
    hasff(v) || error("Vertex does not have feed forward property.")
    # TODO: ff_to_constraint with ext einput
    has_external_input(v) && error("ff_to_constraint for model with external input not implemented.")

    olddim = dim(v)
    odim = outdim(v)
    newf = _newf(fftype(v), v.f, v.g, olddim)
    newg = StateMask(olddim+1:olddim+odim)

    mass_matrix = [v.mass_matrix zeros(olddim, odim)
                   zeros(odim, olddim) zeros(odim, odim)]
    if LinearAlgebra.isdiag(v.mass_matrix)
        mass_matrix = LinearAlgebra.Diagonal(mass_matrix)
    end

    isnothing(v.obsf) || @warn "Observed function might be broke due to ff_to_constraint conversion."
    newsym = vcat(sym(v), outsym(v))
    VertexModel(; f=newf, g=newg, sym=newsym, mass_matrix, name=v.name,
        psym=v.psym, obsf=v.obsf, obssym=v.obssym, insym=v.insym,
        metadata=v.metadata, symmetadata=v.symmetadata)
end

function _newf(::PureFeedForward, f, g, dim)
    @closure (du, u, esum, p, t) -> begin
        if !isnothing(f)
            duf = @view du[1:dim]
            uf = @view u[1:dim]
            f(duf, uf, esum, p, t)
        end
        dug = @view du[dim+1:end]
        g(dug, esum, p, t)
        ug = @view u[dim+1:end]
        dug .= dug .- ug
        nothing
    end
end
function _newf(::FeedForward, f, g, dim)
    @closure (du, u, esum, p, t) -> begin
        if !isnothing(f)
            duf = @view du[1:dim]
            uf = @view u[1:dim]
            f(duf, uf, esum, p, t)
        end
        dug = @view du[dim+1:end]
        g(dug, uf, esum, p, t)
        ug = @view u[dim+1:end]
        dug .= dug .- ug
        nothing
    end
end
