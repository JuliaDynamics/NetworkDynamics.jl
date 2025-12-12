const VALID_NUMERIC_SUB_IDX = Union{Colon, Int, AbstractVector{<:Int}, NTuple{<:Any, Int}}

abstract type NumericSubIndex{T} end
"""
    ParamIdx{T} <: NumericSubIndex{T}

Wrapper type for indexing into parameters by index rather than symbol.
`VPIndex(:name, 1)` is equivalent to `VIndex(:name, ParamIdx(1))` and tells
NetworkDynamics that you mean the first parameter of the component in contrast
to the first state of the component.
Similarly, `ParamIdx(:)`, `ParamIdx(1:3)` or `ParamIdx([1,2,3])` can be used
to access multiple parameters by numeric index at once.

See also: [`StateIdx`](@ref), [`VIndex`](@ref), [`EIndex`](@ref)
"""
struct ParamIdx{T<:VALID_NUMERIC_SUB_IDX} <: NumericSubIndex{T}
    idx::T
    function ParamIdx(i)
        _i = transform_numeric_subidx(i)
        new{typeof(_i)}(_i)
    end
end
"""
    StateIdx{T} <: NumericSubIndex{T}

Wrapper type for indexing into states by index rather than symbol.
`VIndex(:name, 1)` is equivalent to `VIndex(:name, StateIdx(1))` and tells
NetworkDynamics that you mean the first state of the component in contrast
to the first parameter of the component.
Similarly, `StateIdx(:)`, `StateIdx(1:3)` or `StateIdx([1,2,3])` can be used
to access multiple states by numeric index at once.

See also: [`ParamIdx`](@ref), [`VIndex`](@ref), [`EIndex`](@ref)
"""
struct StateIdx{T<:VALID_NUMERIC_SUB_IDX} <: NumericSubIndex{T}
    idx::T
    function StateIdx(i)
        _i = transform_numeric_subidx(i)
        new{typeof(_i)}(_i)
    end
end
function transform_numeric_subidx(i)
    if i isa VALID_NUMERIC_SUB_IDX
        return i
    elseif i isa AbstractVector && all(x -> x isa Int, i)
        return convert(Vector{Int}, i)
    else
        throw(ArgumentError("StateIdx/ParamIdx only supports Colon, Int or collections of int, got $i"))
    end
end
idxtype(::ParamIdx) = ParamIdx
idxtype(::StateIdx) = StateIdx

# ParamIdx(1:2) == ParamIdx([1,2])
function Base.:(==)(a::NumericSubIndex, b::NumericSubIndex)
    idxtype(a) == idxtype(b) && a.idx == b.idx
end

const ITERATABLE_NUMERIC_SUB_IDX = NumericSubIndex{<:Union{AbstractVector, Tuple}}

Base.length(ni::ITERATABLE_NUMERIC_SUB_IDX) = length(ni.idx)
Base.size(ni::ITERATABLE_NUMERIC_SUB_IDX) = (length(ni),)
Base.IteratorSize(ni::ITERATABLE_NUMERIC_SUB_IDX) = Base.HasShape{1}()
Base.broadcastable(ni::ITERATABLE_NUMERIC_SUB_IDX) = ni
Base.ndims(::Type{<:ITERATABLE_NUMERIC_SUB_IDX}) = 1
Base.axes(ni::ITERATABLE_NUMERIC_SUB_IDX) = axes(ni.idx)
Base.getindex(ni::ITERATABLE_NUMERIC_SUB_IDX, i) = idxtype(ni)(ni.idx[i])
Base.eltype(ni::ITERATABLE_NUMERIC_SUB_IDX) = idxtype(ni){eltype(ni.idx)}
function Base.iterate(ni::ITERATABLE_NUMERIC_SUB_IDX, state=nothing)
    it = isnothing(state) ? iterate(ni.idx) : iterate(ni.idx, state)
    isnothing(it) && return nothing
    idxtype(ni)(it[1]), it[2]
end

# dont wrap empty tuples
wrap_sidx(x::Tuple{}) = x
wrap_pidx(x::Tuple{}) = x
wrap_sidx(x::VALID_NUMERIC_SUB_IDX) = StateIdx(x)
wrap_pidx(x::VALID_NUMERIC_SUB_IDX) = ParamIdx(x)
function wrap_sidx(x)
    if (x isa AbstractVector || x isa Tuple) && !isempty(x) && any(i -> i isa Int, x)
        map(wrap_sidx, x)
    else
        x
    end
end
function wrap_pidx(x)
    if (x isa AbstractVector || x isa Tuple) && !isempty(x) && any(i -> i isa Int, x)
        map(wrap_pidx, x)
    else
        x
    end
end

"""
    VIndex{C,S} <: SymbolicIndex{C,S}
    idx = VIndex(comp, sub)

A symbolic index for a vertex variable.
- `comp`: the component index, either int, symbol or a collection
- `sub`: the subindex, either int, symbol or a collection of those.

Symbolic indices are rather flexible and can potentially point to states, parameters, inputs, outputs, or observables.

The most basic form is `VIndex(1, :P)` which points to the variable with the name `:P` in the first vertex model.
The component can be also given by unique name, so `VIndex(:a, :P)` would point to the vertex with unique name `:a`.

It is also possible to have "collections" of indices, such as
```
VIndex(1:5, 1)     # first state of vertices 1 to 5
VIndex(7, (:x,:y)) # states :x and :y of vertex 7
```

They can be used to index into objects supporting the `SymbolicIndexingInterface`,
e.g. [`NWState`](@ref), [`NWParameter`](@ref) or `ODESolution`.

It is also possible to construct vertices without a sub index, in which case they point
to a component rather than a specific variable:
```
VIndex(2)          # references the second vertex model
VIndex(:a)         # references vertex with unique name :a
```
For example `nw[VIndex(2)]` would return the 2nd [`VertexModel`](@ref) in the `Network`.

See also: [`EIndex`](@ref), [`StateIdx`](@ref), [`ParamIdx`](@ref), [`generate_indices`](@ref), [`NWState`](@ref), [`NWParameter`](@ref)
"""
struct VIndex{C,S} <: SymbolicIndex{C,S}
    compidx::C
    subidx::S
    function VIndex(c, s)
        _s = wrap_sidx(s)
        new{typeof(c),typeof(_s)}(c, _s)
    end
end
VIndex(ci) = VIndex(ci, nothing)
"""
    EIndex{C,S} <: SymbolicIndex{C,S}
    idx = EIndex(comp, sub)

A symbolic index for an edge variable.
- `comp`: the component index, either int, symbol, pair or a collection
- `sub`: the subindex, either int, symbol or a collection of those.

Symbolic indices are rather flexible and can potentially point to states, parameters, inputs, outputs, or observables.

The most basic form is `EIndex(1, :P)` which points to the variable with the name `:P` in the first edge model.
The component can be also given by unique name, so `EIndex(:a, :P)` would point to the edge with unique name `:a`.
For edges, the component can also be specified as a source-destination pair:

```
EIndex(1=>2, :P)    # variable :P in edge from vertex 1 to vertex 2
EIndex(:a=>:b, :P)  # variable :P in edge from vertex :a to vertex :b
```

It is also possible to have "collections" of indices, such as
```
EIndex(1:5, 1)     # first state of edges 1 to 5
EIndex(7, (:x,:y)) # states :x and :y of edge 7
```

They can be used to index into objects supporting the `SymbolicIndexingInterface`,
e.g. [`NWState`](@ref), [`NWParameter`](@ref) or `ODESolution`.

It is also possible to construct edges without a sub index, in which case they point
to a component rather than a specific variable:
```
EIndex(2)          # references the second edge model
EIndex(1=>2)       # references edge from v1 to v2
EIndex(:a=>:b)     # references edge from vertex :a to vertex :b
```
For example `nw[EIndex(2)]` would return the 2nd [`EdgeModel`](@ref) in the `Network`.

See also: [`VIndex`](@ref), [`StateIdx`](@ref), [`ParamIdx`](@ref), [`generate_indices`](@ref), [`NWState`](@ref), [`NWParameter`](@ref)
"""
struct EIndex{C,S} <: SymbolicIndex{C,S}
    compidx::C
    subidx::S
    function EIndex(c, s)
        _s = wrap_sidx(s)
        new{typeof(c),typeof(_s)}(c, _s)
    end
end
EIndex(ci) = EIndex(ci, nothing)

# VIndex(1:2, :foo) == VIndex([1,2], :foo) and so on
function Base.:(==)(a::SymbolicIndex, b::SymbolicIndex)
    idxtype(a) == idxtype(b) && a.compidx == b.compidx && a.subidx == b.subidx
end

"""
    VPIndex(c, s)

Backward compatibility constructor for vertex parameter indices.
Wraps integer-like subindices in `ParamIdx`, leaves other types unchanged.
Equivalent to `VIndex(c, ParamIdx(s))` when `s` is integer-like, otherwise `VIndex(c, s)`.

See also: [`VIndex`](@ref), [`ParamIdx`](@ref)
"""
VPIndex(c, s) = VIndex(c, wrap_pidx(s))
VPIndex(c) = VIndex(c)

"""
    EPIndex(c, s)

Backward compatibility constructor for edge parameter indices.
Wraps integer-like subindices in `ParamIdx`, leaves other types unchanged.
Equivalent to `EIndex(c, ParamIdx(s))` when `s` is integer-like, otherwise `EIndex(c, s)`.

See also: [`EIndex`](@ref), [`ParamIdx`](@ref)
"""
EPIndex(c, s) = EIndex(c, wrap_pidx(s))
EPIndex(c) = EIndex(c)

idxtype(s::VIndex) = VIndex
idxtype(s::EIndex) = EIndex

const CONCRETE_COMPIDX = Union{<:Pair, Int, Symbol}
const CONCRETE_SUBIDX = Union{Symbol, NumericSubIndex{Int}}

SII.symbolic_type(::Type{<:SymbolicIndex{<:CONCRETE_COMPIDX,<:CONCRETE_SUBIDX}}) = SII.ScalarSymbolic()
SII.symbolic_type(::Type{<:SymbolicIndex}) = SII.ArraySymbolic()

SII.hasname(::SymbolicIndex) = false
SII.hasname(::SymbolicIndex{<:CONCRETE_COMPIDX,<:CONCRETE_SUBIDX}) = true
function SII.getname(x::VIndex)
    prefix = x.compidx isa Symbol ? Symbol() : :v
    Symbol(prefix, Symbol(x.compidx), :₊, Symbol(x.subidx))
end
function SII.getname(x::EIndex)
    if x.compidx isa Pair
        src, dst = x.compidx
        _src = src isa Int ? Symbol(:v, src) : Symbol(src)
        _dst = dst isa Int ? Symbol(:v, dst) : Symbol(dst)
        Symbol(_src, "ₜₒ", _dst, :₊, _symbol_repr(x.subidx))
    else
        prefix = x.compidx isa Int ? :e : Symbol()
        Symbol(prefix, Symbol(x.compidx), :₊, _symbol_repr(x.subidx))
    end
end
_symbol_repr(x) = Symbol(x)
_symbol_repr(paramidx::ParamIdx) = Symbol(:p, paramidx.idx)
_symbol_repr(Stateidx::StateIdx) = Symbol(:x, Stateidx.idx)

resolvecompidx(nw::Network, sni) = resolvecompidx(nw.im, sni)
resolvecompidx(::IndexManager, sni::SymbolicIndex{Int}) = sni.compidx
function resolvecompidx(im::IndexManager, sni::SymbolicIndex{Symbol})
    dict = sni isa VIndex ? im.unique_vnames : im.unique_enames
    if haskey(dict, sni.compidx)
        return dict[sni.compidx]
    else
        throw(ArgumentError("Could not resolve component index for $sni, the name might not be unique?"))
    end
end
function resolvecompidx(im::IndexManager, sni::EIndex{<:Pair})
    src, dst = sni.compidx

    src_i = try
        resolvecompidx(im, VIndex(src))
    catch
        throw(ArgumentError("Could not resolve edge source $src"))
    end
    dst_i = try
        resolvecompidx(im, VIndex(dst))
    catch
        throw(ArgumentError("Could not resolve edge destination $dst"))
    end

    eidx = findfirst(im.edgevec) do e
        e.src == src_i && e.dst == dst_i
    end
    if isnothing(eidx)
        reverse = findfirst(im.edgevec) do e
            e.src == dst_i && e.dst == src_i
        end
        err = "Invalid Index: Network does not contain edge from $(src) => $(dst)!"
        if !isnothing(reverse)
            err *= " Maybe you meant the reverse edge from $(dst) => $(src)?"
        end
        throw(ArgumentError(err))
    end
    return eidx
end
getcomp(nw::Network, sni) = getcomp(nw.im, sni)
getcomp(im::IndexManager, sni::EIndex) = im.edgem[resolvecompidx(im, sni)]
getcomp(im::IndexManager, sni::VIndex) = im.vertexm[resolvecompidx(im, sni)]

getcomprange(nw::Network, sni) = getcomprange(nw.im, sni)
getcomprange(im::IndexManager, sni::VIndex{<:Union{Symbol,Int}}) = im.v_data[resolvecompidx(im, sni)]
getcomprange(im::IndexManager, sni::EIndex{<:CONCRETE_COMPIDX}) = im.e_data[resolvecompidx(im, sni)]

getcompoutrange(nw::Network, sni) = getcompoutrange(nw.im, sni)
getcompoutrange(im::IndexManager, sni::VIndex{<:Union{Symbol,Int}}) = im.v_out[resolvecompidx(im, sni)]
getcompoutrange(im::IndexManager, sni::EIndex{<:CONCRETE_COMPIDX}) = flatrange(im.e_out[resolvecompidx(im, sni)])

getcompprange(nw::Network, sni::VIndex{<:Union{Symbol,Int}}) = nw.im.v_para[resolvecompidx(nw, sni)]
getcompprange(nw::Network, sni::EIndex{<:CONCRETE_COMPIDX}) = nw.im.e_para[resolvecompidx(nw, sni)]

subsym_has_idx(sym::Symbol, syms) = sym ∈ syms
subsym_has_idx(idx::NumericSubIndex{Int}, syms) = 1 ≤ idx.idx ≤ length(syms)
subsym_to_idx(sym::Symbol, syms) = findfirst(isequal(sym), syms)
subsym_to_idx(idx::NumericSubIndex{Int}, _) = idx.idx

####
#### Iterator/Broadcast interface for ArraySymbolic types
####
# TODO: not broadcasting over idx with colon is weird
# Base.broadcastable(si::SymbolicIndex{<:Union{CONCRETE_COMPIDX,Colon},<:Union{CONCRETE_SUBIDX,Colon}}) = Ref(si)
Base.broadcastable(si::SymbolicIndex{<:CONCRETE_COMPIDX,<:CONCRETE_SUBIDX}) = Ref(si)

const _IterableComponent = SymbolicIndex{<:Union{AbstractVector,Tuple},<:Union{CONCRETE_SUBIDX,Nothing}}
Base.length(si::_IterableComponent) = length(si.compidx)
Base.size(si::_IterableComponent) = (length(si),)
Base.IteratorSize(si::_IterableComponent) = Base.HasShape{1}()
Base.broadcastable(si::_IterableComponent) = si
Base.ndims(::Type{<:_IterableComponent}) = 1
Base.axes(si::_IterableComponent) = axes(si.compidx)
Base.getindex(si::_IterableComponent, i) = idxtype(si)(si.compidx[i], si.subidx)
function Base.eltype(si::_IterableComponent)
    if isconcretetype(eltype(si.compidx))
        idxtype(si){eltype(si.compidx),typeof(si.subidx)}
    else
        Any
    end
end
function Base.iterate(si::_IterableComponent, state=nothing)
    it = isnothing(state) ? iterate(si.compidx) : iterate(si.compidx, state)
    isnothing(it) && return nothing
    idxtype(si)(it[1], si.subidx), it[2]
end

const _IterableSubcomponent = SymbolicIndex{<:CONCRETE_COMPIDX,<:Union{AbstractVector,Tuple,ITERATABLE_NUMERIC_SUB_IDX}}
Base.length(si::_IterableSubcomponent) = length(si.subidx)
Base.size(si::_IterableSubcomponent) = (length(si),)
Base.IteratorSize(si::_IterableSubcomponent) = Base.HasShape{1}()
Base.broadcastable(si::_IterableSubcomponent) = si
Base.ndims(::Type{<:_IterableSubcomponent}) = 1
Base.axes(si::_IterableSubcomponent) = axes(si.subidx)
Base.getindex(si::_IterableSubcomponent, i) = idxtype(si)(si.compidx, si.subidx[i])
function Base.eltype(si::_IterableSubcomponent)
    if isconcretetype(eltype(si.subidx))
        idxtype(si){eltype(si.compidx),eltype(si.subidx)}
    else
        Any
    end
end
function Base.iterate(si::_IterableSubcomponent, state=nothing)
    it = isnothing(state) ? iterate(si.subidx) : iterate(si.subidx, state)
    isnothing(it) && return nothing
    idxtype(si)(si.compidx, it[1]), it[2]
end

_hascolon(i) = false
_hascolon(::SymbolicIndex{C,S}) where {C,S} = C === Colon || S <: NumericSubIndex{Colon}
_resolve_colon(nw::Network, sni::VIndex{Colon, <:CONCRETE_SUBIDX}) = VIndex(1:nv(nw), sni.subidx)
_resolve_colon(nw::Network, sni::EIndex{Colon, <:CONCRETE_SUBIDX}) = EIndex(1:ne(nw), sni.subidx)
function _resolve_colon(nw::Network, sni::SymbolicIndex{<:Union{Symbol,Int},<:StateIdx{Colon}})
    idxtype(sni)(sni.compidx, StateIdx(1:dim(getcomp(nw, sni))))
end
function _resolve_colon(nw::Network, sni::SymbolicIndex{<:Union{Symbol,Int},<:ParamIdx{Colon}})
    idxtype(sni)(sni.compidx, ParamIdx(1:pdim(getcomp(nw, sni))))
end
_resolve_colon(nw::Network, idx) = throw(ArgumentError("Cannot resolve colons for both component and subindex, got $idx"))


#### Implmentation of index provider interface
####
#### Structural things
####
SII.is_independent_variable(nw::Network, sym) = sym == :t
SII.independent_variable_symbols(nw::Network) = [:t]
SII.is_time_dependent(nw::Network) = true
SII.constant_structure(::Network) = true
SII.all_variable_symbols(nw::Network) = vcat(SII.variable_symbols(nw), observed_symbols(nw))
SII.all_symbols(nw::Network) = vcat(SII.all_variable_symbols(nw), SII.parameter_symbols(nw))


####
#### variable indexing
####
const POTENTIAL_SCALAR_SIDX = Union{SymbolicIndex{<:CONCRETE_COMPIDX,<:Union{StateIdx{Int},Symbol}}}
function SII.is_variable(nw::Network, sni)
    if _hascolon(sni)
        SII.is_variable(nw, _resolve_colon(nw,sni))
    elseif SII.symbolic_type(sni) === SII.ArraySymbolic()
        all(Base.Fix1(SII.is_variable, nw), sni)
    else
        _is_variable(nw, sni)
    end
end
_is_variable(nw::Network, sni) = false
function _is_variable(nw::Network, sni::POTENTIAL_SCALAR_SIDX)
    cf = getcomp(nw, sni)
    return subsym_has_idx(sni.subidx, sym(cf))
end

function SII.variable_index(nw::Network, sni)
    if _hascolon(sni)
        SII.variable_index(nw, _resolve_colon(nw,sni))
    elseif sni isa AbstractVector || SII.symbolic_type(sni) === SII.ArraySymbolic()
        SII.variable_index.(nw, sni)
    else
        _variable_index(nw, sni)
    end
end
function _variable_index(nw::Network, sni::POTENTIAL_SCALAR_SIDX)
    cf = getcomp(nw, sni)
    range = getcomprange(nw, sni)
    idx = subsym_to_idx(sni.subidx, sym(cf))
    isnothing(idx) ? nothing : range[idx]
end

"""
    SymbolicIndexingInterface.variable_symbols(nw::Network)

Returns a vector of all symbolic network indices which in the same order as
the flat state vector.

See also: [`NWState`](@ref), [`uflat`](@ref).
"""
function SII.variable_symbols(nw::Network)
    syms = Vector{SymbolicIndex{Int,Symbol}}(undef, dim(nw))
    for (ci, cf) in pairs(nw.im.vertexm)
        syms[nw.im.v_data[ci]] .= VIndex.(ci, sym(cf))
    end
    for (ci, cf) in pairs(nw.im.edgem)
        syms[nw.im.e_data[ci]] .= EIndex.(ci, sym(cf))
    end
    return syms
end


####
#### parameter indexing
####
# when using an number instead of symbol only PIndex is valid
const POTENTIAL_SCALAR_PIDX = Union{SymbolicIndex{<:CONCRETE_COMPIDX,<:Union{ParamIdx{Int},Symbol}}}
function SII.is_parameter(nw::Network, sni)
    if _hascolon(sni)
        SII.is_parameter(nw, _resolve_colon(nw,sni))
    elseif SII.symbolic_type(sni) === SII.ArraySymbolic()
        all(Base.Fix1(SII.is_parameter, nw), sni)
    else
        _is_parameter(nw, sni)
    end
end
_is_parameter(nw::Network, sni) = false
function _is_parameter(nw::Network, sni::POTENTIAL_SCALAR_PIDX)
    cf = getcomp(nw, sni)
    return subsym_has_idx(sni.subidx, psym(cf))
end

function SII.parameter_index(nw::Network, sni)
    if _hascolon(sni)
        SII.parameter_index(nw, _resolve_colon(nw,sni))
    elseif sni isa AbstractVector || SII.symbolic_type(sni) === SII.ArraySymbolic()
        SII.parameter_index.(nw, sni)
    else
        _parameter_index(nw, sni)
    end
end
function _parameter_index(nw::Network, sni::POTENTIAL_SCALAR_PIDX)
    cf = getcomp(nw, sni)
    range = getcompprange(nw, sni)
    idx = subsym_to_idx(sni.subidx, psym(cf))
    isnothing(idx) ? nothing : range[idx]
end

"""
    SymbolicIndexingInterface.parameter_symbols(nw::Network)

Returns a vector of all symbolic network indices which in the same order as
the flat parameter vector.

See also: [`NWParameter`](@ref), [`NWState`](@ref), [`pflat`](@ref).
"""
function SII.parameter_symbols(nw::Network)
    syms = Vector{SymbolicIndex{Int,Symbol}}(undef, pdim(nw))
    for (ci, cf) in pairs(nw.im.vertexm)
        syms[nw.im.v_para[ci]] .= VIndex.(ci, psym(cf))
    end
    for (ci, cf) in pairs(nw.im.edgem)
        syms[nw.im.e_para[ci]] .= EIndex.(ci, psym(cf))
    end
    return syms
end

####
#### Timeseries parameter indexing
####
const DEFAULT_PARA_TS_IDX = 1
SII.is_timeseries_parameter(nw::Network, sni) = SII.is_parameter(nw::Network, sni)
function SII.timeseries_parameter_index(nw::Network, sni)
    # NOTE: ALL parameters are lumped in timeseries with idx 1
    SII.ParameterTimeseriesIndex.(DEFAULT_PARA_TS_IDX, SII.parameter_index.(nw, sni))
end

function SII.get_all_timeseries_indexes(nw::Network, sym)
    # allways return ContTimeseries if sym is random stuff, see
    # https://github.com/SciML/SymbolicIndexingInterface.jl/issues/95

    # if SII.is_variable(nw, sym) || SII.is_independent_variable(nw, sym) || SII.is_observed(nw, sym)
    #     return Set([SII.ContinuousTimeseries()])
    # elseif SII.is_timeseries_parameter(nw, sym)
    #     return Set([SII.timeseries_parameter_index(nw, sym).timeseries_idx])
    # else
    #     return Set()
    # end
    if !iszero(pdim(nw)) && SII.is_timeseries_parameter(nw, sym)
        return Set{Union{Int, SII.ContinuousTimeseries}}([DEFAULT_PARA_TS_IDX])
    elseif !iszero(pdim(nw)) && SII.is_observed(nw, sym)
        return Set{Union{Int, SII.ContinuousTimeseries}}([SII.ContinuousTimeseries(), DEFAULT_PARA_TS_IDX])
    else
        return Set{Union{Int, SII.ContinuousTimeseries}}([SII.ContinuousTimeseries()])
    end
end
function SII.get_all_timeseries_indexes(nw::Network, sym::AbstractArray)
    return mapreduce(Base.Fix1(SII.get_all_timeseries_indexes, nw), union, sym;
        init = Set{Union{Int, SII.ContinuousTimeseries}}())
end

function SII.with_updated_parameter_timeseries_values(nw::Network, params, args::Pair...)
    @assert length(args) == 1 "Did not expect more than 1 timeseries here, please report issue."
    tsidx, p = args[1]
    @assert tsidx == DEFAULT_PARA_TS_IDX "Did not expect the passed timeseries to have other index then 1, please report issue."
    params .= p
end


function SciMLBase.create_parameter_timeseries_collection(nw::Network, p::AbstractVector, tspan)
    data = DiffEqArray(Vector{eltype(p)}[copy(p)], Float64[tspan[begin]])
    tsc = SII.ParameterTimeseriesCollection((data,), copy(p))
    return tsc
end

function SciMLBase.get_saveable_values(nw::Network, p::AbstractVector, timeseries_idx)
    @assert timeseries_idx == DEFAULT_PARA_TS_IDX # nothing else makes sense
    copy(p)
end
"""
    save_parameters!(integrator::SciMLBase.DEIntegrator)

Save the current parameter values in the integrator. Call this function inside callbacks
if the parameter values have changed. This will store a timeseries of said parameters in the
solution object, thus alowing us to recosntruct observables which depend on time-dependet variables.
"""
function save_parameters!(integrator::SciMLBase.DEIntegrator)
    SciMLBase.save_discretes!(integrator, DEFAULT_PARA_TS_IDX)
end

####
#### Observed indexing
####
# TODO: is_observed should probably also handle parameters?
function SII.is_observed(nw::Network, sni)
    if _hascolon(sni)
        SII.is_observed(nw, _resolve_colon(nw,sni))
    elseif SII.symbolic_type(sni) === SII.ArraySymbolic()
        # if has colon check if all are observed OR variables and return true
        # the observed function will handle the whole thing then
        all(s -> SII.is_variable(nw, s) || SII.is_observed(nw, s) || SII.is_parameter(nw, s), sni)
    elseif sni isa AbstractVector || sni isa Tuple
        any(SII.is_observed.(Ref(nw), sni))
    else
        _is_observed(nw, sni)
    end
end
_is_observed(nw::Network, _) = false
function _is_observed(nw::Network, sni::SymbolicIndex{<:CONCRETE_COMPIDX,Symbol})
    cf = getcomp(nw, sni)
    return sni.subidx ∈ obssym_all(cf)
end

function observed_symbols(nw::Network)
    syms = SymbolicIndex{Int,Symbol}[]
    for (ci, cf) in pairs(nw.im.vertexm)
        for s in obssym_all(cf)
            push!(syms, VIndex(ci, s))
        end
    end
    for (ci, cf) in pairs(nw.im.edgem)
        for s in obssym_all(cf)
            push!(syms, EIndex(ci, s))
        end
    end
    return syms
end

const U_TYPE = 1
const P_TYPE = 2
const OUT_TYPE = 3
const AGG_TYPE = 4
const OBS_TYPE = 5
function SII.observed(nw::Network, snis)
    if (snis isa AbstractVector || snis isa Tuple) && any(sni -> sni isa ObservableExpression, snis)
        throw(ArgumentError("Cannot mix normal symbolic indices with @obsex currently!"))
    end

    _is_normalized(snis) || return SII.observed(nw, _expand_and_collect(nw, snis))

    isscalar = snis isa SymbolicIndex
    _snis = isscalar ? (snis,) : snis

    # mapping i -> (U_TYPE, j) / (P_TYPE, j) / (OUT_TYPE, j) / (AGG_TYPE, j) / (OBS_TYPE, j in obsfuns)
    arraymapping = Vector{Tuple{Int, Int}}(undef, length(_snis))
    # vector of obs functions
    obsfuns = Vector{Function}()

    for (i, sni) in enumerate(_snis)
        if SII.is_variable(nw, sni)
            arraymapping[i] = (U_TYPE, SII.variable_index(nw, sni)::Int)
        elseif SII.is_parameter(nw, sni)
            arraymapping[i] = (P_TYPE, SII.parameter_index(nw, sni)::Int)
        else
            cf = getcomp(nw, sni)

            @argcheck sni.subidx isa Symbol "Observed musst be referenced by symbol, got $sni"
            if (idx=findfirst(isequal(sni.subidx), outsym_flat(cf))) != nothing # output
                _range = getcompoutrange(nw, sni)
                arraymapping[i] = (OUT_TYPE, _range[idx])
            elseif (idx=findfirst(isequal(sni.subidx), obssym(cf))) != nothing #found in observed
                _obsf = _get_observed_f(nw, cf, resolvecompidx(nw, sni))
                _newobsfun = let obsidx = idx # otherwise $idx is boxed everywhere in function
                    (u, outbuf, aggbuf, extbuf, p, t) -> _obsf(u, outbuf, aggbuf, extbuf, p, t)[obsidx]
                end
                push!(obsfuns, _newobsfun)
                arraymapping[i] = (OBS_TYPE, length(obsfuns))
            elseif hasinsym(cf) && sni.subidx ∈ insym_all(cf) # found in input
                if sni isa VIndex
                    idx = findfirst(isequal(sni.subidx), insym_all(cf))
                    arraymapping[i] = (AGG_TYPE, nw.im.v_aggr[resolvecompidx(nw, sni)][idx])
                elseif sni isa EIndex
                    edge = nw.im.edgevec[resolvecompidx(nw, sni)]
                    if (idx = findfirst(isequal(sni.subidx), insym(cf).src)) != nothing
                        arraymapping[i] = (OUT_TYPE, nw.im.v_out[edge.src][idx])
                    elseif (idx = findfirst(isequal(sni.subidx), insym(cf).dst)) != nothing
                        arraymapping[i] = (OUT_TYPE, nw.im.v_out[edge.dst][idx])
                    else
                        error()
                    end
                else
                    error()
                end
            else
                throw(ArgumentError("Cannot resolve observable $sni"))
            end
        end
    end
    needsbuf = any(m -> m[1] ∈ (OUT_TYPE, AGG_TYPE, OBS_TYPE), arraymapping)
    obsfunstup = Tuple(obsfuns) # make obsfuns concretely typed

    if isscalar
        (u, p, t) -> begin
            if needsbuf
                outbuf, aggbuf, extbuf = get_buffers(nw, u, p, t; initbufs=true)
            end
            type, idx = only(arraymapping)
            type == U_TYPE   && return u[idx]
            type == P_TYPE   && return p[idx]
            type == OUT_TYPE && return outbuf[idx]
            type == AGG_TYPE && return aggbuf[idx]
            type == OBS_TYPE && return only(obsfunstup)(u, outbuf, aggbuf, extbuf, p, t)::eltype(u)
        end
    else
        (u, p, t, out=similar(u, length(_snis))) -> begin
            if needsbuf
                outbuf, aggbuf, extbuf = get_buffers(nw, u, p, t; initbufs=true)
            end

            for (i, (type, idx)) in pairs(arraymapping)
                if type == U_TYPE
                    out[i] = u[idx]
                elseif type == P_TYPE
                    out[i] = p[idx]
                elseif type == OUT_TYPE
                    out[i] = outbuf[idx]
                elseif type == AGG_TYPE
                    out[i] = aggbuf[idx]
                elseif type == OBS_TYPE
                    out[i] = obsfunstup[idx](u, outbuf, aggbuf, extbuf, p, t)::eltype(u)
                end
            end
            return out
        end
    end
end
function _expand_and_collect(inpr, sni::SymbolicIndex)
    nw = extract_nw(inpr)
    if _hascolon(sni)
        collect(_resolve_colon(nw, sni))
    elseif SII.symbolic_type(sni) === SII.ArraySymbolic()
        collect(sni)
    else
        sni
    end
end
function _expand_and_collect(inpr, snis)
    nw = extract_nw(inpr)
    mapreduce(vcat, snis; init=Vector{eltype(snis)}()) do sni
        _expand_and_collect(nw, sni)
    end
end
function _is_normalized(snis)
    all(sni -> _is_normalized(sni), snis)
end
_is_normalized(snis::SymbolicIndex) = SII.symbolic_type(snis) === SII.ScalarSymbolic()

# function barrier for obsf
_get_observed_f(nw, cf, vidx) = _get_observed_f(nw.im, cf, vidx, obsf(cf))
function _get_observed_f(im::IndexManager, cf::VertexModel, vidx, _obsf::O) where {O}
    N = length(cf.obssym)
    ur   = im.v_data[vidx]
    aggr = im.v_aggr[vidx]
    extr = im.v_ext[vidx]
    pr   = im.v_para[vidx]
    retcache = DiffCache(Vector{Float64}(undef, N))
    _hasext = has_external_input(cf)

    (u, outbuf, aggbuf, extbuf, p, t) -> begin
        ret = PreallocationTools.get_tmp(retcache, first(u))
        ins = if _hasext
            (view(aggbuf, aggr), view(extbuf, extr))
        else
            (view(aggbuf, aggr), )
        end
        _obsf(ret, view(u, ur), ins..., view(p, pr), t)
        ret
    end
end
function _get_observed_f(im::IndexManager, cf::EdgeModel, eidx, _obsf::O) where {O}
    N = length(cf.obssym)
    ur    = im.e_data[eidx]
    esrcr = im.v_out[im.edgevec[eidx].src]
    edstr = im.v_out[im.edgevec[eidx].dst]
    extr  = im.e_ext[eidx]
    pr    = im.e_para[eidx]
    ret = Vector{Float64}(undef, N)
    _hasext = has_external_input(cf)

    (u, outbuf, aggbuf, extbuf, p, t) -> begin
        ins = if _hasext
            (view(outbuf, esrcr), view(outbuf, edstr), view(extbuf, extr))
        else
            (view(outbuf, esrcr), view(outbuf, edstr))
        end
        _obsf(ret, view(u, ur), ins..., view(p, pr), t)
        ret
    end
end


####
#### Default values
####
function SII.default_values(nw::Network)
    aliased_changed(nw; warn=true)
    defs = Dict{SymbolicIndex{Int,Symbol},Float64}()
    for (ci, cf) in pairs(nw.im.vertexm)
        for s in psym(cf)
            has_default_or_init(cf, s) || continue
            defs[VIndex(ci, s)] = get_default_or_init(cf, s)
        end
        for s in sym(cf)
            has_default_or_init(cf, s) || continue
            defs[VIndex(ci, s)] = get_default_or_init(cf, s)
        end
    end
    for (ci, cf) in pairs(nw.im.edgem)
        for s in psym(cf)
            has_default_or_init(cf, s) || continue
            defs[EIndex(ci, s)] = get_default_or_init(cf, s)
        end
        for s in sym(cf)
            has_default_or_init(cf, s) || continue
            defs[EIndex(ci, s)] = get_default_or_init(cf,s)
        end
    end
    return defs
end
