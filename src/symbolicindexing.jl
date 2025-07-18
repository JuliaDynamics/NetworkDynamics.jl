"""
    VIndex{C,S} <: SymbolicStateIndex{C,S}
    idx = VIndex(comp, sub)

A symbolic index for a vertex state variable.
- `comp`: the component index, either int, symbol or a collection
- `sub`: the subindex, either int, symbol or a collection of those.

```
VIndex(1, :P)      # vertex 1, variable :P
VIndex(1:5, 1)     # first state of vertices 1 to 5
VIndex(7, (:x,:y)) # states :x and :y of vertex 7
VIndex(2)          # references the second vertex model
VIndex(:a)         # references vertex with unique name :a
```

Can be used to index into objects supporting the `SymbolicIndexingInterface`,
e.g. [`NWState`](@ref), [`NWParameter`](@ref) or `ODESolution`.

See also: [`EIndex`](@ref), [`VPIndex`](@ref), [`EPIndex`](@ref)
"""
struct VIndex{C,S} <: SymbolicStateIndex{C,S}
    compidx::C
    subidx::S
end
VIndex(ci::Union{Symbol,Int}) = VIndex(ci, nothing)
"""
    EIndex{C,S} <: SymbolicStateIndex{C,S}
    idx = EIndex(comp, sub)

A symbolic index for an edge state variable.
- `comp`: the component index, either int, symbol, pair or a collection
- `sub`: the subindex, either int, symbol or a collection of those.

```
EIndex(1, :P)      # edge 1, variable :P
EIndex(1:5, 1)     # first state of edges 1 to 5
EIndex(7, (:x,:y)) # states :x and :y of edge 7
EIndex(2)          # references the second edge model
EIndex(1=>2)       # references edge from v1 to v2
EIndex(:a=>:b)     # references edge from (uniquely named) vertex :a to :b
```

Can be used to index into objects supporting the `SymbolicIndexingInterface`,
e.g. [`NWState`](@ref), [`NWParameter`](@ref) or `ODESolution`.

See also: [`VIndex`](@ref), [`VPIndex`](@ref), [`EPIndex`](@ref)
"""
struct EIndex{C,S} <: SymbolicStateIndex{C,S}
    compidx::C
    subidx::S
end
EIndex(ci::Union{Symbol,Int}) = EIndex(ci, nothing)
"""
    VPIndex{C,S} <: SymbolicStateIndex{C,S}
    idx = VPIndex(comp, sub)

A symbolic index into the parameter a vertex:
- `comp`: the component index, either int, symbol or a collection
- `sub`: the subindex, either int, symbol or a collection of those.

Can be used to index into objects supporting the `SymbolicIndexingInterface`,
e.g. [`NWParameter`](@ref) or `ODEProblem`.

See also: [`EPIndex`](@ref), [`VIndex`](@ref), [`EIndex`](@ref)
"""
struct VPIndex{C,S} <: SymbolicParameterIndex{C,S}
    compidx::C
    subidx::S
end
"""
    EPIndex{C,S} <: SymbolicStateIndex{C,S}
    idx = VEIndex(comp, sub)

A symbolic index into the parameter of an edge:
- `comp`: the component index, either int, symbol, pair or a collection
- `sub`: the subindex, either int, symbol or a collection of those.

Can be used to index into objects supporting the `SymbolicIndexingInterface`,
e.g. [`NWParameter`](@ref) or `ODEProblem`.

See also: [`VPIndex`](@ref), [`VIndex`](@ref), [`EIndex`](@ref)
"""
struct EPIndex{C,S} <: SymbolicParameterIndex{C,S}
    compidx::C
    subidx::S
end
const SymbolicEdgeIndex{C,S} = Union{EIndex{C,S}, EPIndex{C,S}}
const SymbolicVertexIndex{C,S} = Union{VIndex{C,S}, VPIndex{C,S}}

idxtype(s::VIndex) = VIndex
idxtype(s::EIndex) = EIndex


SII.symbolic_type(::Type{<:SymbolicIndex{<:Union{<:Pair,Symbol,Int},<:Union{Symbol,Int}}}) = SII.ScalarSymbolic()
SII.symbolic_type(::Type{<:SymbolicIndex}) = SII.ArraySymbolic()

SII.hasname(::SymbolicIndex) = false
SII.hasname(::SymbolicIndex{<:Union{<:Pair,Symbol,Int},<:Union{Symbol,Int}}) = true
function SII.getname(x::SymbolicVertexIndex)
    prefix = x.compidx isa Int ? :v : Symbol()
    Symbol(prefix, Symbol(x.compidx), :₊, Symbol(x.subidx))
end
function SII.getname(x::SymbolicEdgeIndex)
    if x.compidx isa Pair
        src, dst = x.compidx
        _src = src isa Int ? Symbol(:v, src) : Symbol(src)
        _dst = dst isa Int ? Symbol(:v, dst) : Symbol(dst)
        Symbol(_src, "ₜₒ", _dst, :₊, Symbol(x.subidx))
    else
        prefix = x.compidx isa Int ? :e : Symbol()
        Symbol(prefix, Symbol(x.compidx), :₊, Symbol(x.subidx))
    end
end

resolvecompidx(nw::Network, sni) = resolvecompidx(nw.im, sni)
resolvecompidx(::IndexManager, sni::SymbolicIndex{Int}) = sni.compidx
function resolvecompidx(im::IndexManager, sni::SymbolicIndex{Symbol})
    dict = sni isa SymbolicVertexIndex ? im.unique_vnames : im.unique_enames
    if haskey(dict, sni.compidx)
        return dict[sni.compidx]
    else
        throw(ArgumentError("Could not resolve component index for $sni, the name might not be unique?"))
    end
end
function resolvecompidx(im::IndexManager, sni::SymbolicEdgeIndex{<:Pair})
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
getcomp(im::IndexManager, sni::SymbolicEdgeIndex) = im.edgem[resolvecompidx(im, sni)]
getcomp(im::IndexManager, sni::SymbolicVertexIndex) = im.vertexm[resolvecompidx(im, sni)]

getcomprange(nw::Network, sni) = getcomprange(nw.im, sni)
getcomprange(im::IndexManager, sni::VIndex{<:Union{Symbol,Int}}) = im.v_data[resolvecompidx(im, sni)]
getcomprange(im::IndexManager, sni::EIndex{<:Union{<:Pair,Symbol,Int}}) = im.e_data[resolvecompidx(im, sni)]

getcompoutrange(nw::Network, sni) = getcompoutrange(nw.im, sni)
getcompoutrange(im::IndexManager, sni::VIndex{<:Union{Symbol,Int}}) = im.v_out[resolvecompidx(im, sni)]
getcompoutrange(im::IndexManager, sni::EIndex{<:Union{<:Pair,Symbol,Int}}) = flatrange(im.e_out[resolvecompidx(im, sni)])

getcompprange(nw::Network, sni::SymbolicVertexIndex{<:Union{Symbol,Int}}) = nw.im.v_para[resolvecompidx(nw, sni)]
getcompprange(nw::Network, sni::SymbolicEdgeIndex{<:Union{<:Pair,Symbol,Int}}) = nw.im.e_para[resolvecompidx(nw, sni)]

subsym_has_idx(sym::Symbol, syms) = sym ∈ syms
subsym_has_idx(idx::Int, syms) = 1 ≤ idx ≤ length(syms)
subsym_to_idx(sym::Symbol, syms) = findfirst(isequal(sym), syms)
subsym_to_idx(idx::Int, _) = idx

####
#### Iterator/Broadcast interface for ArraySymbolic types
####
# TODO: not broadcasting over idx with colon is weird
Base.broadcastable(si::SymbolicIndex{<:Union{Int,Symbol,<:Pair,Colon},<:Union{Int,Symbol,Colon}}) = Ref(si)

const _IterableComponent = SymbolicIndex{<:Union{AbstractVector,Tuple},<:Union{Int,Symbol}}
Base.length(si::_IterableComponent) = length(si.compidx)
Base.size(si::_IterableComponent) = (length(si),)
Base.IteratorSize(si::_IterableComponent) = Base.HasShape{1}()
Base.broadcastable(si::_IterableComponent) = si
Base.ndims(::Type{<:_IterableComponent}) = 1
Base.axes(si::_IterableComponent) = axes(si.compidx)
Base.getindex(si::_IterableComponent, i) = _baseT(si)(si.compidx[i], si.subidx)
function Base.eltype(si::_IterableComponent)
    if isconcretetype(eltype(si.compidx))
        _baseT(si){eltype(si.compidx),typeof(si.subidx)}
    else
        Any
    end
end
function Base.iterate(si::_IterableComponent, state=nothing)
    it = isnothing(state) ? iterate(si.compidx) : iterate(si.compidx, state)
    isnothing(it) && return nothing
    _similar(si, it[1], si.subidx), it[2]
end

const _IterableSubcomponent = SymbolicIndex{<:Union{<:Pair,Symbol,Int},<:Union{AbstractVector,Tuple}}
Base.length(si::_IterableSubcomponent) = length(si.subidx)
Base.size(si::_IterableSubcomponent) = (length(si),)
Base.IteratorSize(si::_IterableSubcomponent) = Base.HasShape{1}()
Base.broadcastable(si::_IterableSubcomponent) = si
Base.ndims(::Type{<:_IterableSubcomponent}) = 1
Base.axes(si::_IterableSubcomponent) = axes(si.subidx)
Base.getindex(si::_IterableSubcomponent, i) = _baseT(si)(si.compidx, si.subidx[i])
function Base.eltype(si::_IterableSubcomponent)
    if isconcretetype(eltype(si.subidx))
        _baseT(si){eltype(si.compidx),eltype(si.subidx)}
    else
        Any
    end
end
function Base.iterate(si::_IterableSubcomponent, state=nothing)
    it = isnothing(state) ? iterate(si.subidx) : iterate(si.subidx, state)
    isnothing(it) && return nothing
    _similar(si, si.compidx, it[1]), it[2]
end
_similar(si::SymbolicIndex, c, s) = _baseT(si)(c,s)
_baseT(::VIndex)  = VIndex
_baseT(::EIndex)  = EIndex
_baseT(::VPIndex) = VPIndex
_baseT(::EPIndex) = EPIndex

_hascolon(i) = false
_hascolon(::SymbolicIndex{C,S}) where {C,S} = C === Colon || S === Colon
_resolve_colon(nw::Network, idx) = idx
_resolve_colon(nw::Network, sni::VIndex{Colon}) = VIndex(1:nv(nw), sni.subidx)
_resolve_colon(nw::Network, sni::EIndex{Colon}) = EIndex(1:ne(nw), sni.subidx)
_resolve_colon(nw::Network, sni::VPIndex{Colon}) = VPIndex(1:nv(nw), sni.subidx)
_resolve_colon(nw::Network, sni::EPIndex{Colon}) = EPIndex(1:ne(nw), sni.subidx)
_resolve_colon(nw::Network, sni::VIndex{<:Union{Symbol,Int},Colon}) = VIndex{Int, UnitRange{Int}}(sni.compidx, 1:dim(getcomp(nw,sni)))
_resolve_colon(nw::Network, sni::EIndex{<:Union{<:Pair,Symbol,Int},Colon}) = EIndex{Int, UnitRange{Int}}(sni.compidx, 1:dim(getcomp(nw,sni)))
_resolve_colon(nw::Network, sni::VPIndex{<:Union{Symbol,Int},Colon}) = VPIndex{Int, UnitRange{Int}}(sni.compidx, 1:pdim(getcomp(nw,sni)))
_resolve_colon(nw::Network, sni::EPIndex{<:Union{<:Pair,Symbol,Int},Colon}) = EPIndex{Int, UnitRange{Int}}(sni.compidx, 1:pdim(getcomp(nw,sni)))


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
const POTENTIAL_SCALAR_SIDX = Union{SymbolicStateIndex{<:Union{<:Pair,Symbol,Int},<:Union{Int,Symbol}}}
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
    syms = Vector{SymbolicStateIndex{Int,Symbol}}(undef, dim(nw))
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
const POTENTIAL_SCALAR_PIDX = Union{
    SymbolicParameterIndex{<:Union{<:Pair,Symbol,Int},<:Union{Int,Symbol}},
    SymbolicIndex{<:Union{<:Pair,Symbol,Int},Symbol}
}
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
    syms = Vector{SymbolicParameterIndex{Int,Symbol}}(undef, pdim(nw))
    for (ci, cf) in pairs(nw.im.vertexm)
        syms[nw.im.v_para[ci]] .= VPIndex.(ci, psym(cf))
    end
    for (ci, cf) in pairs(nw.im.edgem)
        syms[nw.im.e_para[ci]] .= EPIndex.(ci, psym(cf))
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
function SII.is_observed(nw::Network, sni)
    if _hascolon(sni)
        SII.is_observed(nw, _resolve_colon(nw,sni))
    elseif SII.symbolic_type(sni) === SII.ArraySymbolic()
        # if has colon check if all are observed OR variables and return true
        # the observed function will handle the whole thing then
        all(s -> SII.is_variable(nw, s) || SII.is_observed(nw, s), sni)
    elseif sni isa AbstractVector
        any(SII.is_observed.(Ref(nw), sni))
    else
        _is_observed(nw, sni)
    end
end
_is_observed(nw::Network, _) = false
function _is_observed(nw::Network, sni::SymbolicStateIndex{<:Union{<:Pair,Symbol,Int},Symbol})
    cf = getcomp(nw, sni)
    return sni.subidx ∈ obssym_all(cf)
end

function observed_symbols(nw::Network)
    syms = SymbolicStateIndex{Int,Symbol}[]
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
                if sni isa SymbolicVertexIndex
                    idx = findfirst(isequal(sni.subidx), insym_all(cf))
                    arraymapping[i] = (AGG_TYPE, nw.im.v_aggr[resolvecompidx(nw, sni)][idx])
                elseif sni isa SymbolicEdgeIndex
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
    mapreduce(vcat, snis) do sni
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
    extr = im.v_out[vidx]
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
    extr  = im.e_out[eidx]
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
            defs[VPIndex(ci, s)] = get_default_or_init(cf, s)
        end
        for s in sym(cf)
            has_default_or_init(cf, s) || continue
            defs[VIndex(ci, s)] = get_default_or_init(cf, s)
        end
    end
    for (ci, cf) in pairs(nw.im.edgem)
        for s in psym(cf)
            has_default_or_init(cf, s) || continue
            defs[EPIndex(ci, s)] = get_default_or_init(cf, s)
        end
        for s in sym(cf)
            has_default_or_init(cf, s) || continue
            defs[EIndex(ci, s)] = get_default_or_init(cf,s)
        end
    end
    return defs
end


####
#### NWParameter and NWState objects as value provider
####

"""
    NWParameter(nw_or_nw_wraper, pflat)

Indexable wrapper for flat parameter array `pflat`. Needs Network or wrapper of
Network, e.g. `ODEProblem`.

```
p = NWParameter(nw)
p.v[idx, :sym] # get parameter :sym of vertex idx
p.e[idx, :sym] # get parameter :sym of edge idx
p[s::Union{VPIndex, EPIndex}] # get parameter for specific index
```

Get flat array representation using [`pflat`](@ref). The order of parameters in the flat
representation corresponds to the order given by [`parameter_symbols`](@ref).
"""
struct NWParameter{P,NW<:Network}
    nw::NW
    pflat::P
    function NWParameter(thing, pflat)
        nw = extract_nw(thing)
        new{typeof(pflat),typeof(nw)}(nw, pflat)
    end
end

Base.eltype(p::NWParameter) = eltype(p.pflat)
Base.length(s::NWParameter) = length(s.pflat)

"""
    NWParameter(nw_or_nw_wraper;
                ptype=Vector{Float64}, pfill=filltype(ptype), default=true)

Creates "empty" `NWParameter` object for the Network/Wrapper `nw` with flat type `ptype`.
The array will be prefilled with `pfill` (defaults to NaN).

If `default=true` the default parameter values attached to the network components will be loaded.
"""
function NWParameter(thing; ptype=Vector{Float64}, pfill=filltype(ptype), default=true)
    nw = extract_nw(thing)
    pflat = _init_flat(ptype, pdim(nw), pfill)
    p = NWParameter(nw, pflat)
    default || return p
    for (k, v) in SII.default_values(nw)
        k isa Union{VPIndex, EPIndex} || continue
        p[k] = v
    end
    return p
end

"""
    NWParameter(p::NWParameter; ptype=typeof(p.pflat))

Create `NWParameter` based on other parameter object, just convert type.
"""
function NWParameter(p::NWParameter; ptype=typeof(p.pflat))
    NWParameter(p.nw, _convertorcopy(ptype,pflat(p)))
end

"""
    NWParameter(int::SciMLBase.DEIntegrator)

Create `NWParameter` object from `integrator`.
"""
NWParameter(int::SciMLBase.DEIntegrator) = NWParameter(int, int.p)


"""
    NWState(nw_or_nw_wrapper, uflat, [pflat], [t])

Indexable wrapper for flat state & parameter array. Needs Network or wrapper of
Network, e.g. `ODEProblem`.

```
s = NWState(nw)
s.v[idx, :sym] # get state :sym of vertex idx
s.e[idx, :sym] # get state :sym of edge idx
s.p.v[idx, :sym] # get parameter :sym of vertex idx
s.p.e[idx, :sym] # get parameter :sym of edge idx
s[s::Union{VIndex, EIndex, EPIndex, VPIndex}] # get parameter for specific index
```

Get flat array representation using [`uflat`](@ref) and [`pflat`](@ref). The order of states in the flat
representation corresponds to the order given by [`variable_symbols`](@ref), and the order of parameters
corresponds to [`parameter_symbols`](@ref).
"""
struct NWState{U,P,T,NW<:Network}
    nw::NW
    uflat::U
    p::P
    t::T
    function NWState(thing, uflat::AbstractVector, p=nothing, t=nothing)
        nw = extract_nw(thing)
        _p = p isa Union{NWParameter,Nothing} ? p : NWParameter(nw, p)
        s = new{typeof(uflat),typeof(_p),typeof(t),typeof(nw)}(nw,uflat,_p,t)
        @argcheck isnothing(p) || s.nw === s.p.nw
        return s
    end
end


"""
    NWState(nw_or_nw_wrapper;
            utype=Vector{Float64}, ufill=filltype(utype),
            ptype=Vector{Float64}, pfill=filltype(ptype), default=true)

Creates "empty" `NWState` object for the Network/Wrapper `nw` with flat types
`utype` & `ptype`. The arrays will be prefilled with `ufill` and `pfill`
respectively (defaults to NaN).

If `default=true` the default state & parameter values attached to the network
components will be loaded.
"""
function NWState(thing;
                 utype=Vector{Float64}, ufill=filltype(utype),
                 ptype=Vector{Float64}, pfill=filltype(ptype),
                 default=true)
    nw = extract_nw(thing)
    t = nothing
    uflat = _init_flat(utype, dim(nw), ufill)
    p = NWParameter(nw; ptype, pfill, default=false)
    s = NWState(nw,uflat,p,t)
    default || return s
    for (k, v) in SII.default_values(nw)
        s[k] = v
    end
    return s
end

"""
    NWState(p::NWState; utype=typeof(uflat(s)), ptype=typeof(pflat(s)))

Create `NWState` based on other state object, just convert types.
"""
function NWState(s::NWState;
                 utype=typeof(uflat(s)),
                 ptype=typeof(pflat(s)))
    p = NWParameter(s.p; ptype)
    NWState(s.nw, _convertorcopy(utype,uflat(s)), p, s.t)
end

"""
    NWState(p::NWParameter; utype=Vector{Float64}, ufill=filltype(utype), default=true)

Create `NWState` based on existing `NWParameter` object.
"""
function NWState(p::NWParameter; utype=Vector{Float64}, ufill=filltype(utype), default=true)
    t = nothing
    nw = p.nw
    uflat = _init_flat(utype, dim(nw), ufill)
    s = NWState(nw,uflat,p,t)
    default || return s
    for (k, v) in SII.default_values(nw)
        k isa SymbolicParameterIndex && continue
        s[k] = v
    end
    return s
end

"""
    NWState(int::SciMLBase.DEIntegrator)

Create `NWState` object from `integrator`.
"""
NWState(int::SciMLBase.DEIntegrator) = NWState(int, int.u, int.p, int.t)

"""
    NWState(sol::SciMLBase.AbstractODESolution, t)

Create `NWState` object from solution object `sol` for timepoint `t`.
"""
function NWState(sol::SciMLBase.AbstractODESolution, t::Number)
    u = sol(t)
    discs = RecursiveArrayTools.get_discretes(sol)
    para_ts = discs[DEFAULT_PARA_TS_IDX]

    tidx = findfirst(_t -> _t > t, para_ts.t)
    p = if isnothing(tidx)
        para_ts.u[end]
    else
        para_ts.u[tidx-1]
    end

    NWState(sol, u, p, t)
end

# init flat array of type T with length N. Init with nothing if possible, else with zeros
function _init_flat(T, N, fill)
    vec = T(undef, N)
    fill!(vec, fill)
end

"""
    filltype(T)

Return a value which will be used to fill an abstract array of type `T`.

- `nothing` if eltype(T) allows nothing
- `missing` if eltype(T) allows missing
- `NaN` if eltype(T) is a `AbstractFloat
- `zero(eltype(T))` else.
"""
function filltype(T::Type{<:AbstractVector})
    elT = eltype(T)
    if Nothing <: elT
        nothing
    elseif Missing <: elT
        missing
    elseif elT <: AbstractFloat
        NaN
    else
        zero(elT)
    end
end

_convertorcopy(::Type{T}, x::T) where {T} = copy(x)
_convertorcopy(T, x) = convert(T, x)

Base.eltype(s::NWState) = eltype(s.uflat)
Base.length(s::NWState) = length(s.uflat)
"""
    uflat(s::NWState)

Retrieve the wrapped flat array representation of the state. The order of states in this
flat representation corresponds exactly to the order given by [`variable_symbols`](@ref).
"""
uflat(s::NWState) = s.uflat
"""
    pflat(p::NWParameter)
    pflat(s::NWState)

Retrieve the wrapped flat array representation of the parameters. The order of parameters in this
flat representation corresponds exactly to the order given by [`parameter_symbols`](@ref).
"""
pflat(s::NWState) = pflat(s.p)
pflat(p::NWParameter) = p.pflat

Base.broadcastable(s::NWState) = Ref(s)

SII.symbolic_container(s::NWState) = s.nw
SII.symbolic_container(s::NWParameter) = s.nw
SII.state_values(s::NWState) = s.uflat
SII.state_values(s::NWParameter) = error("Parameter type does not hold State values.")
SII.parameter_values(s::NWState) = s.p isa NWParameter ? SII.parameter_values(s.p) : s.p
SII.parameter_values(p::NWParameter) = p.pflat
SII.current_time(s::NWState) = s.t
SII.current_time(s::NWParameter) = error("Parameter type does not holde time value.")

# NWParameter: getindex
Base.getindex(p::NWParameter, ::Colon) = pflat(p)
# HACK: _expand_and_collect to workaround https://github.com/SciML/SymbolicIndexingInterface.jl/issues/94
Base.getindex(p::NWParameter, idx) = SII.getp(p, _expand_and_collect(p, _paraindex(idx)))(p)

# NWParameter: setindex!
function Base.setindex!(p::NWParameter, val, idx)
    setter = SII.setp(p, _paraindex(idx))
    _chk_dimensions(setter, val)
    setter(p, val)
end

# Converts a given index (collection) to a (collection) of parameter index if it is a state index.
_paraindex(idxs) = all(i -> i isa SymbolicParameterIndex, idxs) ? idxs : _paraindex.(idxs)
_paraindex(idx::VIndex) = VPIndex(idx.compidx, idx.subidx)
_paraindex(idx::EIndex) = EPIndex(idx.compidx, idx.subidx)
_paraindex(idx::SymbolicParameterIndex) = idx

# NWState: getindex
Base.getindex(s::NWState, ::Colon) = uflat(s)
# HACK: _expand_and_collect to workaround https://github.com/SciML/SymbolicIndexingInterface.jl/issues/94
Base.getindex(s::NWState, idx::SymbolicParameterIndex) = SII.getp(s, _expand_and_collect(s, idx))(s)
Base.getindex(s::NWState, idx::SymbolicStateIndex) = SII.getu(s, idx)(s)
function Base.getindex(s::NWState, idxs)
    if all(i -> i isa SymbolicParameterIndex, idxs)
        # HACK: _expand_and_collect to workaround https://github.com/SciML/SymbolicIndexingInterface.jl/issues/94
        SII.getp(s, _expand_and_collect(s, idxs))(s)
    else
        SII.getu(s, idxs)(s)
    end
end

# NWState: setindex!
function Base.setindex!(s::NWState, val, idx::SymbolicIndex)
    setter = if idx isa SymbolicParameterIndex
        SII.setp(s, idx)
    else
        SII.setu(s, idx)
    end
    _chk_dimensions(setter, val)
    setter(s, val)
end
function Base.setindex!(s::NWState, val, idxs)
    setter = if all(i -> i isa SymbolicParameterIndex, idxs)
        SII.setp(s, idxs)
    else
        SII.setu(s, idxs)
    end
    _chk_dimensions(setter, val)
    setter(s, val)
end

function _chk_dimensions(::Union{SII.SetParameterIndex,SII.SetStateIndex}, val)
    length(val) == 1 || throw(DimensionMismatch("Cannot set multiple values to single index."))
end
_chk_dimensions(s::SII.ParameterHookWrapper, val) = _chk_dimensions(s.setter, val)
function _chk_dimensions(ms::SII.MultipleSetters, val)
    if size(ms.setters) != size(val)
        throw(DimensionMismatch("Cannot set variables of size $(size(ms.setters)) to values of size $(size(val))."))
    end
end


####
#### Indexing proxys
####
abstract type IndexingProxy end
struct VProxy{S} <: IndexingProxy
    s::S
end
Base.getindex(p::VProxy, comp, state) = getindex(p.s, VIndex(comp, state))
Base.getindex(p::VProxy, ::Colon, state) = getindex(p, 1:nv(extract_nw(p)), state)
Base.setindex!(p::VProxy, val, comp, state) = setindex!(p.s, val, VIndex(comp, state))

struct EProxy{S} <: IndexingProxy
    s::S
end
Base.getindex(p::EProxy, comp, state) = getindex(p.s, EIndex(comp, state))
Base.getindex(p::EProxy, ::Colon, state) = getindex(p, 1:ne(extract_nw(p)), state)
Base.setindex!(p::EProxy, val, comp, state) = setindex!(p.s, val, EIndex(comp, state))

function Base.getproperty(s::Union{NWParameter, NWState}, sym::Symbol)
    if sym === :v
        return VProxy(s)
    elseif sym === :e
        return EProxy(s)
    else
        return getfield(s, sym)
    end
end


####
#### enable broadcasted setindex
#### https://discourse.julialang.org/t/broadcasting-setindex-is-a-noobtrap/94700
####
Base.dotview(s::Union{NWParameter, NWState, VProxy, EProxy}, idxs...) = view(s, idxs...)
Base.view(p::VProxy, comp, state) = view(p.s, VIndex(comp, state))
Base.view(p::EProxy, comp, state) = view(p.s, EIndex(comp, state))

# NWParameter: view
Base.view(s::NWParameter, ::Colon) = s.pflat
function Base.view(p::NWParameter, idx::SymbolicIndex)
    _idx = _paraindex(idx)
    if !SII.is_parameter(p, _idx)
        throw(ArgumentError("Index $idx is not a valid parameter index."))
    end
    view(p.pflat, SII.parameter_index(p, _idx))
end
function Base.view(p::NWParameter, idxs)
    _idxs = _paraindex(idxs)
    if !(all(i -> SII.is_parameter(p, i), _idxs))
        throw(ArgumentError("Index $idxs is not a valid parameter index collection."))
    end
    view(p.pflat, map(i -> SII.parameter_index(p, i), _idxs))
end

# NWState: view
Base.view(s::NWState, ::Colon) = s.uflat
Base.view(s::NWState, idx::SymbolicParameterIndex) = view(s.p, idx)
function Base.view(s::NWState, idx::SymbolicStateIndex)
    if !SII.is_variable(s, idx)
        throw(ArgumentError("Index $idx is not a valid state index."))
    end
    view(uflat(s), SII.variable_index(s, idx))
end
function Base.view(s::NWState, idxs)
    if all(i -> SII.is_parameter(s, i), idxs)
        _viewidx =  map(i -> SII.parameter_index(s, i), idxs)
        return view(pflat(s), _viewidx)
    elseif all(i -> SII.is_variable(s, i), idxs)
        _viewidx =  map(i -> SII.variable_index(s, i), idxs)
        return view(uflat(s), _viewidx)
    else
        throw(ArgumentError("Index $idx is neither a valid parameter nor state index collection."))
    end
end


# TODO: vidx(nw, :, :u) has different semantics from s.v[:, :u] (one searches?)
"""
    vidxs([inpr], components=:, variables=:) :: Vector{VIndex}

Generate vector of symbolic indexes for vertices.

- `inpr`: Only needed for name matching or `:` access. Can be Network, sol, prob, ...
- `components`: Number/Vector, `:`, `Symbol` (name matches), `String`/`Regex` (name contains)
- `variables`: Symbol/Number/Vector, `:`, `String`/`Regex` (all sym containing)

Examples:

    vidxs(nw)                 # all vertex state indices
    vidxs(1:2, :u)            # [VIndex(1, :u), VIndex(2, :u)]
    vidxs(nw, :, [:u, :v])    # [VIndex(i, :u), VIndex(i, :v) for i in 1:nv(nw)]
    vidxs(nw, "ODEVertex", :) # all symbols of all vertices with name containing "ODEVertex"
"""
vidxs(args...) = _idxs(VIndex, args...)
"""
    vpidxs([inpr], components=:, variables=:) :: Vector{VPIndex}

Generate vector of symbolic indexes for parameters. See [`vidxs`](@ref) for more information.
"""
vpidxs(args...) = _idxs(VPIndex, args...)
"""
    vidxs([inpr], components=:, variables=:) :: Vector{EIndex}

Generate vector of symbolic indexes for edges.

- `inpr`: Only needed for name matching or `:` access. Can be Network, sol, prob, ...
- `components`: Number/Vector, `:`, `Symbol` (name matches), `String`/`Regex` (name contains)
- `variables`: Symbol/Number/Vector, `:`, `String`/`Regex` (all sym containing)

Examples:

    eidxs(nw)                # all edge state indices
    eidxs(1:2, :u)           # [EIndex(1, :u), EIndex(2, :u)]
    eidxs(nw, :, [:u, :v])   # [EIndex(i, :u), EIndex(i, :v) for i in 1:ne(nw)]
    eidxs(nw, "FlowEdge", :) # all symbols of all edges with name containing "FlowEdge"
"""
eidxs(args...) = _idxs(EIndex, args...)
"""
    epidxs([inpr], components=:, variables=:) :: Vector{EPIndex}

Generate vector of symbolic indexes for parameters. See [`eidxs`](@ref) for more information.
"""
epidxs(args...) = _idxs(EPIndex, args...)

_idxs(IT, cidxs::Union{Number, AbstractVector}, sidxs) = _idxs(IT, nothing, cidxs, sidxs)
function _idxs(IT, inpr, cidxs=Colon(), sidxs=Colon())
    res = IT[]
    for ci in _make_cidx_iterable(IT, inpr, cidxs)
        for si in _make_sidx_iterable(IT, inpr, ci, sidxs)
            push!(res, IT(ci, si))
        end
    end
    res
end

_make_iterabel(idxs) = idxs
_make_iterabel(idx::Symbol) = Ref(idx)

_make_cidx_iterable(_, _, idx) = _make_iterabel(idx)
_make_cidx_iterable(::Type{<:SymbolicVertexIndex}, inpr, ::Colon) = 1:nv(extract_nw(inpr))
_make_cidx_iterable(::Type{<:SymbolicEdgeIndex}, inpr, ::Colon) = 1:ne(extract_nw(inpr))
function _make_cidx_iterable(IT, inpr, s::Symbol)
    names = getproperty.(_get_components(IT, inpr), :name)
    findall(isequal(s), names)
end
function _make_cidx_iterable(IT, inpr, s::Union{AbstractString,AbstractPattern})
    names = getproperty.(_get_components(IT, inpr), :name)
    findall(sym -> contains(string(sym), s), names)
end

_make_sidx_iterable(IT, inpr, cidx, idx) = _make_iterabel(idx)
function _make_sidx_iterable(IT::Type{<:SymbolicStateIndex}, inpr, cidx, ::Colon)
    comp = _get_components(IT, inpr)[cidx]
    Iterators.flatten((sym(comp), obssym_all(comp)))
end
function _make_sidx_iterable(IT::Type{<:SymbolicParameterIndex}, inpr, cidx, ::Colon)
    psym(_get_components(IT, inpr)[cidx])
end
function _make_sidx_iterable(IT::Type{<:SymbolicStateIndex}, inpr, cidx, s::Union{AbstractString,AbstractPattern})
    syms = _make_sidx_iterable(IT, inpr, cidx, :) # get all possible
    Iterators.filter(sym -> contains(string(sym), s), syms)
end
function _make_sidx_iterable(IT::Type{<:SymbolicParameterIndex}, inpr, cidx, s::Union{AbstractString,AbstractPattern})
    syms = _make_sidx_iterable(IT, inpr, cidx, :) # get all possible
    filter(sym -> contains(string(sym), s), syms)
end

_get_components(::Type{<:SymbolicVertexIndex}, inpr) = extract_nw(inpr).im.vertexm
_get_components(::Type{<:SymbolicEdgeIndex}, inpr) = extract_nw(inpr).im.edgem


"""
    extract_nw(thing)

Try to extract the `Network` object from thing.

Thing can by many things, e.g. `ODEProblem`, `ODESolution`, `Integrator`, `NWState`, `NWParameter`, ...
"""
function extract_nw(inpr)
    sc = SII.symbolic_container(inpr)
    if sc === inpr
       throw(ArgumentError("Cannot extract Network from $(typeof(sc))"))
    end
    extract_nw(sc)
end
extract_nw(nw::Network) = nw
extract_nw(sol::SciMLBase.AbstractSolution) = extract_nw(sol.prob)
extract_nw(prob::SciMLBase.ODEProblem) = extract_nw(prob.f)
extract_nw(f::SciMLBase.ODEFunction) = extract_nw(f.sys)
extract_nw(int::SciMLBase.DEIntegrator) = extract_nw(int.f)
extract_nw(s::NWState) = s.nw
extract_nw(p::NWParameter) = p.nw
extract_nw(p::IndexingProxy) = extract_nw(p.s)
function extract_nw(::Nothing)
    throw(ArgumentError("Needs system context to generate matching indices. Pass Network, sol, prob, ..."))
end

####
#### Observable Expressions
####
struct ObservableExpression{VT,F,Ex,N}
    inputs::Vector{VT}
    f::F
    ex::Ex
    name::N
end
function Base.show(io::IO, mime::MIME"text/plain", obsex::ObservableExpression)
    print(io, "ObservableExpression(")
    isnothing(obsex.name) || print(io, obsex.name, " = ")
    show(io, mime, obsex.ex)
    print(io, ")")
end

"""
    @obsex([name =] expression)

Define observable expressions, which are simple combinations of knonw
states/parameters/observables. `@obsex(...)` returns an `ObservableExpression`
which can be used as an symbolic index. This is mainly intended for quick
plotting or export of common "derived" variables, such as the argument of a
2-component complex state. For example:

    sol(t; idxs=@obsex(arg = atan(VIndex(1,:u_i), VIndex(1,:u_r))]
    sol(t; idxs=@obsex(δrel = VIndex(1,:δ) - VIndex(2,:δ)))

"""
macro obsex(ex)
    generate_observable_expression(ex)
end

function generate_observable_expression(::Any)
    error("@obsex can only be used when Symbolics.jl is loaded.")
end

# define function stub to overload in SymbolicsExt
function collect_symbol! end

function SII.is_observed(nw::Network, obsex::ObservableExpression)
    true
end

SII.symbolic_type(::Type{<:ObservableExpression}) = SII.ScalarSymbolic()
SII.hasname(::ObservableExpression) = true
function SII.getname(obsex::ObservableExpression)
    if obsex.name isa Symbol
        return obsex.name
    else
        io = IOBuffer()
        show(io, MIME"text/plain"(), obsex.ex)
        str = String(take!(io))
        return Symbol(replace(str, " "=>""))
    end
end

function SII.observed(nw::Network, obsex::ObservableExpression)
    inputf = SII.observed(nw, obsex.inputs)
    (u, p, t) -> begin
        input = inputf(u, p, t)
        obsex.f(input)
    end
end
function SII.observed(nw::Network, obsexs::AbstractVector{<:ObservableExpression})
    inputfs = map(obsex -> SII.observed(nw, obsex.inputs), obsexs)
    (u, p, t) -> begin
        map(obsexs, inputfs) do obsex, inputf
            input = inputf(u, p, t)
            obsex.f(input)
        end
    end
end

Base.getindex(s::NWState, idx::ObservableExpression) = SII.getu(s, idx)(s)
Base.getindex(s::NWParameter, idx::ObservableExpression) = SII.getp(s, idx)(s)

# using getindex to access component models
function Base.getindex(nw::Network, i::EIndex{<:Union{<:Pair,Symbol,Int}, Nothing})
    return getcomp(nw, i)
end
function Base.getindex(nw::Network, i::VIndex{<:Union{Symbol,Int}, Nothing})
    return getcomp(nw, i)
end
