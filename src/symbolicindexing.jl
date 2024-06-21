abstract type SymbolicIndex{C,S} end
abstract type SymbolicStateIndex{C,S} <: SymbolicIndex{C,S} end
abstract type SymbolicParameterIndex{C,S} <: SymbolicIndex{C,S} end
"""
    VIndex{C,S} <: SymbolicStateIndex{C,S}
    idx = VIndex(comp, sub)

A symbolic index for a vertex state variable.
- `comp`: the component index, either int or a collection of ints
- `sub`: the subindex, either int, symbol or a collection of those.

```
VIndex(1, :P)      # vertex 1, variable :P
VIndex(1:5, 1)     # first state of vertices 1 to 5
VIndex(7, (:x,:y)) # states :x and :y of vertex 7
```

Can be used to index into objects supporting the `SymbolicIndexingInterface`,
e.g. [`NWState`](@ref), [`NWParameter`](@ref) or `ODESolution`.

See also: [`EIndex`](@ref), [`VPIndex`](@ref), [`EPIndex`](@ref)
"""
struct VIndex{C,S} <: SymbolicStateIndex{C,S}
    compidx::C
    subidx::S
end
"""
    EIndex{C,S} <: SymbolicStateIndex{C,S}
    idx = EIndex(comp, sub)

A symbolic index for an edge state variable.
- `comp`: the component index, either int or a collection of ints
- `sub`: the subindex, either int, symbol or a collection of those.

```
EIndex(1, :P)      # edge 1, variable :P
EIndex(1:5, 1)     # first state of edges 1 to 5
EIndex(7, (:x,:y)) # states :x and :y of edge 7
```

Can be used to index into objects supporting the `SymbolicIndexingInterface`,
e.g. [`NWState`](@ref), [`NWParameter`](@ref) or `ODESolution`.

See also: [`VIndex`](@ref), [`VPIndex`](@ref), [`EPIndex`](@ref)
"""
struct EIndex{C,S} <: SymbolicStateIndex{C,S}
    compidx::C
    subidx::S
end
"""
    VPIndex{C,S} <: SymbolicStateIndex{C,S}
    idx = VPIndex(comp, sub)

A symbolic index into the parameter a vertex:
- `comp`: the component index, either int or a collection of ints
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
    VEIndex{C,S} <: SymbolicStateIndex{C,S}
    idx = VEIndex(comp, sub)

A symbolic index into the parameter a vertex:
- `comp`: the component index, either int or a collection of ints
- `sub`: the subindex, either int, symbol or a collection of those.

Can be used to index into objects supporting the `SymbolicIndexingInterface`,
e.g. [`NWParameter`](@ref) or `ODEProblem`.

See also: [`VPIndex`](@ref), [`VIndex`](@ref), [`EIndex`](@ref)
"""
struct EPIndex{C,S} <: SymbolicParameterIndex{C,S}
    compidx::C
    subidx::S
end
const SymbolicEdgeIndex = Union{EIndex, EPIndex}
const SymbolicVertexIndex = Union{VIndex, VPIndex}

#=
SciMLBase gets the index provider from ODEFunction.sys which defaults to f.sys so we provide it...
SSI Maintainer assured that f.sys is really only used for symbolic indexig so method seems legit
=#
SciMLBase.__has_sys(nw::Network) = true
Base.getproperty(nw::Network, s::Symbol) = s===:sys ? nw : getfield(nw, s)

SII.symbolic_type(::Type{<:SymbolicIndex{Int,<:Union{Symbol,Int}}}) = SII.ScalarSymbolic()
SII.symbolic_type(::Type{<:SymbolicIndex}) = SII.ArraySymbolic()

SII.hasname(::SymbolicIndex) = false
SII.hasname(::SymbolicIndex{Int,<:Union{Symbol,Int}}) = true
SII.getname(x::SymbolicVertexIndex) = Symbol("v$(x.compidx)₊$(x.subidx)")
SII.getname(x::SymbolicEdgeIndex) = Symbol("e$(x.compidx)₊$(x.subidx)")

getcomp(nw::Network, sni::Union{EIndex{Int},EPIndex{Int}}) = nw.im.edgef[sni.compidx]
getcomp(nw::Network, sni::Union{VIndex{Int},VPIndex{Int}}) = nw.im.vertexf[sni.compidx]
getcomprange(nw::Network, sni::VIndex{Int}) = nw.im.v_data[sni.compidx]
getcomprange(nw::Network, sni::EIndex{Int}) = nw.im.e_data[sni.compidx]
getcompprange(nw::Network, sni::VPIndex{Int}) = nw.im.v_para[sni.compidx]
getcompprange(nw::Network, sni::EPIndex{Int}) = nw.im.e_para[sni.compidx]

subsym_has_idx(sym::Symbol, syms) = sym ∈ syms
subsym_has_idx(idx::Int, syms) = 1 ≤ idx ≤ length(syms)
subsym_to_idx(sym::Symbol, syms) = findfirst(isequal(sym), syms)
subsym_to_idx(idx::Int, _) = idx

####
#### Iterator/Broadcast interface for ArraySymbolic types
####
Base.broadcastable(si::SymbolicIndex{<:Union{Int,Colon},<:Union{Int,Symbol,Colon}}) = Ref(si)

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

const _IterableSubcomponent = SymbolicIndex{Int,<:Union{AbstractVector,Tuple}}
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
_resolve_colon(nw::Network, sni::VIndex{Int,Colon}) = VIndex{Int, UnitRange{Int}}(sni.compidx, 1:dim(getcomp(nw,sni)))
_resolve_colon(nw::Network, sni::EIndex{Int,Colon}) = EIndex{Int, UnitRange{Int}}(sni.compidx, 1:dim(getcomp(nw,sni)))
_resolve_colon(nw::Network, sni::VPIndex{Int,Colon}) = VPIndex{Int, UnitRange{Int}}(sni.compidx, 1:pdim(getcomp(nw,sni)))
_resolve_colon(nw::Network, sni::EPIndex{Int,Colon}) = EPIndex{Int, UnitRange{Int}}(sni.compidx, 1:pdim(getcomp(nw,sni)))


#### Implmentation of index provider interface
####
#### Structural things
####
SII.symbolic_container(nw::Network) = nw
SII.is_independent_variable(nw::Network, sym) = sym == :t
SII.independent_variable_symbols(nw::Network) = [:t]
SII.is_time_dependent(nw::Network) = true
SII.constant_structure(::Network) = true
SII.all_variable_symbols(nw::Network) = vcat(SII.variable_symbols(nw), observed_symbols(nw))
SII.all_symbols(nw::Network) = vcat(SII.all_variable_symbols(nw), SII.parameter_symbols(nw))


####
#### variable indexing
####
function SII.is_variable(nw::Network, sni)
    if _hascolon(sni)
        SII.is_variable(nw, _resolve_colon(nw,sni))
    elseif SII.symbolic_type(sni) === SII.ArraySymbolic()
        all(s -> SII.is_variable(nw, s), sni)
    else
        _is_variable(nw, sni)
    end
end
_is_variable(nw::Network, sni) = false
function _is_variable(nw::Network, sni::SymbolicStateIndex{Int,<:Union{Int,Symbol}})
    cf = getcomp(nw, sni)
    return isdynamic(cf) && subsym_has_idx(sni.subidx, sym(cf))
end

function SII.variable_index(nw::Network, sni)
    if _hascolon(sni)
        SII.variable_index(nw, _resolve_colon(nw,sni))
    elseif SII.symbolic_type(sni) === SII.ArraySymbolic()
        SII.variable_index.(nw, sni)
    else
        _variable_index(nw, sni)
    end
end
function _variable_index(nw::Network, sni::SymbolicStateIndex{Int,<:Union{Int,Symbol}})
    cf = getcomp(nw, sni)
    range = getcomprange(nw, sni)
    range[subsym_to_idx(sni.subidx, sym(cf))]
end

function SII.variable_symbols(nw::Network)
    syms = Vector{SymbolicStateIndex{Int,Symbol}}(undef, dim(nw))
    for (ci, cf) in pairs(nw.im.vertexf)
        isdynamic(cf) || continue
        syms[nw.im.v_data[ci]] .= VIndex.(ci, sym(cf))
    end
    for (ci, cf) in pairs(nw.im.edgef)
        isdynamic(cf) || continue
        syms[nw.im.e_data[ci]] .= EIndex.(ci, sym(cf))
    end
    return syms
end


####
#### parameter indexing
####
function SII.is_parameter(nw::Network, sni)
    if _hascolon(sni)
        SII.is_parameter(nw, _resolve_colon(nw,sni))
    elseif SII.symbolic_type(sni) === SII.ArraySymbolic()
        all(s -> SII.is_parameter(nw, s), sni)
    else
        _is_parameter(nw, sni)
    end
end
_is_parameter(nw::Network, sni) = false
function _is_parameter(nw::Network,
                          sni::SymbolicParameterIndex{Int,<:Union{Int,Symbol}})
    cf = getcomp(nw, sni)
    return subsym_has_idx(sni.subidx, psym(cf))
end

# SII.is_timeseries_parameter(nw::Network, sym) = false

function SII.parameter_index(nw::Network, sni)
    if _hascolon(sni)
        SII.parameter_index(nw, _resolve_colon(nw,sni))
    elseif SII.symbolic_type(sni) === SII.ArraySymbolic()
        SII.parameter_index.(nw, sni)
    else
        _parameter_index(nw, sni)
    end
end
function _parameter_index(nw::Network,
                             sni::SymbolicParameterIndex{Int,<:Union{Int,Symbol}})
    cf = getcomp(nw, sni)
    range = getcompprange(nw, sni)
    range[subsym_to_idx(sni.subidx, psym(cf))]
end

function SII.parameter_symbols(nw::Network)
    syms = Vector{SymbolicParameterIndex{Int,Symbol}}(undef, pdim(nw))
    i = 1
    for (ci, cf) in pairs(nw.im.vertexf)
        syms[nw.im.v_para[ci]] .= VPIndex.(ci, psym(cf))
    end
    for (ci, cf) in pairs(nw.im.edgef)
        syms[nw.im.e_para[ci]] .= EPIndex.(ci, psym(cf))
    end
    return syms
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
    else
        _is_observed(nw, sni)
    end
end
_is_observed(nw::Network, _) = false
function _is_observed(nw::Network, sni::SymbolicStateIndex{Int,<:Union{Int,Symbol}})
    cf = getcomp(nw, sni)

    if isdynamic(cf)
        return sni.subidx ∈ obssym(cf) # only works if symbol is given
    else
        return sni.subidx ∈ obssym(cf) || subsym_has_idx(sni.subidx, sym(cf))
    end
    error("Could not handle $sni. Should never be reached.")
end

function observed_symbols(nw::Network)
    syms = SymbolicStateIndex{Int,Symbol}[]
    for (ci, cf) in pairs(nw.im.vertexf)
        if !isdynamic(cf)
            for s in sym(cf)
                push!(syms, VIndex(ci, s))
            end
        end
        for s in obssym(cf)
            push!(syms, VIndex(ci, s))
        end
    end
    for (ci, cf) in pairs(nw.im.edgef)
        if !isdynamic(cf)
            for s in sym(cf)
                push!(syms, EIndex(ci, s))
            end
        end
        for s in obssym(cf)
            push!(syms, EIndex(ci, s))
        end
    end
    return syms
end

#=
Types of observable calls:
- observable of ODEVertex
  - recalculate incident static edges
  - aggregate locally
  - execut observable
- State of StaticVertex
  - aggregate locally
  - execute vertex
- observable of StaticVertex
  - aggregate locallay
  - execute vertex
  - execute observable

- observable of ODEEdge
  - if incident vertex is static: locally aggregate & execute vertex function
  - execute observable
- State of StaticEdge
  - execute edge (cannot depend on static vertices)
- Observable of static edge
  - execute edge
  - execute observable
=#

function SII.observed(nw::Network, snis)
    _snis = _expand_and_collect(nw, snis)

    # First: resolve everything in fullstate (static and dynamic states)
    flatidxs = broadcast(_snis) do sni
        if SII.is_variable(nw, sni)
            SII.variable_index(nw, sni)
        else
            cf = getcomp(nw, sni)
            if !isdynamic(cf) && subsym_has_idx(sni.subidx, sym(cf))
                _range = getcomprange(nw, sni)
                _range[subsym_to_idx(sni.subidx, sym(cf))]
            else
                error("Observalbe mechanism cannot handle $sni yet.")
            end
        end
    end

    let _nw=nw, _flatidxs=flatidxs
        function(u, p, t)
            du = _nw.cachepool[u]
            _nw(du, u, p, t)
            # XXX: split coreloop in static/dynamic parts and don't rely on the same buffer
            _u = _nw.cachepool[u, _nw.im.lastidx_static]
            _u[_flatidxs]
        end
    end
end
function _expand_and_collect(nw::Network, sni::SymbolicIndex)
    if _hascolon(sni)
        collect(_resolve_colon(nw, sni))
    elseif SII.symbolic_type(sni) === SII.ArraySymbolic()
        collect(sni)
    else
        sni
    end
end
function _expand_and_collect(nw::Network, snis)
    mapreduce(vcat, snis) do sni
        _expand_and_collect(nw, sni)
    end
end


####
#### Default values
####
function SII.default_values(nw::Network)
    defs = Dict{SymbolicIndex{Int,Symbol},Float64}()
    for (ci, cf) in pairs(nw.im.vertexf)
        for (s, def) in zip(psym(cf), pdef(cf))
            isnothing(def) && continue
            defs[VPIndex(ci, s)] = def
        end
        for (s, def) in zip(sym(cf), def(cf))
            isnothing(def) && continue
            defs[VIndex(ci, s)] = def
        end
    end
    for (ci, cf) in pairs(nw.im.edgef)
        for (s, def) in zip(psym(cf), pdef(cf))
            isnothing(def) && continue
            defs[EPIndex(ci, s)] = def
        end
        for (s, def) in zip(sym(cf), def(cf))
            isnothing(def) && continue
            defs[EIndex(ci, s)] = def
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

Get flat array representation using `pflat(p)`.
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

Get flat array representation using `uflat(s)` and `pflat(s)`.
"""
struct NWState{U,P,T,NW<:Network}
    nw::NW
    uflat::U
    p::P
    t::T
    function NWState(thing, uflat, p=nothing, t=nothing)
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
uflat(s::NWState) = s.uflat
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
Base.getindex(p::NWParameter, idx) = SII.getp(p, _paraindex(idx))(p)

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
Base.getindex(s::NWState, idx::SymbolicParameterIndex) = SII.getp(s, idx)(s)
Base.getindex(s::NWState, idx::SymbolicStateIndex) = SII.getu(s, idx)(s)
function Base.getindex(s::NWState, idxs)
    if all(i -> i isa SymbolicParameterIndex, idxs)
        SII.getp(s, idxs)(s)
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
    _get_components(IT, inpr)[cidx].sym
end
function _make_sidx_iterable(IT::Type{<:SymbolicParameterIndex}, inpr, cidx, ::Colon)
    _get_components(IT, inpr)[cidx].psym
end
function _make_sidx_iterable(IT::Type{<:SymbolicStateIndex}, inpr, cidx, s::Union{AbstractString,AbstractPattern})
    comp = _get_components(IT, inpr)[cidx]
    syms = vcat(comp.sym, comp.obssym)
    filter(sym -> contains(string(sym), s), syms)
end
function _make_sidx_iterable(IT::Type{<:SymbolicParameterIndex}, inpr, cidx, s::Union{AbstractString,AbstractPattern})
    syms = _get_components(IT, inpr)[cidx].psym
    filter(sym -> contains(string(sym), s), syms)
end

_get_components(::Type{<:SymbolicVertexIndex}, inpr) = extract_nw(inpr).im.vertexf
_get_components(::Type{<:SymbolicEdgeIndex}, inpr) = extract_nw(inpr).im.edgef


"""
    extract_nw(thing)

Try to extract the `Network` object from thing.
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

#=
nds = wrap(nd, u, [p]) -> NWState (contains nw para, optional für observables/static)
ndp = wrap(nd, p) -> NWPara

u = unwrap(nds)

s = NWState(nd) -> initial guess aus den component functions
NWPara(nd) -> default parameter

s = NWState(nd, u, [p]) # p nur für observables/static
s.e[j, idx/:i_r]
s.v[i, idx/:u_r]
s[:edge, j, idx/:i_r]
s[:vertex, i, idx/:u_r]


p = NWPara(nd, [p])
p.pv[i, idx/:M]
p.pe[i, idx/:Y]



s = State(nd, u, p)
s.e[j, idx/:i_r]
s.v[i, idx/:u_r]
s.pv[i, idx/:M]
s.v.p[i, idx/:M] ?
s.pe[i, idx/:Y]
s.e.p[i, idx/:Y] ?

i_r_j(t) = State(nd, sol(t), pr(t)).e[j, :i_r]


State(nd, sol(t), pr(t))

nds = NDSol(sol, pr)
s.e[j, idx/:i_r](t)
s.v[i, idx/:u_r](t)
s.pv[i, idx/:M](t)
s.pe[i, idx/:Y](t)

s = nds(t)
=#
