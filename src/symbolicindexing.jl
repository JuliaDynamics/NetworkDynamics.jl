abstract type SymbolicIndex{C,S} end
abstract type SymbolicStateIndex{C,S} <: SymbolicIndex{C,S} end
abstract type SymbolicParameterIndex{C,S} <: SymbolicIndex{C,S} end
struct VIndex{C,S} <: SymbolicStateIndex{C,S}
    compidx::C
    subidx::S
end
struct EIndex{C,S} <: SymbolicStateIndex{C,S}
    compidx::C
    subidx::S
end
struct VPIndex{C,S} <: SymbolicParameterIndex{C,S}
    compidx::C
    subidx::S
end
struct EPIndex{C,S} <: SymbolicParameterIndex{C,S}
    compidx::C
    subidx::S
end

#=
XXX: SciMLBase Issue regarding f.sys
SciMLBase gets the index provider from ODEFunction.sys which defaults to f.sys so I provide it...
# SII.symbolic_container(odef::SciMLBase.ODEFunction{<:Any,<:Any,<:Network}) = odef.f
=#
SciMLBase.__has_sys(nw::Network) = true
Base.getproperty(nw::Network, s::Symbol) = s===:sys ? nw : getfield(nw, s)

SII.symbolic_type(::Type{<:SymbolicIndex{Int,<:Union{Symbol,Int}}}) = SII.ScalarSymbolic()
SII.symbolic_type(::Type{<:SymbolicIndex}) = SII.ArraySymbolic()


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
#### State as value provider
####

struct NWParameter{P,NW}
    nw::NW
    pflat::P
end

Base.eltype(p::NWParameter) = eltype(p.pflat)
Base.length(s::NWParameter) = length(s.pflat)

function NWParameter(nw::Network; ptype=Vector{Union{Float64,Nothing}}, default=true)
    pflat = _flat(ptype, pdim(nw))
    p = NWParameter(nw,pflat)
    default || return p
    for (k, v) in SII.default_values(nw)
        k isa Union{VPIndex, EPIndex} || continue
        p[k] = v
    end
    return p
end

function NWParameter(p::NWParameter; ptype=typeof(p.pflat))
    NWParameter(p.nw, _convertorcopy(ptype,pflat(p)))
end

struct NWState{U,P,T,NW}
    nw::NW
    uflat::U
    p::P
    t::T
    function NWState(nw, uflat, p=nothing, t=nothing)
        _p = (!indexable(p) || p isa NWParameter) ? p : NWParameter(nw,p)
        s = new{typeof(uflat),typeof(_p),typeof(t),typeof(nw)}(nw,uflat,_p,t)
        @argcheck !indexable(p) || s.nw === s.p.nw
        return s
    end
end

function NWState(nw::Network;
                 utype=Vector{Union{Float64,Nothing}},
                 ptype=Vector{Union{Float64,Nothing}},
                 default=true)
    t = nothing
    uflat = _flat(utype, dim(nw))
    p = NWParameter(nw; ptype, default=false)
    s = NWState(nw,uflat,p,t)
    default || return s
    for (k, v) in SII.default_values(nw)
        s[k] = v
    end
    return s
end

function NWState(s::NWState;
                 utype=typeof(utype(s)),
                 ptype=typeof(pflat(s)))
    p = NWParameter(s.p; ptype)
    NWState(s.nw, _convertorcopy(utype,uflat(s)), p, s.t)
end

function _flat(T, N)
    if Nothing <: eltype(T)
        vec = T(undef, N)
        fill!(vec, nothing)
    else
        zeros(eltype(T), N)
    end
end

_convertorcopy(::Type{T}, x::T) where {T} = copy(x)
function _convertorcopy(T, x)
    # _x = similar(T, length(x))
    # _x .= x
    # T(undef, length(x)) .= x
    convert(T, x)
end

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

Base.getindex(p::NWParameter, idx::SymbolicStateIndex) = getindex(p, _paraindex(idx))
Base.setindex!(p::NWParameter, val, idx::SymbolicStateIndex) = setindex!(p, val, _paraindex(idx))
Base.getindex(p::NWParameter, idx::SymbolicParameterIndex) = SII.getp(p, idx)(p)
Base.setindex!(p::NWParameter, val, idx::SymbolicParameterIndex) = SII.setp(p, idx)(p, val)
_paraindex(idx::VIndex) = VPIndex(idx.compidx, idx.subidx)
_paraindex(idx::EIndex) = EPIndex(idx.compidx, idx.subidx)

Base.getindex(s::NWState, idx::SymbolicParameterIndex) = SII.getp(s, idx)(s)
Base.setindex!(s::NWState, val, idx::SymbolicParameterIndex) = SII.setp(s, idx)(s, val)
Base.getindex(s::NWState, idx::SymbolicStateIndex) = SII.getu(s, idx)(s)
Base.setindex!(s::NWState, val, idx::SymbolicStateIndex) = SII.setu(s, idx)(s, val)
function Base.getindex(s::NWState, idxs)
    et = eltype(idxs)
    if et <: SymbolicStateIndex
        SII.getu(s, idxs)(s)
    elseif et <: SymbolicParameterIndex
        SII.getp(s, idxs)(s)
    else
        getindex.(Ref(s), idxs)
    end
end

Base.getindex(s::NWState, ::Colon) = uflat(s)
Base.getindex(p::NWParameter, ::Colon) = pflat(p)
Base.setindex!(s::NWState, val, ::Colon) = uflat(s) .= val
Base.setindex!(p::NWParameter, val, ::Colon) = pflat(p) .= val

struct VProxy{S} s::S end
Base.getindex(p::VProxy, comp, state) = getindex(p.s, VIndex(comp, state))
Base.setindex!(p::VProxy, val, comp, state) = setindex!(p.s, val, VIndex(comp, state))
struct EProxy{S} s::S end
Base.getindex(p::EProxy, comp, state) = getindex(p.s, EIndex(comp, state))
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

#=
https://discourse.julialang.org/t/broadcasting-setindex-is-a-noobtrap/94700
its hard to overload .= braodcasting so lets error if somebody tries
=#
Base.dotview(s::Union{NWParameter, NWState, VProxy, EProxy}, idxs...) = view(s, idxs...)
function Base.view(s::Union{NWParameter, NWState, VProxy, EProxy}, idxs...)
    error("Cannot create view into for indice $idxs")
end
Base.view(s::NWState, ::Colon) = s.uflat
Base.view(s::NWParameter, ::Colon) = s.pflat

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
