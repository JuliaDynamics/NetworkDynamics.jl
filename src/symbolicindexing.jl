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

#
Base.broadcastable(si::SymbolicIndex{<:Union{Int,Colon},<:Union{Int,Symbol,Colon}}) = Ref(si)

Base.length(si::SymbolicIndex{<:Union{AbstractVector,Tuple},<:Union{Int,Symbol}}) = length(si.compidx)
function Base.eltype(si::SymbolicIndex{<:Union{AbstractVector,Tuple},<:Union{Int,Symbol}})
    if isconcretetype(eltype(si.compidx))
        _baseT(si){eltype(si.compidx),typeof(si.subidx)}
    else
        Any
    end
end
function Base.iterate(si::SymbolicIndex{<:Union{AbstractVector,Tuple},<:Union{Int,Symbol}}, state=nothing)
    it = isnothing(state) ? iterate(si.compidx) : iterate(si.compidx, state)
    isnothing(it) && return nothing
    _similar(si, it[1], si.subidx), it[2]
end

Base.length(si::SymbolicIndex{Int,<:Union{AbstractVector,Tuple}}) = length(si.subidx)
function Base.eltype(si::SymbolicIndex{Int,<:Union{AbstractVector,Tuple}})
    if isconcretetype(eltype(si.subidx))
        _baseT(si){eltype(si.compidx),eltype(si.subidx)}
    else
        Any
    end
end
function Base.iterate(si::SymbolicIndex{Int,<:Union{AbstractVector,Tuple}}, state=nothing)
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
_resolve_colon(nw::Network, sni::VIndex{Int,Colon}) = VIndex(sni.compidx, 1:dim(getcomp(nw,sni)))
_resolve_colon(nw::Network, sni::EIndex{Int,Colon}) = EIndex(sni.compidx, 1:dim(getcomp(nw,sni)))
_resolve_colon(nw::Network, sni::VPIndex{Int,Colon}) = VPIndex(sni.compidx, 1:pdim(getcomp(nw,sni)))
_resolve_colon(nw::Network, sni::EPIndex{Int,Colon}) = EPIndex(sni.compidx, 1:pdim(getcomp(nw,sni)))


#### Implmentation of index provider interface
####
#### Structural things
####
SII.symbolic_container(nw::Network) = nw
SII.is_independent_variable(nw::Network, _) = false
SII.independent_variable_symbols(nw::Network) = []
SII.is_time_dependent(nw::Network) = true
SII.constant_structure(::Network) = true
SII.all_variable_symbols(nw::Network) = vcat(SII.variable_symbols(nw), observed_symbols(nw))
SII.all_symbols(nw::Network) = vcat(SII.all_variable_symbols(nw), SII.parameter_symbols(nw))


####
#### variable indexing
####
function SII.is_variable(nw::Network, sni)
    _sni = _resolve_colon.(nw,sni)
    all(_is_variable.(nw, _sni))
end
_is_variable(nw::Network, sni) = false
function _is_variable(nw::Network, sni::SymbolicStateIndex{Int,<:Union{Int,Symbol}})
    cf = getcomp(nw, sni)
    return isdynamic(cf) && subsym_has_idx(sni.subidx, cf.sym)
end

function SII.variable_index(nw::Network, sni)
    _sni = _resolve_colon.(nw, sni)
     _variable_index.(nw, _sni)
end
function _variable_index(nw::Network, sni::SymbolicStateIndex{Int,<:Union{Int,Symbol}})
    cf = getcomp(nw, sni)
    range = getcomprange(nw, sni)
    range[subsym_to_idx(sni.subidx, cf.sym)]
end

function SII.variable_symbols(nw::Network)
    syms = Vector{SymbolicStateIndex{Int,Symbol}}(undef, dim(nw))
    for (ci, cf) in pairs(nw.im.vertexf)
        isdynamic(cf) || continue
        syms[nw.im.v_data[ci]] .= VIndex.(ci, cf.sym)
    end
    for (ci, cf) in pairs(nw.im.edgef)
        isdynamic(cf) || continue
        syms[nw.im.e_data[ci]] .= EIndex.(ci, cf.sym)
    end
    return syms
end


####
#### parameter indexing
####
function SII.is_parameter(nw::Network, sni)
    _sni = _resolve_colon.(nw,sni)
    all(_is_parameter.(nw, _sni))
end
_is_parameter(nw::Network, sni) = false
function _is_parameter(nw::Network,
                          sni::SymbolicParameterIndex{Int,<:Union{Int,Symbol}})
    cf = getcomp(nw, sni)
    return subsym_has_idx(sni.subidx, cf.psym)
end

function SII.parameter_index(nw::Network, sni)
    _sni = _resolve_colon.(nw, sni)
     _parameter_index.(nw, _sni)
end
function _parameter_index(nw::Network,
                             sni::SymbolicParameterIndex{Int,<:Union{Int,Symbol}})
    cf = getcomp(nw, sni)
    range = getcompprange(nw, sni)
    range[subsym_to_idx(sni.subidx, cf.psym)]
end

function SII.parameter_symbols(nw::Network)
    syms = Vector{SymbolicParameterIndex{Int,Symbol}}(undef, pdim(nw))
    i = 1
    for (ci, cf) in pairs(nw.im.vertexf)
        syms[nw.im.v_para[ci]] .= VPIndex.(ci, cf.psym)
    end
    for (ci, cf) in pairs(nw.im.edgef)
        syms[nw.im.e_para[ci]] .= EPIndex.(ci, cf.psym)
    end
    return syms
end

####
#### Observed indexing
function SII.is_observed(nw::Network, sni)
    if !_hascolon(sni)
        all(_is_observed.(nw, sni))
    else
        # if has colon check if all are observed OR variables and return true
        # the observed function will handle the whole thing then
        _sni = _resolve_colon.(nw,sni)
        all(s -> _is_observed(nw,s) || _is_variable(nw,s), _sni)
    end
end
_is_observed(nw::Network, _) = false
function _is_observed(nw::Network, sni::SymbolicStateIndex{Int,<:Union{Int,Symbol}})
    cf = getcomp(nw, sni)

    if isdynamic(cf)
        return sni.subidx ∈ cf.obssym # only works if symbol is given
    else
        return sni.subidx ∈ cf.obssym || subsym_has_idx(sni.subidx, cf.sym)
    end
    error("Could not handle $sni. Should never be reached.")
end

function observed_symbols(nw::Network)
    syms = SymbolicStateIndex{Int,Symbol}[]
    for (ci, cf) in pairs(nw.im.vertexf)
        if !isdynamic(cf)
            for s in cf.sym
                push!(syms, VIndex(ci, s))
            end
        end
        for s in cf.obssym
            push!(syms, VIndex(ci, s))
        end
    end
    for (ci, cf) in pairs(nw.im.edgef)
        if !isdynamic(cf)
            for s in cf.sym
                push!(syms, EIndex(ci, s))
            end
        end
        for s in cf.obssym
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
    _snis = _resolve_colon.(nw, snis)

    # First: resolve everything in fullstate (static and dynamic states)
    flatidxs = broadcast(_snis) do sni
        if SII.is_variable(nw, sni)
            SII.variable_index(nw, sni)
        else
            cf = getcomp(nw, sni)
            if !isdynamic(cf) && subsym_has_idx(sni.subidx, cf.sym)
                _range = getcomprange(nw, sni)
                _range[subsym_to_idx(sni.subidx, cf.sym)]
            else
                error("Observalbe mechanism cannot handle $sni yet.")
            end
        end
    end

    function(u, p, t)
        du = nw.cachepool[u]
        nw(du, u, p, t)
        # XXX: split coreloop in static/dynamic parts and don't rely on the same buffer
        _u = nw.cachepool[u, nw.im.lastidx_static]
        _u[flatidxs]
    end
end


####
#### Default values
####
function SII.default_values(nw::Network)
    defs = Dict{SymbolicIndex{Int,Symbol},Float64}()
    for (ci, cf) in pairs(nw.im.vertexf)
        for (s, def) in zip(cf.psym, cf.pdef)
            isnothing(def) && continue
            defs[VPIndex(ci, s)] = def
        end
        for (s, def) in zip(cf.sym, cf.def)
            isnothing(def) && continue
            defs[VIndex(ci, s)] = def
        end
    end
    for (ci, cf) in pairs(nw.im.edgef)
        for (s, def) in zip(cf.psym, cf.pdef)
            isnothing(def) && continue
            defs[EPIndex(ci, s)] = def
        end
        for (s, def) in zip(cf.sym, cf.def)
            isnothing(def) && continue
            defs[EIndex(ci, s)] = def
        end
    end
    return defs
end


####
#### State as value provider
####
struct State{U,P,NW,T}
    nw::NW
    uflat::U
    pflat::P
    t::T
end

function State(nw::Network)
    t = nothing
    uflat = Union{Nothing, Float64}[nothing for _ in 1:dim(nw)]
    pflat = Union{Nothing, Float64}[nothing for _ in 1:pdim(nw)]
    s = State(nw,uflat,pflat,t)
    for (k, v) in SII.default_values(nw)
        s[k] = v
    end
    s
end

Base.broadcastable(s::State) = Ref(s)

SII.symbolic_container(s::State) = s.nw
SII.state_values(s::State) = s.uflat
SII.parameter_values(s::State) = s.pflat
SII.current_time(s::State) = s.t

Base.getindex(s::State, idx::SymbolicParameterIndex) = SII.getp(s, idx)(s)
Base.setindex!(s::State, val, idx::SymbolicParameterIndex) = SII.setp(s, idx)(s, val)
Base.getindex(s::State, idx::SymbolicStateIndex) = SII.getu(s, idx)(s)
Base.setindex!(s::State, val, idx::SymbolicStateIndex) = SII.setu(s, idx)(s, val)
function Base.getindex(s::State, idxs)
    et = eltype(idxs)
    if et <: SymbolicStateIndex
        SII.getu(s, idxs)(s)
    elseif et <: SymbolicParameterIndex
        SII.getp(s, idxs)(s)
    else
        getindex.(Ref(s), idxs)
    end
end

struct VStateProxy{S} s::S end
Base.getindex(p::VStateProxy, comp, state) = getindex(p.s, VIndex(comp, state))
Base.setindex!(p::VStateProxy, val, comp, state) = setindex!(p.s, val, VIndex(comp, state))
struct EStateProxy{S} s::S end
Base.getindex(p::EStateProxy, comp, state) = getindex(p.s, EIndex(comp, state))
Base.setindex!(p::EStateProxy, val, comp, state) = setindex!(p.s, val, EIndex(comp, state))
struct VParameterProxy{S} s::S end
Base.getindex(p::VParameterProxy, comp, state) = getindex(p.s, VPIndex(comp, state))
Base.setindex!(p::VParameterProxy, val, comp, state) = setindex!(p.s, val, VPIndex(comp, state))
struct EParameterProxy{S} s::S end
Base.getindex(p::EParameterProxy, comp, state) = getindex(p.s, EPIndex(comp, state))
Base.setindex!(p::EParameterProxy, val, comp, state) = setindex!(p.s, val, EPIndex(comp, state))

function Base.getproperty(s::State, sym::Symbol)
    if sym === :v
        return VStateProxy(s)
    elseif sym === :e
        return EStateProxy(s)
    elseif sym === :pv
        return VParameterProxy(s)
    elseif sym === :pe
        return EParameterProxy(s)
    else
        return getfield(s, sym)
    end
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
