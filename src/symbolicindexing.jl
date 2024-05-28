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
=#
SciMLBase.__has_sys(nw::Network) = true
Base.getproperty(nw::Network, s::Symbol) = s===:sys ? nw : getfield(nw, s)

SII.symbolic_type(::Type{<:SymbolicIndex{Int,<:Union{Symbol,Int}}}) = SII.ScalarSymbolic()
SII.symbolic_type(::Type{<:SymbolicIndex}) = SII.ArraySymbolic()

function Base.collect(sni::SymbolicIndex)
    multiple_comp = sni.compidx isa Union{AbstractVector,Tuple}
    multiple_sym = sni.subidx isa Union{AbstractVector,Tuple}
    if multiple_comp && multiple_sym
        error("Symbolic network indexing with multiple indices and multiple symbols not implemented!")
    end
    if multiple_comp
        return [_similar(sni, i, sni.subidx) for i in sni.compidx]
    end
    if multiple_sym
        return [_similar(sni, sni.compidx, i) for i in sni.subidx]
    end
    error("Shouldn't be reached, $sni does not seam to be a ArraySymbolic.")
end
_similar(::VIndex, c, s)  = VIndex(c,s)
_similar(::EIndex, c, s)  = EIndex(c,s)
_similar(::VPIndex, c, s) = VPIndex(c,s)
_similar(::EPIndex, c, s) = EPIndex(c,s)

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
SII.is_variable(nw::Network, sni) = false
function SII.is_variable(nw::Network,
                         sni::SymbolicStateIndex{Int,<:Union{Int,Symbol}})
    cf = getcomp(nw, sni)
    return isdynamic(cf) && subsym_has_idx(sni.subidx, cf.sym)
end

function SII.variable_index(nw::Network,
                            sni::SymbolicStateIndex{Int,<:Union{Int,Symbol}})
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
SII.is_parameter(nw::Network, sni) = false
function SII.is_parameter(nw::Network,
                          sni::SymbolicParameterIndex{Int,<:Union{Int,Symbol}})
    cf = getcomp(nw, sni)
    return subsym_has_idx(sni.subidx, cf.psym)
end

function SII.parameter_index(nw::Network,
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
####
SII.is_observed(nw::Network, _) = false
function SII.is_observed(nw::Network, sni::SymbolicStateIndex{Int,<:Union{Int,Symbol}})
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

SII.observed(nw::Network, sni) = error("Cant do for $sni")
function SII.observed(nw::Network, sni::SymbolicStateIndex{Int,<:Union{Int,Symbol}})
    SII.observed(nw, StaticArrays.SA[sni], true)
end
function SII.observed(nw::Network, snis::AbstractVector{<:SymbolicStateIndex}, singlevalue=false)
    flatidxs = map(snis) do sni
        cf = getcomp(nw, sni)
        if isdynamic(cf) || !(subsym_has_idx(sni.subidx, cf.sym))
            error("Currently the Observable mecahnism only supports the retrivel of static states.")
        end
        _range = getcomprange(nw, sni)
        _range[subsym_to_idx(sni.subidx, cf.sym)]
    end
    _flatidxs = singlevalue ? only(flatidxs) : flatidxs

    function(u, p, t)
        du = nw.cachepool[u]
        nw(du, u, p, t)
        # XXX: split coreloop in static/dynamic parts and don't rely on the same buffer
        _u = nw.cachepool[u, nw.im.lastidx_static]
        _u[_flatidxs]
    end
end


####
#### Default values
####
function SII.default_values(nw::Network)
    defs = Dict{SymbolicStateIndex{Int,Symbol},Float64}()
    for (ci, cf) in pairs(nw.im.vertexf)
        for (s, def) in zip(cf.psym, cf.pdef)
            isnothing(def) && continue
            defs[VIndex(ci, s)] = def
        end
        for (s, def) in zip(cf.sym, cf.def)
            isnothing(def) && continue
            defs[VIndex(ci, s)] = def
        end
    end
    for (ci, cf) in pairs(nw.im.edgef)
        for (s, def) in zip(cf.psym, cf.pdef)
            isnothing(def) && continue
            defs[EIndex(ci, s)] = def
        end
        for (s, def) in zip(cf.sym, cf.def)
            isnothing(def) && continue
            defs[EIndex(ci, s)] = def
        end
    end
    return defs
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
