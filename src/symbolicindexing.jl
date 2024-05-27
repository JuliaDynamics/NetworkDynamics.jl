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
getcomprange(nw::Network, sni::VPIndex{Int}) = nw.im.v_para[sni.compidx]
getcomprange(nw::Network, sni::EPIndex{Int}) = nw.im.e_para[sni.compidx]

subsym_has_idx(sym::Symbol, syms) = sym ∈ syms
subsym_has_idx(idx::Int, syms) = 1 ≤ idx ≤ length(syms)
subsym_to_idx(sym::Symbol, syms) = findfirst(isequal(sym), syms)
subsym_to_idx(idx::Int, _) = idx

SII.symbolic_container(nw::Network) = nw

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
    i = 1
    for (ci, cf) in pairs(nw.im.vertexf)
        isdynamic(cf) || continue
        for s in cf.sym
            syms[i] = VIndex(ci, s)
            i = i + 1
        end
    end
    for (ci, cf) in pairs(nw.im.edgef)
        isdynamic(cf) || continue
        for s in cf.sym
            syms[i] = EIndex(ci, s)
            i = i + 1
        end
    end
    return syms
end

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
    syms = Vector{SymbolicStateIndex{Int,Symbol}}(undef, pdim(nw))
    i = 1
    for (ci, cf) in pairs(nw.im.vertexf)
        for s in cf.psym
            syms[i] = VIndex(ci, s)
            i = i + 1
        end
    end
    for (ci, cf) in pairs(nw.im.edgef)
        for s in cf.psym
            syms[i] = EIndex(ci, s)
            i = i + 1
        end
    end
    return syms
end

SII.is_independent_variable(nw::Network, _) = false

SII.independent_variable_symbols(nw::Network) = []

SII.is_time_dependent(nw::Network) = true

SII.constant_structure(::Network) = true

SII.is_observed(nw::Network, _) = false

# function SII.all_solvable_symbols(nw::Network)
#   return vcat(
#     collect(keys(sys.state_index)),
#     collect(keys(sys.observed)),
#   )
# end

# function SII.all_symbols(nw::Network)
#   return vcat(
#     all_solvable_symbols(sys),
#     collect(keys(sys.parameter_index)),
#     sys.independent_variable === nothing ? Symbol[] : sys.independent_variable
#   )
# end

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
