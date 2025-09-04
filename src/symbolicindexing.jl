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
    defaults = SII.default_values(nw)
    for (i, sym) in enumerate(SII.parameter_symbols(nw))
        if haskey(defaults, sym)
            pflat[i] = defaults[sym]
        end
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
    defaults = SII.default_values(nw)
    for (i, k) in enumerate(SII.variable_symbols(nw))
        if haskey(defaults, k)
            uflat[i] = defaults[k]
        end
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

Copies both u and p objects so you cannot mess with the solution through NWState object.
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

    NWState(sol, copy(u), copy(p), t)
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
function Base.getindex(p::NWParameter, idx::Union{SymbolicIndex, AbstractVector, Tuple})
    SII.getp(p, _expand_and_collect(p, idx))(p)
end

# NWParameter: setindex!
function Base.setindex!(p::NWParameter, val, idx)
    setter = SII.setp(p, idx)
    _chk_dimensions(setter, val)
    setter(p, val)
end

# NWState: getindex
Base.getindex(s::NWState, ::Colon) = uflat(s)
function Base.getindex(s::NWState, idx::Union{SymbolicIndex, AbstractVector, Tuple})
    # for non-vector we can use getu regardless
    # for vector, SII will dispatch to getp when all are parameters, which requires _exapand_and_collect
    if SII.symbolic_type(idx) != SII.ScalarSymbolic() && SII.is_parameter(s, idx)
        getindex(s.p, idx) # this will expand and collect
    else
        SII.getu(s, idx)(s)
    end
end

# NWState: setindex!
function Base.setindex!(s::NWState, val, idx)
    setter = SII.setu(s, idx)
    _chk_dimensions(setter, val)
    setter(s, val)
end

function _chk_dimensions(::Union{SII.SetParameterIndex,SII.SetStateIndex}, val)
    length(val) == 1 || throw(DimensionMismatch("Cannot set a single index to multiple values $val."))
end
_chk_dimensions(s::SII.ParameterHookWrapper, val) = _chk_dimensions(s.setter, val)
function _chk_dimensions(ms::SII.MultipleSetters, val)
    if length(ms.setters) != length(val)
        msg = if length(val) == 1
            "Cannot set multiple indices ($(length(ms.setters))) to single value $val. Consider using broadcasted .= syntax."
        else
            "Cannot set multiple indices ($(length(ms.setters))) to $(length(val)) values $val."
        end
        throw(DimensionMismatch(msg))
    end
end


####
#### Comparison of states/parameters
####
function Base.:(==)(s1::NWState, s2::NWState)
    s1 === s2 && return true
    SII.variable_symbols(s1) == SII.variable_symbols(s2) || return false
    uflat(s1) == uflat(s2) || return false
    s1.p == s2.p
end
function Base.:(==)(p1::NWParameter, p2::NWParameter)
    p1 === p2 && return true
    SII.parameter_symbols(p1) == SII.parameter_symbols(p2) || return false
    pflat(p1) == pflat(p2)
end
function Base.isapprox(s1::NWState, s2::NWState; kwargs...)
    s1 === s2 && return true
    SII.variable_symbols(s1) == SII.variable_symbols(s2) || return false
    isapprox(uflat(s1), uflat(s2); kwargs...) || return false
    isapprox(s1.p, s2.p; kwargs...)
end
function Base.isapprox(p1::NWParameter, p2::NWParameter; kwargs...)
    p1 === p2 && return true
    SII.parameter_symbols(p1) == SII.parameter_symbols(p2) || return false
    isapprox(pflat(p1), pflat(p2); kwargs...)
end


####
#### enable broadcasted setindex
#### https://discourse.julialang.org/t/broadcasting-setindex-is-a-noobtrap/94700
####
Base.dotview(s::Union{NWParameter, NWState}, idxs...) = view(s, idxs...)

# NWParameter: view
Base.view(s::NWParameter, ::Colon) = s.pflat
function Base.view(p::NWParameter, idx::SymbolicIndex)
    if !SII.is_parameter(p, idx)
        throw(ArgumentError("Cannot create view into NWParameter: Index $idx is not a valid \
            parameter index (states and observables are not allowed)."))
    end
    view(p.pflat, SII.parameter_index(p, idx))
end
function Base.view(p::NWParameter, idxs)
    if !(all(i -> SII.is_parameter(p, i), idxs))
        throw(ArgumentError("Cannot create view into NWParameter: Index $idxs is not a valid \
            parameter index collection (states and observables are not allowed)."))
    end
    view(p.pflat, map(i -> SII.parameter_index(p, i), idxs))
end

# NWState: view
Base.view(s::NWState, ::Colon) = s.uflat
function Base.view(s::NWState, idx::SymbolicIndex)
    if SII.is_variable(s, idx)
        return view(uflat(s), SII.variable_index(s, idx))
    elseif SII.is_parameter(s, idx)
        isnothing(s.p) && throw(ArgumentError("Cannot create view into NWState: No parameter array present."))
        return view(pflat(s), SII.parameter_index(s, idx))
    else
        throw(ArgumentError("Cannot create view into NWState: Index $idx is neither state \
            nor parameter index (observables are not allowed)."))
    end
    view(uflat(s), SII.variable_index(s, idx))
end
function Base.view(s::NWState, idxs)
    if isnothing(s.p) && any(i -> SII.is_parameter(s, i))
        throw(ArgumentError("Cannot create view into NWState: No parameter array present."))
    end
    if all(i -> SII.is_parameter(s, i), idxs)
        _viewidx =  map(i -> SII.parameter_index(s, i), idxs)
        return view(pflat(s), _viewidx)
    elseif all(i -> SII.is_variable(s, i), idxs)
        _viewidx =  map(i -> SII.variable_index(s, i), idxs)
        return view(uflat(s), _viewidx)
    else
        throw(ArgumentError("Cannot create view into NWState: Index $idxs is neither a valid \
            parameter nor state index collection (observables are not allowed)."))
    end
end

###
### Indexing generation
###
struct AllVertices end
struct AllEdges end

generate_indices(inpr, s::SymbolicIndex; kwargs...) = generate_indices(inpr, idxtype(s)(s.compidx), s.subidx; kwargs...)
function generate_indices(inpr, compfilter=nothing, varfilter=nothing; s=true, p=true, out=true, obs=false, in=false,
    return_types=false, just_count=false)
    nw = extract_nw(inpr)
    indices = !just_count ? SymbolicIndex[] : nothing
    types = return_types ? Symbol[] : nothing

    # potentially wrap varfilter
    if varfilter isa Union{Int, Colon, AbstractVector{Int}, NTuple{<:Any, Int}}
        if s
            varfilter = StateIdx(varfilter)
        elseif p && !s
            varfilter = ParamIdx(varfilter)
        else
            throw(ArgumentError("varfilter as Int or collection of Ints only makes sense if either s or p is true."))
        end
    end

    switches = (; s, p, out, obs, in)
    count = _collect_indices!(nw, indices, types, compfilter, varfilter, switches)
    @assert isnothing(indices) || count == length(indices)

    if just_count
        return count
    elseif return_types
        return indices, types
    else
        return indices
    end
end
function _collect_indices!(nw, indices, types, compfilter, varfilter, switches)
    all_verts = VIndex(1:nv(nw))
    all_edges = EIndex(1:ne(nw))
    vidxs = Iterators.filter(compidx -> match_compfilter(nw, compfilter, compidx), all_verts)
    eidxs = Iterators.filter(compidx -> match_compfilter(nw, compfilter, compidx), all_edges)
    count_v = _collect_symbols!(nw, indices, types, vidxs, varfilter, switches)
    count_e = _collect_symbols!(nw, indices, types, eidxs, varfilter, switches)
    count_v + count_e
end
function _collect_symbols!(nw, indices, types, compidxs, varfilter, switches)
    counter = 0
    for cidx in compidxs
        comp = getcomp(nw, cidx)
        if switches.s
            for (i, s) in enumerate(sym(comp))
                if match_varfilter(comp, varfilter, s, StateIdx(i))
                    !isnothing(indices) && push!(indices, idxtype(cidx)(cidx.compidx, s))
                    !isnothing(types) && push!(types, :State)
                    counter += 1
                end
            end
        end
        if switches.p
            for (i, s) in enumerate(psym(comp))
                if match_varfilter(comp, varfilter, s, ParamIdx(i))
                    !isnothing(indices) && push!(indices, idxtype(cidx)(cidx.compidx, s))
                    !isnothing(types) && push!(types, :Parameter)
                    counter += 1
                end
            end
        end
        if switches.out
            for s in outsym_flat(comp)
                switches.s && s ∈ sym(comp) && continue
                if match_varfilter(comp, varfilter, s)
                    !isnothing(indices) && push!(indices, idxtype(cidx)(cidx.compidx, s))
                    !isnothing(types) && push!(types, :Output)
                    counter += 1
                end
            end
        end
        # insym slighly unstable for unknown comp, tried with function barrier but didnt help
        if switches.in && hasinsym(comp)
            for s::Symbol in insym_all(comp)
                if match_varfilter(comp, varfilter, s)
                    !isnothing(indices) && push!(indices, idxtype(cidx)(cidx.compidx, s))
                    !isnothing(types) && push!(types, :Input)
                    counter += 1
                end
            end
        end
        if switches.obs
            for s in obssym(comp)
                if match_varfilter(comp, varfilter, s)
                    !isnothing(indices) && push!(indices, idxtype(cidx)(cidx.compidx, s))
                    !isnothing(types) && push!(types, :Observable)
                    counter += 1
                end
            end
        end
    end
    counter
end

match_compfilter(nw, filter::Nothing, idx) = true
function match_compfilter(nw, filter::Union{AbstractVector, Tuple}, idx)
    any(f -> match_compfilter(nw, f, idx), filter)
end
function match_compfilter(nw, filter::SymbolicIndex{<:Union{AbstractVector, Tuple}}, idx)
    _comptypematch(idx, filter) && any(f -> match_compfilter(nw, f, idx), filter)
end
match_compfilter(nw, filter::AllVertices, idx) = idx isa VIndex
match_compfilter(nw, filter::AllEdges, idx) = idx isa EIndex
match_compfilter(nw, filter::SymbolicIndex{Colon}, idx) = _comptypematch(idx, filter)
match_compfilter(nw, filter::SymbolicIndex{Int}, idx) = _comptypematch(idx, filter) && idx.compidx == filter.compidx
match_compfilter(nw, filter::SymbolicIndex{Symbol}, idx) = _comptypematch(idx, filter) && getcomp(nw, idx).name == filter.compidx
function match_compfilter(nw, filter::SymbolicIndex{<:Union{AbstractString,Regex}}, idx)
    _comptypematch(idx, filter) && occursin(filter.compidx, string(getcomp(nw, idx).name))
end
function match_compfilter(nw, filter::EIndex{<:Pair}, idx)
    _comptypematch(idx, filter) || return false
    src_idx_filter = VIndex(filter.compidx.first)
    dst_idx_filter = VIndex(filter.compidx.second)

    simpleedge = nw.im.edgevec[idx.compidx]
    src_match = match_compfilter(nw, src_idx_filter, VIndex(simpleedge.src))
    dst_match = match_compfilter(nw, dst_idx_filter, VIndex(simpleedge.dst))
    src_match && dst_match
end
_comptypematch(idx, filter) = idxtype(idx) == idxtype(filter)

match_varfilter(comp, filter, sym) = match_varfilter(comp, filter, sym, nothing)
match_varfilter(comp, filter::Nothing, sym, _) = true
match_varfilter(comp, filter::Union{AbstractVector, Tuple}, sym, _) = any(f -> match_varfilter(comp, f, sym), filter)
match_varfilter(comp, filter::Symbol, sym, _) = sym == filter
match_varfilter(comp, filter::Union{AbstractString,Regex}, sym, _) = occursin(filter, string(sym))
function match_varfilter(comp, filter::NumericSubIndex, _, idx)
    idxtype(idx) == idxtype(filter) && (filter.idx == Colon() || idx.idx ∈ filter.idx)
end


# shortcuts
function vidxs(nw, cf=:, vf=nothing; kwargs...)
    generate_indices(nw, VIndex(cf), vf; s=true, p=false, in=false, out=true, obs=true, kwargs...)
end
function eidxs(nw, cf=:, vf=nothing; kwargs...)
    generate_indices(nw, EIndex(cf), vf; s=true, p=false, in=false, out=true, obs=true, kwargs...)
end
function vpidxs(nw, cf=:, vf=nothing; kwargs...)
    generate_indices(nw, VIndex(cf), vf; s=false, p=true, in=false, out=false, obs=false, kwargs...)
end
function epidxs(nw, cf=:, vf=nothing; kwargs...)
    generate_indices(nw, EIndex(cf), vf; s=false, p=true, in=false, out=false, obs=false, kwargs...)
end

# deprecated methods
vidxs(compfilter::Union{Int,AbstractVector{Int},NTuple{<:Any,Int}}, varfilter) = _deprecated_idxs_gen(VIndex, compfilter, varfilter)
eidxs(compfilter::Union{Int,AbstractVector{Int},NTuple{<:Any,Int}}, varfilter) = _deprecated_idxs_gen(EIndex, compfilter, varfilter)
vpidxs(compfilter::Union{Int,AbstractVector{Int},NTuple{<:Any,Int}}, varfilter) = _deprecated_idxs_gen(VPIndex, compfilter, varfilter, VIndex)
epidxs(compfilter::Union{Int,AbstractVector{Int},NTuple{<:Any,Int}}, varfilter) = _deprecated_idxs_gen(EPIndex, compfilter, varfilter, EIndex)

function _deprecated_idxs_gen(constructor, compiter, variter, type=constructor)
    @warn "*idxs(compfilter, varfilter) methods are deprecated. Use *idxs(nw, compfilter, varfilter) instead."
    indices = type[]
    for ci in compfilter
        for si in varfilter
            push!(indices, constructor(ci, si))
        end
    end
    indices
end


####
#### Indexing proxys
####
function Base.getproperty(s::Union{NWParameter, NWState}, sym::Symbol)
    if sym === :v
        return FilteringProxy(s).v
    elseif sym === :e
        return FilteringProxy(s).e
    else
        return getfield(s, sym)
    end
end

abstract type IndexingProxy{S} end

struct FilteringProxy{S,CF,VF} <: IndexingProxy{S}
    data::S # the underlying state/parameter object
    compfilter::CF
    varfilter::VF
    states::Bool
    parameters::Bool
    outputs::Bool
    inputs::Bool
    observables::Bool
end

function (f::FilteringProxy)(args...; kwargs...)
    if !isempty(args)
        f = refine_filter(f, args...)
    end
    if !isempty(kwargs)
        f = FilteringProxy(f; kwargs...)
    end
    return f
end
Base.keys(f::FilteringProxy) = generate_indices(f)
Base.values(f::FilteringProxy) = f.data[generate_indices(f)]

# "copy" constructor
function FilteringProxy(
    f::FilteringProxy;
    compfilter=f.compfilter, varfilter=f.varfilter,
    s=f.states, p=f.parameters, out=f.outputs,
    in=f.inputs, obs=f.observables
)
    FilteringProxy(f.data, compfilter, varfilter, s, p, out, in, obs)
end
FilteringProxy(s::NWState) = FilteringProxy(s, nothing, nothing, true, !isnothing(s.p), false, false, false)
FilteringProxy(s::NWParameter) = FilteringProxy(s, nothing, nothing, false, true, false, false, false)
function generate_indices(f::FilteringProxy; kwargs...)
    generate_indices(f.data, f.compfilter, f.varfilter;
        s=f.states,
        p=f.parameters,
        out=f.outputs,
        in=f.inputs,
        obs=f.observables,
        kwargs...)
end
function resolve_to_index(f::FilteringProxy)
    idxs = generate_indices(f)
    if length(idxs) == 1
        only(idxs)
    else
        idxs
    end
end
resolve_to_value(f::FilteringProxy) = f.data[resolve_to_index(f)]
# shortcut route if the "filter" is just a concrete index
function resolve_to_value(f::FilteringProxy{<:Any, <:SymbolicIndex{<:CONCRETE_COMPIDX, Nothing}, <:CONCRETE_SUBIDX})
    idx = idxtype(f.compfilter)(f.compfilter.compidx, f.varfilter)
    f.data[idx]
end

function is_fully_refined(f::FilteringProxy)
    f.compfilter isa SymbolicIndex && !isnothing(f.varfilter)
end

# acessing the underlying data directly with concrete index
function Base.getindex(f::FilteringProxy, idx::SymbolicIndex)
    if !isnothing(idx.subidx)
        getindex(f.data, idx)
    else # if idx.subidx is nothing its just a refinement
        new_filter = refine_filter(f, idx)
        is_fully_refined(new_filter) ? resolve_to_value(new_filter) : new_filter
    end
end
function Base.getindex(f::FilteringProxy, ref)
    new_filter = refine_filter(f, ref)
    is_fully_refined(new_filter) ? resolve_to_value(new_filter) : new_filter
end
Base.getindex(f::FilteringProxy, a, b) = refine_filter(f, a)[b]
# vector valued things
function Base.getindex(f::FilteringProxy, idxs::Union{AbstractVector,Tuple})
    if all(i -> i isa SymbolicIndex && !isnothing(i.subidx), idxs)
        return getindex(f.data, idxs)
    else
        new_filter = refine_filter(f, idxs)
        return is_fully_refined(new_filter) ? resolve_to_value(new_filter) : new_filter
    end
end

function Base.getproperty(f::FilteringProxy, sym::Symbol)
    if sym === :v
        if f.compfilter == nothing
            FilteringProxy(f; compfilter=AllVertices())
        else
            throw(ArgumentError("Cannot apply filter \".v\": already filtering for $(f.compfilter)."))
        end
    elseif sym === :e
        if f.compfilter == nothing
            FilteringProxy(f; compfilter=AllEdges())
        else
            throw(ArgumentError("Cannot apply filter \".e\": already filtering for $(f.compfilter)."))
        end
    elseif sym === :p
        FilteringProxy(f; p=true, s=false, out=false, in=false, obs=false)
    elseif sym === :s
        FilteringProxy(f; p=false, s=true, out=false, in=false, obs=false)
    elseif sym === :out
        FilteringProxy(f; p=false, s=false, out=true, in=false, obs=false)
    elseif sym === :in
        FilteringProxy(f; p=false, s=false, out=false, in=true, obs=false)
    elseif sym === :obs
        FilteringProxy(f; p=false, s=false, out=false, in=false, obs=true)
    elseif sym === :all
        FilteringProxy(f; p=true, s=true, out=true, in=true, obs=true)
    else
        getfield(f, sym)
    end
end

refine_filter(f::FilteringProxy, a, b) = refine_filter(refine_filter(f, a), b)

function refine_filter(f::FilteringProxy{<:Any, <:AllVertices}, idxs::Union{AbstractVector,Tuple})
    return FilteringProxy(f, compfilter=VIndex(idxs))
end
function refine_filter(f::FilteringProxy{<:Any, <:AllEdges}, idxs::Union{AbstractVector,Tuple})
    return FilteringProxy(f, compfilter=EIndex(idxs))
end
function refine_filter(f::FilteringProxy{<:Any, <:SymbolicIndex}, idxs::Union{AbstractVector,Tuple})
    return FilteringProxy(f, varfilter=idxs)
end

####
#### refinement by SymbolicIndex(something, nothing) is allways compfilter update
#### refinement by SymolicIndex(something, something) is allways comp + varfilter update
####
function refine_filter(f::FilteringProxy{<:Any}, s::SymbolicIndex{<:Any, <:Any})
    isnothing(s.compidx) && error("Index $s has component nothing, which is not a valid refinement")
    if isnothing(s.subidx)
        FilteringProxy(f; compfilter=s)
    else
        FilteringProxy(f; compfilter=idxtype(s)(s.compidx), varfilter=s.subidx)
    end
end
####
#### if unspecific compfilter is set (i.e. AllVertices/AllEdges by .v .e property)
#### we refine the comfilter by that type
####
function refine_filter(f::FilteringProxy{<:Any, AllVertices}, idx::Union{Colon,Int,Symbol,Pair,String,Regex})
    FilteringProxy(f; compfilter=VIndex(idx))
end
function refine_filter(f::FilteringProxy{<:Any, AllEdges}, idx::Union{Colon,Int,Symbol,Pair,String,Regex})
    FilteringProxy(f; compfilter=EIndex(idx))
end
####
#### From a specific compfilter we can only refine the varfilter
####
function refine_filter(f::FilteringProxy{<:Any, <:SymbolicIndex}, varfilter)
    if varfilter isa Union{StateIdx, ParamIdx, Symbol, String, Regex,
                           Colon, Int, AbstractVector{Int}, NTuple{<:Any,<:Int}}
        FilteringProxy(f; varfilter=varfilter)
    else
        msg = "Index $varfilter is not a valid variable refinement on $(f.compfilter)."
        throw(ArgumentError(msg))
    end
end

function Base.setindex!(f::FilteringProxy, val, idx...)
    refined_f = refine_filter(f, idx...)
    allindices = resolve_to_index(refined_f)
    f.data[allindices] = val
end

#### enable broadcasted setindex
#### https://discourse.julialang.org/t/broadcasting-setindex-is-a-noobtrap/94700
function Base.dotview(f::FilteringProxy, idxs...)
    refined_f = refine_filter(f, idxs...)
    allindices = resolve_to_index(refined_f)
    view(f.data, allindices)
end

# methods below are used to allow for s.v[1].p .= 1.0 style broadcasting
# which his not covered by "dotview" becaus it is a full proxy object
Base.size(f::FilteringProxy) = length(resolve_to_index(f))
Base.ndims(::Type{<:FilteringProxy}) = 1
function Base.copyto!(f::FilteringProxy, bc::Broadcast.Broadcasted{<:Base.Broadcast.DefaultArrayStyle})
    allindices = resolve_to_index(f)
    v = view(f.data, allindices)
    Base.copyto!(v, bc)
    return f
end



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
extract_nw(p::IndexingProxy) = extract_nw(p.data)
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
function Base.getindex(nw::Network, i::EIndex{<:CONCRETE_COMPIDX, Nothing})
    return getcomp(nw, i)
end
function Base.getindex(nw::Network, i::VIndex{<:Union{Symbol,Int}, Nothing})
    return getcomp(nw, i)
end
function Base.getindex(nw::Network, i::EIndex{Colon, Nothing})
    return copy(nw.im.edgem)
end
function Base.getindex(nw::Network, i::VIndex{Colon, Nothing})
    return copy(nw.im.vertexm)
end
function Base.getindex(nw::Network, collection::Union{AbstractArray,Tuple})
    getindex.(Ref(nw), collection)
end
