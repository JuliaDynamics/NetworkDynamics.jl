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
Base.getindex(s::NWState, idx) = SII.getu(s, idx)(s)

# NWState: setindex!
function Base.setindex!(s::NWState, val, idx)
    setter = SII.setu(s, idx)
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
#### Indexing proxys
####
abstract type IndexingProxy{S} end
struct VProxy{S} <: IndexingProxy{S}
    s::S
end
Base.getindex(p::VProxy, comp, state) = getindex(p.s, VIndex(comp, _wrap_proxy_subidx(p, state)))
Base.getindex(p::VProxy, ::Colon, state) = getindex(p, 1:nv(extract_nw(p)), _wrap_proxy_subidx(p, state))
Base.setindex!(p::VProxy, val, comp, state) = setindex!(p.s, val, VIndex(comp, _wrap_proxy_subidx(p, state)))

struct EProxy{S} <: IndexingProxy{S}
    s::S
end
Base.getindex(p::EProxy, comp, state) = getindex(p.s, EIndex(comp, _wrap_proxy_subidx(p, state)))
Base.getindex(p::EProxy, ::Colon, state) = getindex(p, 1:ne(extract_nw(p)), _wrap_proxy_subidx(p, state))
Base.setindex!(p::EProxy, val, comp, state) = setindex!(p.s, val, EIndex(comp, _wrap_proxy_subidx(p, state)))

_wrap_proxy_subidx(::IndexingProxy{<:NWState}, sub) = wrap_sidx(sub)
_wrap_proxy_subidx(::IndexingProxy{<:NWParameter}, sub) = wrap_pidx(sub)

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
    if !SII.is_parameter(p, idx)
        throw(ArgumentError("Index $idx is not a valid parameter index."))
    end
    view(p.pflat, SII.parameter_index(p, idx))
end
function Base.view(p::NWParameter, idxs)
    if !(all(i -> SII.is_parameter(p, i), idxs))
        throw(ArgumentError("Index $idxs is not a valid parameter index collection."))
    end
    view(p.pflat, map(i -> SII.parameter_index(p, i), idxs))
end

# NWState: view
Base.view(s::NWState, ::Colon) = s.uflat
function Base.view(s::NWState, idx::SymbolicIndex)
    if SII.is_variable(s, idx)
        return view(uflat(s), SII.variable_index(s, idx))
    elseif SII.is_parameter(s, idx)
        return view(pflat(s), SII.parameter_index(s, idx))
    else
        throw(ArgumentError("Index $idx is neither state nor parameter index."))
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
        throw(ArgumentError("Index $idxs is neither a valid parameter nor state index collection."))
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
    res = SymbolicIndex[]
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
_make_cidx_iterable(::Union{Type{<:VIndex}, <:typeof(VPIndex)}, inpr, ::Colon) = 1:nv(extract_nw(inpr))
_make_cidx_iterable(::Union{Type{<:EIndex}, <:typeof(EPIndex)}, inpr, ::Colon) = 1:ne(extract_nw(inpr))
function _make_cidx_iterable(IT, inpr, s::Symbol)
    names = getproperty.(_get_components(IT, inpr), :name)
    findall(isequal(s), names)
end
function _make_cidx_iterable(IT, inpr, s::Union{AbstractString,AbstractPattern})
    names = getproperty.(_get_components(IT, inpr), :name)
    findall(sym -> contains(string(sym), s), names)
end

_make_sidx_iterable(IT, inpr, cidx, idx) = _make_iterabel(idx)
function _make_sidx_iterable(IT::Type{<:SymbolicIndex}, inpr, cidx, ::Colon)
    comp = _get_components(IT, inpr)[cidx]
    Iterators.flatten((sym(comp), obssym_all(comp)))
end
function _make_sidx_iterable(IT::Union{typeof(VPIndex), typeof(EPIndex)}, inpr, cidx, ::Colon)
    psym(_get_components(IT, inpr)[cidx])
end
function _make_sidx_iterable(IT::Type{<:SymbolicIndex}, inpr, cidx, s::Union{AbstractString,AbstractPattern})
    syms = _make_sidx_iterable(IT, inpr, cidx, :) # get all possible
    Iterators.filter(sym -> contains(string(sym), s), syms)
end
function _make_sidx_iterable(IT::Union{typeof(VPIndex), typeof(EPIndex)}, inpr, cidx, s::Union{AbstractString,AbstractPattern})
    syms = _make_sidx_iterable(IT, inpr, cidx, :) # get all possible
    filter(sym -> contains(string(sym), s), syms)
end

_get_components(::Union{Type{<:VIndex}, <:typeof(VPIndex)}, inpr) = extract_nw(inpr).im.vertexm
_get_components(::Union{Type{<:EIndex}, <:typeof(EPIndex)}, inpr) = extract_nw(inpr).im.edgem


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
