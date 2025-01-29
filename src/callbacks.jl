"""
    abstract type ComponentCallback end

Abstract type for a component based callback. A component callback
bundles a [`ComponentCondition`](@ref) as well as a [`ComponentAffect`](@ref)
which can be then tied to a component model using [`add_callback!`](@ref) or
[`set_callback!`](@ref).

On a Network level, you can automaticially create network wide `CallbackSet`s using
[`get_callbacks`](@ref).

See
[`ContinousComponentCallback`](@ref) and [`VectorContinousComponentCallback`](@ref) for concrete
implemenations of this abstract type.
"""
abstract type ComponentCallback end

"""
    ComponentCondition(f::Function, sym, psym)

Creates a callback condition for a [`ComponentCallback`].
- `f`: The condition function. Must be a function of the form `out=f(u, p, t)`
  when used for [`ContinousComponentCallback`](@ref) or
  [`DiscreteComponentCallback`](@ref) and `f!(out, u, p, t)` when used for
  [`VectorContinousComponentCallback`](@ref).
  - Arguments of `f`
    - `u`: The current value of the selecte `sym` states, provided as a [`SymbolicView`](@ref) object.
    - `p`: The current value of the selected `psym` parameters.
    - `t`: The current simulation time.
- `sym`: A vector or tuple of symbols, which represent **states** (including
  inputs, outputs, observed) of the component model. Determines, which states will
  be available thorugh parameter `u` in the callback condition function `f`.
- `psym`: A vector or tuple of symbols, which represetn **parameters** of the component mode.
  Determines, which parameters will be available in the condition function `f`

# Example
Consider a component model with states `[:u1, :u2]`, inputs `[:i]`, outputs
`[:o]` and parameters `[:p1, :p2]`.

    ComponentCondition([:u1, :o], [:p1]) do (u, p, t)
        # access states symbolicially or via int index
        u[:u1] == u[1]
        u[:o] == u[2]
        p[:p1] == p[1]
        # the states/prameters `:u2`, `:i` and `:p2` are not available as
        # they are not listed in the `sym` and `psym` arguments.
    end
"""
struct ComponentCondition{C,DIM,PDIM}
    f::C
    sym::NTuple{DIM,Symbol}
    psym::NTuple{PDIM,Symbol}
    function ComponentCondition(f, sym, psym)
        new{typeof(f), length(sym), length(psym)}(f, Tuple(sym), Tuple(psym))
    end
end

"""
    ComponentAffect(f::Function, sym, psym)

Creates a callback condition for a [`ComponentCallback`].
- `f`: The affect function. Must be a function of the form `f(u, p, [event_idx], ctx)` where `event_idx`
  is only available in [`VectorContinousComponentCallback`](@ref).
  - Arguments of `f`
    - `u`: The current (mutable) value of the selected `sym` states, provided as a [`SymbolicView`](@ref) object.
    - `p`: The current (mutalbe) value of the selected `psym` parameters.
    - `event_idx`: The current event index, i.e. which `out` element triggerd in case of [`VectorContinousComponentCallback`](@ref).
    - `ctx::NamedTuple` a named tuple with context variables.
       - `ctx.model`: a referenc to the ocmponent model
       - `ctx.vidx`/`ctx.eidx`: The index of the vertex/edge model.
       - `ctx.src`/`ctx.dst`: src and dst indices (only for edge models).
       - `ctx.integrator`: The integrator object.
       - `ctx.t=ctx.integrator.t`: The current simulation time.
- `sym`: A vector or tuple of symbols, which represent **states** (**excluding**
  inputs, outputs, observed) of the component model. Determines, which states will
  be available thorugh parameter `u` in the callback condition function `f`.
- `psym`: A vector or tuple of symbols, which represetn **parameters** of the component mode.
  Determines, which parameters will be available in the condition function `f`

# Example
Consider a component model with states `[:u1, :u2]`, inputs `[:i]`, outputs
`[:o]` and parameters `[:p1, :p2]`.

    ComponentAffect([:u1, :o], [:p1]) do (u, p, ctx)
        u[:u1] = 0 # change the state
        p[:p1] = 1 # change the parameter
        @info "Changed :u1 and :p1 on vertex \$(ctx.vidx)" # access context
    end
"""
struct ComponentAffect{A,DIM,PDIM}
    f::A
    sym::NTuple{DIM,Symbol}
    psym::NTuple{PDIM,Symbol}
    function ComponentAffect(f, sym, psym)
        new{typeof(f), length(sym), length(psym)}(f, Tuple(sym), Tuple(psym))
    end
end

"""
    ContinousComponentCallback(condition, affect; kwargs...)

Connect a [`ComponentCondition`](@ref) and a [`ComponentAffect`)[@ref] to a
continous callback which can be attached to a component model using
[`add_callback!`](@ref) or [`set_callback!`](@ref).

The `kwargs` will be forwarded to the `VectorContinuousCallback` when the component based
callbacks are collected for the whole network using `get_callbacks`.
[`DiffEq.jl docs`](https://docs.sciml.ai/DiffEqDocs/stable/features/callback_functions/)
for available options.
"""
struct ContinousComponentCallback{C,A,CDIM,CPDIM,ADIM,APDIM} <: ComponentCallback
    condition::ComponentCondition{C,CDIM,CPDIM}
    affect::ComponentAffect{A,ADIM,APDIM}
    kwargs::NamedTuple
end
function ContinousComponentCallback(condition, affect; kwargs...)
    if haskey(kwargs, :affect_neg!)
        throw(ArgumentError("affect_neg! not supported yet. Please raise issue."))
    end
    ContinousComponentCallback(condition, affect, NamedTuple(kwargs))
end

"""
    VectorContinousComponentCallback(condition, affect, len; kwargs...)

Connect a [`ComponentCondition`](@ref) and a [`ComponentAffect`)[@ref] to a
continous callback which can be attached to a component model using
[`add_callback!`](@ref) or [`set_callback!`](@ref). This vector version allows
for `condions` which have `len` output dimensions.
The `affect` will be triggered with the additional `event_idx` argument to know in which
dimension the zerocrossing was detected.

The `kwargs` will be forwarded to the `VectorContinuousCallback` when the component based
callbacks are collected for the whole network using [`get_callbacks(::Network)`](@ref).
[`DiffEq.jl docs`](https://docs.sciml.ai/DiffEqDocs/stable/features/callback_functions/)
for available options.
"""
struct VectorContinousComponentCallback{C,A,CDIM,CPDIM,ADIM,APDIM} <: ComponentCallback
    condition::ComponentCondition{C,CDIM,CPDIM}
    affect::ComponentAffect{A,ADIM,APDIM}
    len::Int
    kwargs::NamedTuple
end
function VectorContinousComponentCallback(condition, affect, len; kwargs...)
    if haskey(kwargs, :affect_neg!)
        throw(ArgumentError("affect_neg! not supported yet. Please raise issue."))
    end
    VectorContinousComponentCallback(condition, affect, len, NamedTuple(kwargs))
end

"""
    DiscreteComponentCallback(condition, affect; kwargs...)

Connect a [`ComponentCondition`](@ref) and a [`ComponentAffect`)[@ref] to a
discrete callback which can be attached to a component model using
[`add_callback!`](@ref) or [`set_callback!`](@ref).

Note that the `condition` function returns a boolean value, as the discrete
callback perform no rootfinding.

The `kwargs` will be forwarded to the `DiscreteCallback` when the component based
callbacks are collected for the whole network using [`get_callbacks(::Network)`](@ref).
[`DiffEq.jl docs`](https://docs.sciml.ai/DiffEqDocs/stable/features/callback_functions/)
for available options.
"""
struct DiscreteComponentCallback{C<:ComponentCondition,A<:ComponentAffect} <: ComponentCallback
    condition::C
    affect::A
    kwargs::NamedTuple
end
function DiscreteComponentCallback(condition, affect; kwargs...)
    DiscreteComponentCallback(condition, affect, NamedTuple(kwargs))
end

"""
    PresetTimeComponentCallback(ts, affect; kwargs...)

Tirgger a [`ComponentAffect`](@ref) at given timesteps `ts` in discrete
callback, which can be attached to a component model using
[`add_callback!`](@ref) or [`set_callback!`](@ref).

The `kwargs` will be forwarded to the [`PresetTimeCallback`](@extref DiffEqCallbacks.PresetTimeCallback)
when the component based callbacks are collected for the whole network using
[`get_callbacks(::Network)`](@ref).

The `PresetTimeCallback` will take care of adding the timesteps to the solver, ensuring to
exactly trigger at the correct times.
"""
struct PresetTimeComponentCallback{T,A} <: ComponentCallback
    ts::T
    affect::A
    kwargs::NamedTuple
end
function PresetTimeComponentCallback(ts, affect; kwargs...)
    PresetTimeComponentCallback(ts, affect, NamedTuple(kwargs))
end

"""
    get_callbacks(nw::Network)::CallbackSet

Returns a `CallbackSet` composed of all the "component-based" callbacks in the metadata of the
Network components.
"""
function get_callbacks(nw::Network)
    aliased_changed(nw; warn=true)
    cbbs = wrap_component_callbacks(nw)
    if isempty(cbbs)
        return nothing
    elseif length(cbbs) == 1
        return to_callback(only(cbbs))
    else
        CallbackSet(to_callback.(cbbs)...)
    end
end

####
#### identifying callbacks which can be combined into batches
####
function wrap_component_callbacks(nw)
    components = SymbolicIndex[]
    callbacks = ComponentCallback[]
    for (i, v) in pairs(nw.im.vertexm)
        has_callback(v) || continue
        for cb in get_callbacks(v)
            push!(components, VIndex(i, nothing))
            push!(callbacks, cb)
        end
    end
    for (i, v) in pairs(nw.im.edgem)
        has_callback(v) || continue
        for cb in get_callbacks(v)
            push!(components, EIndex(i, nothing))
            push!(callbacks, cb)
        end
    end
    # group the callbacks such that they are in groups which are "batchequal"
    # batchequal groups can be wrapped into a single callback
    idx_per_type = _find_identical(callbacks, 1:length(components))
    batches = []
    for typeidx in idx_per_type
        batchcomps = components[typeidx]
        batchcbs = callbacks[typeidx]
        if first(batchcbs) isa Union{ContinousComponentCallback, VectorContinousComponentCallback}
            cb = ContinousCallbackWrapper(nw, batchcomps, batchcbs)
        elseif only(batchcbs) isa Union{DiscreteComponentCallback, PresetTimeComponentCallback}
            cb = DiscreteCallbackWrapper(nw, only(batchcomps), only(batchcbs))
        else
            error("Unknown callback type, should never be reached. Please report this issue.")
        end
        push!(batches, cb)
    end
    return batches
end
function batchequal(a::ContinousComponentCallback, b::ContinousComponentCallback)
    batchequal(a.condition, b.condition) || return false
    batchequal(a.kwargs, b.kwargs)       || return false
    return true
end
function batchequal(a::VectorContinousComponentCallback, b::VectorContinousComponentCallback)
    batchequal(a.condition, b.condition) || return false
    batchequal(a.kwargs, b.kwargs)       || return false
    a.len == b.len                       || return false
    return true
end
function batchequal(a::T, b::T) where {T <: Union{ComponentCondition, ComponentAffect}}
    typeof(a) == typeof(b)
end
function batchequal(a::NamedTuple, b::NamedTuple)
    length(a) == length(b) || return false
    for (k, v) in a
        haskey(b, k) || return false
        v == b[k] || return false
    end
    return true
end

# a callback wrapper is a container, which wraps a network, a component callback
# and a component index. It is used for bookkeeping to know to which component
# each callback belongs to
abstract type CallbackWrapper end

####
#### wrapping of continous callbacks
####
struct ContinousCallbackWrapper{T<:ComponentCallback,C,ST<:SymbolicIndex} <: CallbackWrapper
    nw::Network
    components::Vector{ST}
    callbacks::Vector{T}
    sublen::Int # length of each callback
    condition::C
end
function ContinousCallbackWrapper(nw, components, callbacks)
    if !isconcretetype(eltype(components))
        components = [c for c in components]
    end
    if !isconcretetype(eltype(callbacks))
        callbacks = [cb for cb in callbacks]
    end
    sublen = eltype(callbacks) <: ContinousComponentCallback ? 1 : first(callbacks).len
    condition = first(callbacks).condition.f
    ContinousCallbackWrapper(nw, components, callbacks, sublen, condition)
end

Base.length(ccw::ContinousCallbackWrapper) = length(ccw.callbacks)
cbtype(ccw::ContinousCallbackWrapper{T}) where {T} = T

condition_dim(ccw)  = first(ccw.callbacks).condition.sym  |> length
condition_pdim(ccw) = first(ccw.callbacks).condition.psym |> length
affect_dim(ccw,i)  = ccw.callbacks[i].affect.sym  |> length
affect_pdim(ccw,i) = ccw.callbacks[i].affect.psym |> length

condition_urange(ccw, i) = (1 + (i-1)*condition_dim(ccw))  : i*condition_dim(ccw)
condition_prange(ccw, i) = (1 + (i-1)*condition_pdim(ccw)) : i*condition_pdim(ccw)
affect_urange(ccw, i) = (1 + (i-1)*affect_dim(ccw,i) ) : i*affect_dim(ccw,i)
affect_prange(ccw, i) = (1 + (i-1)*affect_pdim(ccw,i)) : i*affect_pdim(ccw,i)

condition_outrange(ccw, i) = (1 + (i-1)*ccw.sublen) : i*ccw.sublen

cbidx_from_outidx(ccw, outidx) = div(outidx-1, ccw.sublen) + 1
subidx_from_outidx(ccw, outidx) = mod(outidx, 1:ccw.sublen)

function collect_c_or_a_indices(ccw::ContinousCallbackWrapper, c_or_a, u_or_p)
    sidxs = SymbolicIndex[]
    for (component, cb) in zip(ccw.components, ccw.callbacks)
        syms = getproperty(getproperty(cb, c_or_a), u_or_p)
        symidxtype = if component isa VIndex
            u_or_p == :sym ? VIndex : VPIndex
        else
            u_or_p ==:sym ? EIndex : EPIndex
        end
        sidx = collect(symidxtype(component.compidx, syms))
        append!(sidxs, sidx)
    end
    sidxs
end

# generate VectorContinuousCallback from a ContinousCallbackWrapper
function to_callback(ccw::ContinousCallbackWrapper)
    kwargs = first(ccw.callbacks).kwargs
    cond = _batch_condition(ccw)
    affect = _batch_affect(ccw)
    len = ccw.sublen * length(ccw.callbacks)
    VectorContinuousCallback(cond, affect, len; kwargs...)
end
function _batch_condition(ccw::ContinousCallbackWrapper)
    usymidxs = collect_c_or_a_indices(ccw, :condition, :sym)
    psymidxs = collect_c_or_a_indices(ccw, :condition, :psym)
    ucache = DiffCache(zeros(length(usymidxs)), 12)

    obsf = SII.observed(ccw.nw, usymidxs)
    pidxs = SII.parameter_index.(Ref(ccw.nw), psymidxs)

    (out, u, t, integrator) -> begin
        us = PreallocationTools.get_tmp(ucache, u)
        obsf(u, integrator.p, t, us) # fills us inplace

        for i in 1:length(ccw)
            # symbolic view into u
            uv = view(us, condition_urange(ccw, i))
            _u = SymbolicView(uv, ccw.callbacks[i].condition.sym)

            # symbolic view into p
            pidxsv = view(pidxs, condition_prange(ccw, i))
            pv = view(integrator.p, pidxsv)
            _p = SymbolicView(pv, ccw.callbacks[i].condition.psym)

            if cbtype(ccw) <: ContinousComponentCallback
                oidx = only(condition_outrange(ccw, i))
                out[oidx] = ccw.condition(_u, _p, t)
            elseif cbtype(ccw) <: VectorContinousComponentCallback
                @views _out = out[condition_outrange(ccw, i)]
                ccw.condition(_out, _u, _p, t)
            else
                error()
            end
        end
        nothing
    end
end
function _batch_affect(ccw::ContinousCallbackWrapper)
    usymidxs = collect_c_or_a_indices(ccw, :affect, :sym)
    psymidxs = collect_c_or_a_indices(ccw, :affect, :psym)

    uidxs = SII.variable_index.(Ref(ccw.nw), usymidxs)
    pidxs = SII.parameter_index.(Ref(ccw.nw), psymidxs)

    (integrator, outidx) -> begin
        i = cbidx_from_outidx(ccw, outidx)

        uidxsv = view(uidxs, affect_urange(ccw, i))
        uv = view(integrator.u, uidxsv)
        _u = SymbolicView(uv, ccw.callbacks[i].affect.sym)

        pidxsv = view(pidxs, affect_prange(ccw, i))
        pv = view(integrator.p, pidxsv)
        _p = SymbolicView(pv, ccw.callbacks[i].affect.psym)

        ctx = get_ctx(integrator, ccw.components[i])

        uhash = hash(uv)
        phash = hash(pv)
        if cbtype(ccw) <: ContinousComponentCallback
            ccw.callbacks[i].affect.f(_u, _p, ctx)
        elseif cbtype(ccw) <: VectorContinousComponentCallback
            num = subidx_from_outidx(ccw, outidx)
            ccw.callbacks[i].affect.f(_u, _p, num, ctx)
        else
            error()
        end
        pchanged = hash(pv) != phash
        uchanged = hash(uv) != uhash

        (pchanged || uchanged) && SciMLBase.auto_dt_reset!(integrator)
        pchanged && save_parameters!(integrator)
    end
end

####
#### wrapping of discrete callbacks
####
struct DiscreteCallbackWrapper{N,ST,T} <: CallbackWrapper
    nw::N
    component::ST
    callback::T
    function DiscreteCallbackWrapper(nw, component, callback)
        @assert nw isa Network
        @assert component isa SymbolicIndex
        @assert callback isa Union{DiscreteComponentCallback, PresetTimeComponentCallback}
        new{typeof(nw),typeof(component),typeof(callback)}(nw, component, callback)
    end
end

# generate a DiscreteCallback from a DiscreteCallbackWrapper
function to_callback(dcw::DiscreteCallbackWrapper)
    kwargs = dcw.callback.kwargs
    cond = _batch_condition(dcw)
    affect = _batch_affect(dcw)
    DiscreteCallback(cond, affect; kwargs...)
end
# generate a PresetTimeCallback from a DiscreteCallbackWrapper
function to_callback(dcw::DiscreteCallbackWrapper{<:Any,<:Any,<:PresetTimeComponentCallback})
    kwargs = dcw.callback.kwargs
    affect = _batch_affect(dcw)
    ts = dcw.callback.ts
    DiffEqCallbacks.PresetTimeCallback(ts, affect; kwargs...)
end
function _batch_condition(dcw::DiscreteCallbackWrapper)
    uidxtype = dcw.component isa EIndex ? EIndex : VIndex
    pidxtype = dcw.component isa EIndex ? EPIndex : VPIndex
    usymidxs = uidxtype(dcw.component.compidx, dcw.callback.condition.sym)
    psymidxs = pidxtype(dcw.component.compidx, dcw.callback.condition.psym)
    ucache = DiffCache(zeros(length(usymidxs)), 12)

    obsf = SII.observed(dcw.nw, usymidxs)
    pidxs = SII.parameter_index.(Ref(dcw.nw), psymidxs)

    (u, t, integrator) -> begin
        us = PreallocationTools.get_tmp(ucache, u)
        obsf(u, integrator.p, t, us) # fills us inplace
        _u = SymbolicView(u, dcw.callback.condition.sym)
        pv = view(integrator.p, pidxs)
        _p = SymbolicView(pv, dcw.callback.condition.psym)
        dcw.callback.condition.f(_u, _p, t)
    end
end
function _batch_affect(dcw::DiscreteCallbackWrapper)
    uidxtype = dcw.component isa EIndex ? EIndex : VIndex
    pidxtype = dcw.component isa EIndex ? EPIndex : VPIndex
    usymidxs = uidxtype(dcw.component.compidx, dcw.callback.affect.sym)
    psymidxs = pidxtype(dcw.component.compidx, dcw.callback.affect.psym)

    uidxs = SII.variable_index.(Ref(dcw.nw), usymidxs)
    pidxs = SII.parameter_index.(Ref(dcw.nw), psymidxs)

    (integrator) -> begin
        uv = view(integrator.u, uidxs)
        _u = SymbolicView(uv, dcw.callback.affect.sym)
        pv = view(integrator.p, pidxs)
        _p = SymbolicView(pv, dcw.callback.affect.psym)
        ctx = get_ctx(integrator, dcw.component)

        uhash = hash(uv)
        phash = hash(pv)
        dcw.callback.affect.f(_u, _p, ctx)
        pchanged = hash(pv) != phash
        uchanged = hash(uv) != uhash

        (pchanged || uchanged) && SciMLBase.auto_dt_reset!(integrator)
        pchanged && save_parameters!(integrator)
    end
end


####
#### generate the context for the callback effects
####
function get_ctx(integrator, sym::VIndex)
    nw = extract_nw(integrator)
    idx = sym.compidx
    (; integrator, t=integrator.t, model=nw[sym], vidx=idx)
end
function get_ctx(integrator, sym::EIndex)
    nw = extract_nw(integrator)
    idx = sym.compidx
    edge = nw.im.edgevec[idx]
    (; integrator, t=integrator.t, model=nw[sym], eidx=idx, src=edge.src, dst=edge.dst)
end

####
#### SymbolicView helper type
####
"""
    SymbolicView{N,VT} <: AbstractVetor{VT}

Is a (smallish) fixed size vector type with named dimensions.
Its main purpose is to allow named acces to variables in
[`ComponentCondition`](@ref) and [`ComponentAffect`](@ref) functions.

I.e. when the `ComponentAffect` declared `sym=[:x, :y]`, you can
acces `u[:x]` and `u[:y]` inside the condition function.
"""
struct SymbolicView{N,VT} <: AbstractVector{VT}
    v::VT
    syms::NTuple{N,Symbol}
end
Base.IteratorSize(::Type{SymbolicView}) = IteratorSize(x.v)
Base.IteratorEltype(::Type{SymbolicView}) = IteratorEltype(x.v)
Base.eltype(::Type{SymbolicView}) = eltype(x.v)
Base.size(x::SymbolicView) = size(x.v)
Base.firstindex(x::SymbolicView) = firstindex(x.v)
Base.lastindex(x::SymbolicView) = lastindex(x.v)
Base.iterate(x::SymbolicView) = iterate(x.v)
Base.iterate(x::SymbolicView, state) = iterate(x.v, state)
Base.length(x::SymbolicView{N}) where {N} = N
Base.IndexStyle(::Type{SymbolicView}) = IndexLinear()
Base.getindex(x::SymbolicView, index) = x.v[_sym_to_int(x, index)]
Base.setindex!(x::SymbolicView, value, index) = x.v[_sym_to_int(x, index)] = value
function _sym_to_int(x::SymbolicView, sym::Symbol)
    idx = findfirst(isequal(sym), x.syms)
    isnothing(idx) && throw(ArgumentError("SymbolError: try to access SymbolicView($(x.syms)) with symbol $sym"))
    idx
end
_sym_to_int(x::SymbolicView, idx::Int) = idx
_sym_to_int(x::SymbolicView, idx) = _sym_to_int.(Ref(x), idx)

####
#### Internal function to check cb compat when added as metadata
####
assert_cb_compat(comp::ComponentModel, t::Tuple) = assert_cb_compat.(Ref(comp), t)
function assert_cb_compat(comp::ComponentModel, cb)
    insym = hasinsym(comp) ? insym_all(comp) : []
    all_obssym = Set(sym(comp)) ∪ Set(comp.obssym) ∪ insym ∪ outsym_flat(comp)
    pcond = s -> s in comp.psym
    ucond_cond = s -> s in all_obssym
    ucond_affect = s -> s in comp.sym

    if !(cb isa PresetTimeComponentCallback)
        if !(all(ucond_cond, cb.condition.sym))
            invalid = filter(!ucond_cond, cb.condition.sym)
            throw(ArgumentError("All u symbols in the callback condition must be observed or variable. Found invalid $invalid !⊆ $all_obssym."))
        end
        if !(all(pcond, cb.condition.psym))
            invalid = filter(!pcond, cb.condition.psym)
            throw(ArgumentError("All p symbols in the callback condition must be parameters. Found invalid $invalid !⊆ $(comp.psym)."))
        end
    end
    if !(all(ucond_affect, cb.affect.sym))
        invalid = filter(!ucond_affect, cb.affect.sym)
        throw(ArgumentError("All u symbols in the callback affect must be variables (in contrast to condition, observables are not allowed here). Found invalid $invalid !⊆ $(comp.sym)."))
    end
    if !(all(pcond, cb.affect.psym))
        invalid = filter(!pcond, cb.affect.psym)
        throw(ArgumentError("All p symbols in the callback affect must be parameters. Found invalid $invalid !⊆ $(comp.psym)."))
    end
    cb
end
