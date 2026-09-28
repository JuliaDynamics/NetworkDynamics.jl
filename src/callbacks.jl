"""
    abstract type ComponentCallback end

Abstract type for a component based callback. A component callback
bundles a [`ComponentCondition`](@ref) as well as a [`ComponentAffect`](@ref)
which can be then tied to a component model using [`add_callback!`](@ref) or
[`set_callback!`](@ref).

On a Network level, you can automatically create network wide `CallbackSet`s using
[`get_callbacks`](@ref).

See
[`ContinuousComponentCallback`](@ref) and [`VectorContinuousComponentCallback`](@ref) for concrete
implementations of this abstract type.
"""
abstract type ComponentCallback end

"""
    ComponentCondition(f::Function, sym)

Creates a callback condition for a [`ComponentCallback`].
- `f`: The condition function. Must be a function of the form `out=f(u, t)`
  when used for [`ContinuousComponentCallback`](@ref) or
  [`DiscreteComponentCallback`](@ref) and `f!(out, u, t)` when used for
  [`VectorContinuousComponentCallback`](@ref).
  - Arguments of `f`
    - `u`: The current values of the selected `sym` symbols, provided as a [`SymbolicView`](@ref) object.
    - `t`: The current simulation time.
- `sym`: A vector or tuple of symbols naming **states, parameters, inputs, outputs or observed**
  of the component model. Determines which values will be available through `u` in `f`.

# Example
Consider a component model with states `[:u1, :u2]`, inputs `[:i]`, outputs
`[:o]` and parameters `[:p1, :p2]`.

    ComponentCondition([:u1, :o, :p1]) do u, t
        # access values symbolically or via int index
        u[:u1] == u[1]
        u[:o] == u[2]
        u[:p1] == u[3]
        # `:u2`, `:i` and `:p2` are not available as they are not listed in `sym`.
    end

The legacy form `ComponentCondition(f, sym, psym)` with `f(u, p, t)` is still supported.
"""
struct ComponentCondition{C,DIM}
    f::C
    sym::NTuple{DIM,Symbol}
    function ComponentCondition(f, sym)
        if !hasmethod(f, Tuple{SymbolicView, Float64}) &&
           !hasmethod(f, Tuple{Vector{Float64}, SymbolicView, Float64})
            throw(ArgumentError(
                "The provided condition function has no method defined for (u, t). Available method signatures are:\n$(methods(f))\n"
            ))
        end
        new{typeof(f), length(sym)}(f, Tuple(sym))
    end
end

"""
    ComponentAffect(f::Function, sym)

Creates a callback affect for a [`ComponentCallback`].
- `f`: The affect function. Must be a function of the form `f(u, [event_signs], ctx)` where `event_signs`
  is only available in [`VectorContinuousComponentCallback`](@ref).
  - Arguments of `f`
    - `u`: The current values of the selected `sym` symbols, provided as a [`SymbolicView`](@ref) object.
      Entries which are states or parameters of the component can be written to, all other entries
      (inputs, outputs, observed) are read only.
    - `event_signs`: Only for [`VectorContinuousComponentCallback`](@ref): a length-`len` vector of
      `Int8`s encoding, for each condition output `i`, whether it crossed (`0` no crossing, `+1` upcrossing,
      `-1` downcrossing). The affect resolves the direction and any simultaneous crossings itself.
    - `ctx::NamedTuple` a named tuple with context variables.
       - `ctx.model`: a reference to the component model
       - `ctx.vidx`/`ctx.eidx`: The index of the vertex/edge model.
       - `ctx.src`/`ctx.dst`: src and dst indices (only for edge models).
       - `ctx.integrator`: The integrator object. Use [`extract_nw`](@ref) to obtain the network.
       - `ctx.t=ctx.integrator.t`: The current simulation time.
       - `ctx.dt_reset::Ref{Bool}`: Set `ctx.dt_reset[] = false` to skip the automatic
         [`SciMLBase.auto_dt_reset!`](@extref) after a change of `u`. Meant for affects which
         only store a value and don't introduce a discontinuity; parameter changes are saved either way.
         If several affects fire at the same time, one of them asking for the reset is enough.
- `sym`: A vector or tuple of symbols naming **states, parameters, inputs, outputs or observed**
  of the component model. Determines which values will be available through `u` in `f`.

The values in `u` are a snapshot taken when the affect fires. Writing a state or parameter updates
the integrator immediately, but observed entries which depend on it are not recomputed.

# Example
Consider a component model with states `[:u1, :u2]`, inputs `[:i]`, outputs
`[:o]` and parameters `[:p1, :p2]`.

    ComponentAffect([:u1, :p1, :o]) do u, ctx
        u[:u1] = 0 # change the state
        u[:p1] = u[:o] # change the parameter based on the observed
        @info "Changed :u1 and :p1 on vertex \$(ctx.vidx)" # access context
    end

The legacy form `ComponentAffect(f, sym, psym)` with `f(u, p, ctx)` is still supported.
"""
struct ComponentAffect{A,DIM}
    f::A
    sym::NTuple{DIM,Symbol}
    function ComponentAffect(f, sym)
        if !hasmethod(f, Tuple{SymbolicView, NamedTuple}) &&
           !hasmethod(f, Tuple{SymbolicView, AbstractVector{Int8}, NamedTuple})
            throw(ArgumentError(
                "The provided affect function has no method defined for (u, ctx). Available method signatures are:\n$(methods(f))\n"
            ))
        end
        new{typeof(f), length(sym)}(f, Tuple(sym))
    end
end

"""
    ContinuousComponentCallback(condition, affect; affect_neg! = affect, kwargs...)

Connect a [`ComponentCondition`](@ref) and a [`ComponentAffect`](@ref) to a
continuous callback which can be attached to a component model using
[`add_callback!`](@ref) or [`set_callback!`](@ref).

The `affect_neg!` is also a `ComponentAffect` but will be triggered on downcrossing.
It defaults to the same `affect` as on upcrossing.

The `kwargs` will be forwarded to the `VectorContinuousCallback` when the component based
callbacks are collected for the whole network using `get_callbacks`.
[`DiffEq.jl docs`](https://docs.sciml.ai/DiffEqDocs/stable/features/callback_functions/)
for available options.
"""
struct ContinuousComponentCallback{
    C  <: ComponentCondition,
    A  <: ComponentAffect,
    An <: Union{ComponentAffect,Nothing}
} <: ComponentCallback
    condition::C
    affect::A
    affect_neg::An
    kwargs::NamedTuple
end
function ContinuousComponentCallback(condition, affect; affect_neg! = affect, kwargs...)
    _assert_scalar_signature(condition) # (u, t)
    _assert_scalar_signature(affect)    # (u, ctx)
    isnothing(affect_neg!) || _assert_scalar_signature(affect_neg!) # (u, ctx)
    ContinuousComponentCallback(condition, affect, affect_neg!, NamedTuple(kwargs))
end

"""
    VectorContinuousComponentCallback(condition, affect, len; kwargs...)

Connect a [`ComponentCondition`](@ref) and a [`ComponentAffect`](@ref) to a
continuous callback which can be attached to a component model using
[`add_callback!`](@ref) or [`set_callback!`](@ref). This vector version allows
for `conditions` which have `len` output dimensions.

Mirroring `VectorContinuousCallback` from DiffEq, this callback has no separate
`affect_neg!`. Instead, the single `affect` is triggered with an additional
`event_signs` argument: a length-`len` vector of `Int8`s where entry `i` encodes
what happened to the `i`-th condition output:
- `0`: no zerocrossing detected in this dimension,
- `+1`: upcrossing (condition went from negative to positive),
- `-1`: downcrossing (condition went from positive to negative).

The `affect` fires once whenever *any* dimension crossed and is responsible for
resolving the crossing direction (and any simultaneous crossings) itself. See
[`ComponentAffect`](@ref) for the affect signature.

The `kwargs` will be forwarded to the `VectorContinuousCallback` when the component based
callbacks are collected for the whole network using [`get_callbacks(::Network)`](@ref).
[`DiffEq.jl docs`](https://docs.sciml.ai/DiffEqDocs/stable/features/callback_functions/)
for available options.
"""
struct VectorContinuousComponentCallback{
    C  <: ComponentCondition,
    A  <: ComponentAffect,
} <: ComponentCallback
    condition::C
    affect::A
    len::Int
    kwargs::NamedTuple
end
function VectorContinuousComponentCallback(condition, affect, len; kwargs...)
    _assert_vector_signature(condition) # (out, u, t)
    _assert_vector_signature(affect)    # (u, event_signs, ctx)
    VectorContinuousComponentCallback(condition, affect, len, NamedTuple(kwargs))
end

"""
    DiscreteComponentCallback(condition, affect; kwargs...)

Connect a [`ComponentCondition`](@ref) and a [`ComponentAffect`](@ref) to a
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
    _assert_scalar_signature(condition) # (u, t)
    _assert_scalar_signature(affect)    # (u, ctx)
    DiscreteComponentCallback(condition, affect, NamedTuple(kwargs))
end

"""
    PresetTimeComponentCallback(ts, affect; kwargs...)

Trigger a [`ComponentAffect`](@ref) at given timesteps `ts` in discrete
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
    _assert_scalar_signature(affect) # (u, ctx)
    PresetTimeComponentCallback(ts, affect, NamedTuple(kwargs))
end

# The condition/affect constructors accept both the scalar and the vector arity, so a function
# with the old `(u, p, t)` shape would silently pass as a vector `(out, u, t)`. The callback
# constructors know which arity they need and catch that here.
function _assert_scalar_signature(c::ComponentCondition)
    hasmethod(c.f, Tuple{SymbolicView, Float64}) && return
    throw(ArgumentError("The condition function must have the signature f(u, t) for this callback \
        type. Got a function with signatures:\n$(methods(c.f))\n\
        If it has the form f(u, p, t), use ComponentCondition(f, sym, psym) or move the parameters into the symbol list."))
end
function _assert_vector_signature(c::ComponentCondition)
    hasmethod(c.f, Tuple{Vector{Float64}, SymbolicView, Float64}) && return
    throw(ArgumentError("The condition function must have the signature f!(out, u, t) for a vector \
        callback. Got a function with signatures:\n$(methods(c.f))\n"))
end
function _assert_scalar_signature(a::ComponentAffect)
    hasmethod(a.f, Tuple{SymbolicView, NamedTuple}) && return
    throw(ArgumentError("The affect function must have the signature f(u, ctx) for this callback \
        type. Got a function with signatures:\n$(methods(a.f))\n\
        If it has the form f(u, p, ctx), use ComponentAffect(f, sym, psym) or move the parameters into the symbol list."))
end
function _assert_vector_signature(a::ComponentAffect)
    hasmethod(a.f, Tuple{SymbolicView, AbstractVector{Int8}, NamedTuple}) && return
    throw(ArgumentError("The affect function must have the signature f(u, event_signs, ctx) for a vector \
        callback. Got a function with signatures:\n$(methods(a.f))\n"))
end

# accessors
getcondition(cb::DiscreteComponentCallback) = cb.condition
getcondition(cb::ContinuousComponentCallback) = cb.condition
getcondition(cb::VectorContinuousComponentCallback) = cb.condition
getaffect(cb::ComponentCallback) = cb.affect
getaffect_neg(cb::ContinuousComponentCallback) = cb.affect_neg

"""
    get_callbacks(nw::Network, additional_callbacks=Dict())::CallbackSet

Returns a `CallbackSet` composed of all the "component-based" callbacks in the metadata of the
Network components.

You can inject additional callbacks at that stage by passing

    get_callbacks(nw, VIndex(7) => cb)
    get_callbacks(nw, Dict(VIndex(1)=>cb1, EIndex(2)=>cb2))

which won't be stored in the metadata of the component.
"""
function get_callbacks(nw::Network, additional_callbacks=Dict())
    aliased_changed(nw; warn=true)
    cbbs = wrap_component_callbacks(nw, additional_callbacks)
    if isempty(cbbs)
        return nothing
    elseif length(cbbs) == 1
        return to_callback(only(cbbs))
    else
        # we split in discrete and continuous manually, otherwise the CallbackSet
        # construction can take forever
        discrete_cb = []
        continuous_cb = []
        for batch in cbbs
            cb = to_callback(batch)
            if cb isa SciMLBase.AbstractContinuousCallback
                push!(continuous_cb, cb)
            elseif cb isa SciMLBase.AbstractDiscreteCallback
                push!(discrete_cb, cb)
            else
                error("Unknown callback type, should never be reached. Please report this issue.")
            end
        end
        CallbackSet(Tuple(continuous_cb), Tuple(discrete_cb));
    end
end

####
#### identifying callbacks which can be combined into batches
####
wrap_component_callbacks(nw, additional_cb::Pair) = wrap_component_callbacks(nw, Dict(additional_cb))
function wrap_component_callbacks(nw, additional_callbacks=Dict())
    components = SymbolicIndex{Int,Nothing}[]
    callbacks = ComponentCallback[]
    for (comp, cb) in additional_callbacks
        @assert isnothing(comp.subidx)
        _comp = idxtype(comp)(resolvecompidx(nw, comp))
        if cb isa AbstractVector || cb isa Tuple
            for _cb in cb
                push!(components, _comp)
                push!(callbacks, _cb)
            end
        else
            push!(components, _comp)
            push!(callbacks, cb)
        end
    end
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
    idx_per_type = find_identical(callbacks; equality=_batchequal)
    batches = []
    for typeidx in idx_per_type
        batchcomps = components[typeidx]
        batchcbs = callbacks[typeidx]
        if first(batchcbs) isa Union{ContinuousComponentCallback, VectorContinuousComponentCallback}
            cb = ContinuousCallbackWrapper(nw, batchcomps, batchcbs)
        elseif first(batchcbs) isa DiscreteComponentCallback
            cb = DiscreteCallbackWrapper(nw, batchcomps, batchcbs)
        elseif first(batchcbs) isa PresetTimeComponentCallback
            # PresetTimeCallbacks cannot be batched - must be single component
            @assert length(batchcbs) == 1 "PresetTimeComponentCallback cannot be batched"
            cb = PresetTimeCallbackWrapper(nw, only(batchcomps), only(batchcbs))
        else
            error("Unknown callback type, should never be reached. Please report this issue.")
        end
        push!(batches, cb)
    end
    return batches
end
_batchequal(a, b) = false
function _batchequal(a::ContinuousComponentCallback, b::ContinuousComponentCallback)
    _batchequal(a.condition, b.condition) || return false
    _batchequal(a.kwargs, b.kwargs)       || return false
    return true
end
function _batchequal(a::VectorContinuousComponentCallback, b::VectorContinuousComponentCallback)
    _batchequal(a.condition, b.condition) || return false
    _batchequal(a.kwargs, b.kwargs)       || return false
    a.len == b.len                       || return false
    return true
end
function _batchequal(a::DiscreteComponentCallback, b::DiscreteComponentCallback)
    _batchequal(a.condition, b.condition) || return false
    _batchequal(a.kwargs, b.kwargs)       || return false
    return true
end
function _batchequal(a::ComponentCondition, b::ComponentCondition)
    typeof(a) == typeof(b) || return false
    a.f === b.f
end
function _batchequal(a::NamedTuple, b::NamedTuple)
    length(a) == length(b) || return false
    for (k, v) in pairs(a)
        haskey(b, k) || return false
        v == b[k] || return false
    end
    return true
end

# a callback wrapper is a container, which wraps a network, a component callback
# and a component index. It is used for bookkeeping to know to which component
# each callback belongs to
abstract type CallbackWrapper end

# Generic functions for all CallbackWrappers with components and callbacks fields
Base.length(cw::CallbackWrapper) = length(cw.callbacks)
cbtype(cw::CallbackWrapper) = eltype(cw.callbacks)

@inline condition_dim(cw::CallbackWrapper) = first(cw.callbacks).condition.sym |> length

@inline affect_dim(cw::CallbackWrapper, aff_or_cond, i) = aff_or_cond(cw.callbacks[i]).sym |> length
@inline function affect_dim(cw::CallbackWrapper, _::typeof(getaffect_neg), i)
    affect_neg = getaffect_neg(cw.callbacks[i])
    isnothing(affect_neg) ? 0 : length(affect_neg.sym)
end

@inline condition_urange(cw::CallbackWrapper, i) = (1 + (i-1)*condition_dim(cw)) : i*condition_dim(cw)
@inline function affect_urange(cw::CallbackWrapper, aff_or_affneg, i)
    offset = sum(j -> affect_dim(cw, aff_or_affneg, j), 1:(i-1), init=0) # collect dimension before
    offset + 1 : offset + affect_dim(cw, aff_or_affneg, i)
end

# flat list of symbolic indices for all members of the batch, in slot order
function collect_c_or_a_indices(cw::CallbackWrapper, accessor)
    sidxs = SymbolicIndex[]
    for (component, cb) in zip(cw.components, cw.callbacks)
        if accessor == getaffect_neg && accessor(cb) === nothing
            continue
        end
        append!(sidxs, _symidxs(component, accessor(cb).sym))
    end
    sidxs
end
_symidxs(component::SymbolicIndex, syms) = collect(idxtype(component)(component.compidx, syms))

####
#### observed function shared by the callback wrappers
####
# The observed function of a callback batch. It sits in an untyped field, so the callback type
# doesn't carry the network type. The call through it is dynamic, and a dynamic call has to
# allocate a box for a plain Float64 argument. So the time is passed in a `Ref` that is allocated
# once, and the closure reads it back out.
struct CallbackObsf
    f::Any
    tref::Base.RefValue{Float64}
end
function CallbackObsf(obsf)
    f = (u, p, tref, out) -> obsf(u, p, tref[], out)
    CallbackObsf(f, Ref(0.0))
end
function (o::CallbackObsf)(u, p, t, out)
    if t isa Float64
        o.tref[] = t
        o.f(u, p, o.tref, out)
    else
        o.f(u, p, Ref(t), out)
    end
end

####
#### gather and write-through for affects
####
# An affect gets its `u` from a scratch buffer which is filled by the observed function of all
# requested symbols. Each slot knows where a write should land in `integrator.u` or
# `integrator.p`, so states and parameters are writable while everything else is read only.
struct AffectAccess{DC}
    obsf::CallbackObsf
    cache::DC
    uidx::Vector{Int} # position in u, 0 if the slot is no state
    pidx::Vector{Int} # position in p, 0 if the slot is no parameter
    changed::Vector{Bool} # [state written, parameter written], reset per fire
    dt_reset::Base.RefValue{Bool} # exposed as ctx.dt_reset, affects may opt out of the step reset
end
function AffectAccess(nw, symidxs)
    missing = filter(s -> !(SII.is_variable(nw, s) || SII.is_parameter(nw, s) || SII.is_observed(nw, s)), symidxs)
    if !isempty(missing)
        throw(ArgumentError("Cannot build callback as it contains references to undefined symbols: $(missing)"))
    end
    obsf = CallbackObsf(SII.observed(nw, symidxs))
    cache = DiffCache(zeros(length(symidxs)), ad_chunksize(nw.im))
    uidx = Int[something(SII.variable_index(nw, s), 0) for s in symidxs]
    pidx = Int[something(SII.parameter_index(nw, s), 0) for s in symidxs]
    AffectAccess(obsf, cache, uidx, pidx, [false, false], Ref(true))
end
function _gather(acc::AffectAccess, integrator)
    scratch = PreallocationTools.get_tmp(acc.cache, integrator.u)
    acc.obsf(integrator.u, integrator.p, integrator.t, scratch)
    scratch
end
function _affect_view(acc::AffectAccess, scratch, integrator, range, syms)
    fill!(acc.changed, false)
    acc.dt_reset[] = true
    wt = WriteThrough(view(scratch, range), integrator.u, integrator.p,
                      view(acc.uidx, range), view(acc.pidx, range), syms, acc.changed)
    SymbolicView(wt, syms)
end
# Several members of a batch may fire at the same event time. Each affect reports whether it
# wants a step reset and whether it changed a parameter, the batch collects the flags and
# acts once per event: one member asking for the reset is enough.
function _affect_flags(u::SymbolicView{<:Any,<:WriteThrough}, ctx)
    wt = u.v
    dt_reset = ctx.dt_reset[] && (uchanged(wt) || pchanged(wt))
    (dt_reset, pchanged(wt))
end
function _finish_affects!(integrator, dt_reset::Bool, pchanged::Bool)
    dt_reset && SciMLBase.auto_dt_reset!(integrator)
    pchanged && save_parameters!(integrator)
    nothing
end

####
#### wrapping of continuous callbacks
####
struct ContinuousCallbackWrapper{T<:ComponentCallback,C,ST<:SymbolicIndex} <: CallbackWrapper
    nw::Network
    components::Vector{ST}
    callbacks::Vector{T}
    sublen::Int # length of each callback
    condition::C
end
function ContinuousCallbackWrapper(nw, components, callbacks)
    if !isconcretetype(eltype(components))
        components = [c for c in components]
    end
    if !isconcretetype(eltype(callbacks))
        callbacks = [cb for cb in callbacks]
    end
    sublen = eltype(callbacks) <: ContinuousComponentCallback ? 1 : first(callbacks).len
    condition = first(callbacks).condition.f
    ContinuousCallbackWrapper(nw, components, callbacks, sublen, condition)
end

# Continuous-specific functions (for vector callbacks)
condition_outrange(ccw::ContinuousCallbackWrapper, i) = (1 + (i-1)*ccw.sublen) : i*ccw.sublen

# generate VectorContinuousCallback from a ContinuousCallbackWrapper
#
# SciMLBase>=3 dropped `affect_neg!` from `VectorContinuousCallback`; instead the
# affect receives an `event_signs::Vector{Int8}` of length `len` (`+1` upcrossing,
# `-1` downcrossing, `0` no event for each output). Component callbacks are always
# batched into such a `VectorContinuousCallback`, but the two component-callback types
# expose different user interfaces (mirroring DiffEq's own `ContinuousCallback` vs
# `VectorContinuousCallback`):
# - `ContinuousComponentCallback` keeps a separate `affect`/`affect_neg`; each callback
#   occupies a single output slot and we dispatch per crossing direction.
# - `VectorContinuousComponentCallback` has a single `affect` which receives the slice of
#   `event_signs` belonging to the component and resolves direction/simultaneity itself.
function to_callback(ccw::ContinuousCallbackWrapper)
    kwargs = first(ccw.callbacks).kwargs
    cond = _batch_condition(ccw)

    len = ccw.sublen * length(ccw.callbacks)
    affect = if cbtype(ccw) <: ContinuousComponentCallback
        _batch_scalar_affect(ccw)
    else # VectorContinuousComponentCallback
        _batch_vector_affect(ccw)
    end
    VectorContinuousCallback(cond, affect, len; kwargs...)
end
function _batch_condition(ccw::ContinuousCallbackWrapper)
    symidxs = collect_c_or_a_indices(ccw, getcondition)
    ucache = DiffCache(zeros(length(symidxs)), ad_chunksize(ccw.nw.im))
    obsf = CallbackObsf(SII.observed(ccw.nw, symidxs))

    (out, u, t, integrator) -> begin
        us = PreallocationTools.get_tmp(ucache, u)
        obsf(u, integrator.p, t, us) # fills us inplace

        for i in 1:length(ccw)
            uv = view(us, condition_urange(ccw, i))
            _u = SymbolicView(uv, ccw.callbacks[i].condition.sym)

            if cbtype(ccw) <: ContinuousComponentCallback
                oidx = only(condition_outrange(ccw, i))
                out[oidx] = ccw.condition(_u, t)
            elseif cbtype(ccw) <: VectorContinuousComponentCallback
                @views _out = out[condition_outrange(ccw, i)]
                ccw.condition(_out, _u, t)
            else
                error()
            end
        end
        nothing
    end
end
# affect builder for `ContinuousComponentCallback` batches. Each member owns one output slot, so
# the index into `event_signs` is the member index and the sign picks its up or down affect.
# The snapshots are taken once per event, so all members see the state before any affect ran.
function _batch_scalar_affect(ccw::ContinuousCallbackWrapper)
    pos_acc, pos_affect = _scalar_member_affect(ccw, getaffect)
    neg_acc, neg_affect = _scalar_member_affect(ccw, getaffect_neg)

    (integrator, event_signs) -> begin
        pos_scratch = any(>(0), event_signs) ? _gather(pos_acc, integrator) : nothing
        neg_scratch = any(<(0), event_signs) ? _gather(neg_acc, integrator) : nothing
        any_dt_reset = false
        any_pchanged = false
        for i in eachindex(event_signs)
            s = event_signs[i]
            dt_reset, pchanged = if s > 0
                pos_affect(integrator, pos_scratch, i)
            elseif s < 0
                neg_affect(integrator, neg_scratch, i)
            else
                (false, false)
            end
            any_dt_reset |= dt_reset
            any_pchanged |= pchanged
        end
        _finish_affects!(integrator, any_dt_reset, any_pchanged)
    end
end
# returns the access object and a per-member affect; the caller gathers the scratch once per event
function _scalar_member_affect(ccw::ContinuousCallbackWrapper, aff_or_affneg::F) where {F}
    acc = AffectAccess(ccw.nw, collect_c_or_a_indices(ccw, aff_or_affneg))

    affect_fn = (integrator, scratch, i) -> begin
        affect = aff_or_affneg(ccw.callbacks[i])
        isnothing(affect) && return (false, false) # affect_neg may be absent

        _u = _affect_view(acc, scratch, integrator, affect_urange(ccw, aff_or_affneg, i), affect.sym)
        ctx = get_ctx(integrator, ccw.components[i], acc)
        affect.f(_u, ctx)
        _affect_flags(_u, ctx)
    end
    acc, affect_fn
end

# affect builder for `VectorContinuousComponentCallback` batches. Each member owns `sublen`
# output slots, several of which may cross at once. Its affect is still called only once and
# receives the whole slice of `event_signs` (`0`/`+1`/`-1`) to sort out directions itself.
function _batch_vector_affect(ccw::ContinuousCallbackWrapper)
    acc = AffectAccess(ccw.nw, collect_c_or_a_indices(ccw, getaffect))

    (integrator, event_signs) -> begin
        scratch = _gather(acc, integrator)
        any_dt_reset = false
        any_pchanged = false
        for i in 1:length(ccw)
            outrange = condition_outrange(ccw, i)
            any(oidx -> !iszero(event_signs[oidx]), outrange) || continue

            affect = getaffect(ccw.callbacks[i])
            _u = _affect_view(acc, scratch, integrator, affect_urange(ccw, getaffect, i), affect.sym)
            ctx = get_ctx(integrator, ccw.components[i], acc)
            signs = view(event_signs, outrange)
            affect.f(_u, signs, ctx)
            dt_reset, pchanged = _affect_flags(_u, ctx)
            any_dt_reset |= dt_reset
            any_pchanged |= pchanged
        end
        _finish_affects!(integrator, any_dt_reset, any_pchanged)
    end
end

####
#### wrapping of discrete callbacks
####
struct DiscreteCallbackWrapper{ST,T,C} <: CallbackWrapper
    nw::Network
    components::Vector{ST}  # Changed to support batching
    callbacks::Vector{T}    # Changed to support batching
    condition::C            # Store condition function to avoid dynamic dispatch
end
function DiscreteCallbackWrapper(nw, components, callbacks)
    @assert nw isa Network
    @assert all(c -> c isa SymbolicIndex, components)
    @assert all(cb -> cb isa DiscreteComponentCallback, callbacks)  # Only DiscreteComponentCallback
    if !isconcretetype(eltype(components))
        components = [c for c in components]
    end
    if !isconcretetype(eltype(callbacks))
        callbacks = [cb for cb in callbacks]
    end
    # Extract condition function - all callbacks in batch have identical conditions
    condition = first(callbacks).condition.f
    DiscreteCallbackWrapper{eltype(components),eltype(callbacks),typeof(condition)}(nw, components, callbacks, condition)
end

# generate a DiscreteCallback from a DiscreteCallbackWrapper
#
# A `DiscreteCallback` condition is a single Bool, so the batch fires if any member does. The
# solver calls the affect right after a true condition on the same state, so the condition
# records which members fired in `fired` and the affect just reads it back.
function to_callback(dcw::DiscreteCallbackWrapper)
    kwargs = first(dcw.callbacks).kwargs
    fired = fill(false, length(dcw))
    cond = _batch_condition(dcw, fired)
    affect = _batch_affect(dcw, fired)
    DiscreteCallback(cond, affect; kwargs...)
end
function _batch_condition(dcw::DiscreteCallbackWrapper, fired)
    symidxs = collect_c_or_a_indices(dcw, getcondition)
    ucache = DiffCache(zeros(length(symidxs)), ad_chunksize(dcw.nw.im))
    obsf = CallbackObsf(SII.observed(dcw.nw, symidxs))

    (u, t, integrator) -> begin
        us = PreallocationTools.get_tmp(ucache, u)
        obsf(u, integrator.p, t, us) # fills us inplace

        for i in 1:length(dcw)
            uv = view(us, condition_urange(dcw, i))
            _u = SymbolicView(uv, dcw.callbacks[i].condition.sym)
            fired[i] = dcw.condition(_u, t)
        end
        return any(fired)
    end
end
function _batch_affect(dcw::DiscreteCallbackWrapper, fired)
    acc = AffectAccess(dcw.nw, collect_c_or_a_indices(dcw, getaffect))

    (integrator) -> begin
        scratch = _gather(acc, integrator)
        any_dt_reset = false
        any_pchanged = false
        for i in 1:length(dcw)
            fired[i] || continue

            affect = getaffect(dcw.callbacks[i])
            _u = _affect_view(acc, scratch, integrator, affect_urange(dcw, getaffect, i), affect.sym)
            ctx = get_ctx(integrator, dcw.components[i], acc)
            affect.f(_u, ctx)
            dt_reset, pchanged = _affect_flags(_u, ctx)
            any_dt_reset |= dt_reset
            any_pchanged |= pchanged
        end
        _finish_affects!(integrator, any_dt_reset, any_pchanged)
    end
end

####
#### wrapping of preset time callbacks
####
struct PresetTimeCallbackWrapper{ST,T}
    nw::Network
    component::ST   # Single component - PresetTime callbacks cannot be batched
    callback::T     # Single callback - PresetTime callbacks cannot be batched
    function PresetTimeCallbackWrapper(nw, component::SymbolicIndex, callback::PresetTimeComponentCallback)
        @assert nw isa Network
        @assert component isa SymbolicIndex
        @assert callback isa PresetTimeComponentCallback
        # PresetTimeCallbacks cannot be batched, so always single component/callback
        new{typeof(component), typeof(callback)}(nw, component, callback)
    end
end

# generate a PresetTimeCallback from a PresetTimeCallbackWrapper
function to_callback(ptcw::PresetTimeCallbackWrapper)
    callback = ptcw.callback
    component = ptcw.component
    kwargs = callback.kwargs
    ts = callback.ts

    affect = getaffect(callback)
    symidxs = _symidxs(component, affect.sym)
    acc = AffectAccess(ptcw.nw, symidxs)

    affect_fn = (integrator) -> begin
        scratch = _gather(acc, integrator)
        _u = _affect_view(acc, scratch, integrator, 1:length(symidxs), affect.sym)
        ctx = get_ctx(integrator, component, acc)
        affect.f(_u, ctx)
        _finish_affects!(integrator, _affect_flags(_u, ctx)...)
    end

    DiffEqCallbacks.PresetTimeCallback(ts, affect_fn; kwargs...)
end


####
#### generate the context for the callback effects
####
function get_ctx(integrator, sym::VIndex, acc::AffectAccess)
    nw = extract_nw(integrator)
    idx = sym.compidx
    (; integrator, t=integrator.t, model=nw[sym], vidx=idx, dt_reset=acc.dt_reset)
end
function get_ctx(integrator, sym::EIndex, acc::AffectAccess)
    nw = extract_nw(integrator)
    idx = sym.compidx
    edge = nw.im.edgevec[idx]
    (; integrator, t=integrator.t, model=nw[sym], eidx=idx, src=edge.src, dst=edge.dst, dt_reset=acc.dt_reset)
end

####
#### Internal function to check cb compat when added as metadata
####
assert_cb_compat(comp::ComponentModel, t::Tuple) = assert_cb_compat.(Ref(comp), t)
function assert_cb_compat(comp::ComponentModel, cb)
    insym = hasinsym(comp) ? insym_all(comp) : []
    named = Set(sym(comp)) ∪ Set(psym(comp)) ∪ Set(comp.obssym) ∪ insym ∪ outsym_flat(comp)

    hints = String[]
    check = (what, syms) -> begin
        invalid = filter(∉(named), syms)
        isempty(invalid) && return
        push!(hints, "All symbols in the callback $what must be states, parameters, inputs, outputs or observed of the component. Found invalid $invalid !⊆ $named.")
    end
    cb isa PresetTimeComponentCallback || check("condition", cb.condition.sym)
    check("affect", getaffect(cb).sym)
    if cb isa ContinuousComponentCallback && !isnothing(getaffect_neg(cb)) && getaffect_neg(cb) != getaffect(cb)
        check("affect_neg!", getaffect_neg(cb).sym)
    end
    if !isempty(hints)
        pushfirst!(hints, "The callback is not compatible with the component model $(comp). Issues found:")
        throw(ArgumentError(join(hints, "\n - ")))
    end

    cb
end

