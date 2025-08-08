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
    ComponentCondition(f::Function, sym, psym)

Creates a callback condition for a [`ComponentCallback`].
- `f`: The condition function. Must be a function of the form `out=f(u, p, t)`
  when used for [`ContinuousComponentCallback`](@ref) or
  [`DiscreteComponentCallback`](@ref) and `f!(out, u, p, t)` when used for
  [`VectorContinuousComponentCallback`](@ref).
  - Arguments of `f`
    - `u`: The current value of the selected `sym` states, provided as a [`SymbolicView`](@ref) object.
    - `p`: The current value of the selected `psym` parameters.
    - `t`: The current simulation time.
- `sym`: A vector or tuple of symbols, which represent **states** (including
  inputs, outputs, observed) of the component model. Determines, which states will
  be available through parameter `u` in the callback condition function `f`.
- `psym`: A vector or tuple of symbols, which represent **parameters** of the component mode.
  Determines, which parameters will be available in the condition function `f`

# Example
Consider a component model with states `[:u1, :u2]`, inputs `[:i]`, outputs
`[:o]` and parameters `[:p1, :p2]`.

    ComponentCondition([:u1, :o], [:p1]) do u, p, t
        # access states symbolically or via int index
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
  is only available in [`VectorContinuousComponentCallback`](@ref).
  - Arguments of `f`
    - `u`: The current (mutable) value of the selected `sym` states, provided as a [`SymbolicView`](@ref) object.
    - `p`: The current (mutable) value of the selected `psym` parameters.
    - `event_idx`: The current event index, i.e. which `out` element triggered in case of [`VectorContinuousComponentCallback`](@ref).
    - `ctx::NamedTuple` a named tuple with context variables.
       - `ctx.model`: a reference to the component model
       - `ctx.vidx`/`ctx.eidx`: The index of the vertex/edge model.
       - `ctx.src`/`ctx.dst`: src and dst indices (only for edge models).
       - `ctx.integrator`: The integrator object. Use [`extract_nw`](@ref) to obtain the network.
       - `ctx.t=ctx.integrator.t`: The current simulation time.
- `sym`: A vector or tuple of symbols, which represent **states** (**excluding**
  inputs, outputs, observed) of the component model. Determines, which states will
  be available through parameter `u` in the callback condition function `f`.
- `psym`: A vector or tuple of symbols, which represent **parameters** of the component mode.
  Determines, which parameters will be available in the condition function `f`

# Example
Consider a component model with states `[:u1, :u2]`, inputs `[:i]`, outputs
`[:o]` and parameters `[:p1, :p2]`.

    ComponentAffect([:u1, :o], [:p1]) do u, p, ctx
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
    ContinuousComponentCallback(condition, affect; kwargs...)

Connect a [`ComponentCondition`](@ref) and a [`ComponentAffect`](@ref) to a
continuous callback which can be attached to a component model using
[`add_callback!`](@ref) or [`set_callback!`](@ref).

The `kwargs` will be forwarded to the `VectorContinuousCallback` when the component based
callbacks are collected for the whole network using `get_callbacks`.
[`DiffEq.jl docs`](https://docs.sciml.ai/DiffEqDocs/stable/features/callback_functions/)
for available options.
"""
struct ContinuousComponentCallback{C,A,CDIM,CPDIM,ADIM,APDIM} <: ComponentCallback
    condition::ComponentCondition{C,CDIM,CPDIM}
    affect::ComponentAffect{A,ADIM,APDIM}
    kwargs::NamedTuple
end
function ContinuousComponentCallback(condition, affect; kwargs...)
    if haskey(kwargs, :affect_neg!)
        throw(ArgumentError("affect_neg! not supported yet. Please raise issue."))
    end
    ContinuousComponentCallback(condition, affect, NamedTuple(kwargs))
end

"""
    VectorContinuousComponentCallback(condition, affect, len; kwargs...)

Connect a [`ComponentCondition`](@ref) and a [`ComponentAffect`](@ref) to a
continuous callback which can be attached to a component model using
[`add_callback!`](@ref) or [`set_callback!`](@ref). This vector version allows
for `conditions` which have `len` output dimensions.
The `affect` will be triggered with the additional `event_idx` argument to know in which
dimension the zerocrossing was detected.

The `kwargs` will be forwarded to the `VectorContinuousCallback` when the component based
callbacks are collected for the whole network using [`get_callbacks(::Network)`](@ref).
[`DiffEq.jl docs`](https://docs.sciml.ai/DiffEqDocs/stable/features/callback_functions/)
for available options.
"""
struct VectorContinuousComponentCallback{C,A,CDIM,CPDIM,ADIM,APDIM} <: ComponentCallback
    condition::ComponentCondition{C,CDIM,CPDIM}
    affect::ComponentAffect{A,ADIM,APDIM}
    len::Int
    kwargs::NamedTuple
end
function VectorContinuousComponentCallback(condition, affect, len; kwargs...)
    if haskey(kwargs, :affect_neg!)
        throw(ArgumentError("affect_neg! not supported yet. Please raise issue."))
    end
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

# Generic functions for all CallbackWrappers with components and callbacks fields
Base.length(cw::CallbackWrapper) = length(cw.callbacks)
cbtype(cw::CallbackWrapper) = eltype(cw.callbacks)

condition_dim(cw::CallbackWrapper)  = first(cw.callbacks).condition.sym  |> length
condition_pdim(cw::CallbackWrapper) = first(cw.callbacks).condition.psym |> length
affect_dim(cw::CallbackWrapper, i)  = cw.callbacks[i].affect.sym  |> length
affect_pdim(cw::CallbackWrapper, i) = cw.callbacks[i].affect.psym |> length

condition_urange(cw::CallbackWrapper, i) = (1 + (i-1)*condition_dim(cw))  : i*condition_dim(cw)
condition_prange(cw::CallbackWrapper, i) = (1 + (i-1)*condition_pdim(cw)) : i*condition_pdim(cw)
affect_urange(cw::CallbackWrapper, i) = (1 + (i-1)*affect_dim(cw,i) ) : i*affect_dim(cw,i)
affect_prange(cw::CallbackWrapper, i) = (1 + (i-1)*affect_pdim(cw,i)) : i*affect_pdim(cw,i)

function collect_c_or_a_indices(cw::CallbackWrapper, c_or_a, u_or_p)
    sidxs = SymbolicIndex[]
    for (component, cb) in zip(cw.components, cw.callbacks)
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

cbidx_from_outidx(ccw::ContinuousCallbackWrapper, outidx) = div(outidx-1, ccw.sublen) + 1
subidx_from_outidx(ccw::ContinuousCallbackWrapper, outidx) = mod(outidx, 1:ccw.sublen)

# generate VectorContinuousCallback from a ContinuousCallbackWrapper
function to_callback(ccw::ContinuousCallbackWrapper)
    kwargs = first(ccw.callbacks).kwargs
    cond = _batch_condition(ccw)
    affect = _batch_affect(ccw)
    len = ccw.sublen * length(ccw.callbacks)
    VectorContinuousCallback(cond, affect, len; kwargs...)
end
function _batch_condition(ccw::ContinuousCallbackWrapper)
    usymidxs = collect_c_or_a_indices(ccw, :condition, :sym)
    psymidxs = collect_c_or_a_indices(ccw, :condition, :psym)
    ucache = DiffCache(zeros(length(usymidxs)), 12)

    obsf = SII.observed(ccw.nw, usymidxs)
    pidxs = SII.parameter_index.(Ref(ccw.nw), psymidxs)

    if any(isnothing, pidxs)
        nidxs = findall(isnothing, pidxs)
        missing_p = psymidxs[nidxs]
        throw(ArgumentError("Cannot build callback as it contains refrences to undefined parameters $(missing_p)"))
    end

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

            if cbtype(ccw) <: ContinuousComponentCallback
                oidx = only(condition_outrange(ccw, i))
                out[oidx] = ccw.condition(_u, _p, t)
            elseif cbtype(ccw) <: VectorContinuousComponentCallback
                @views _out = out[condition_outrange(ccw, i)]
                ccw.condition(_out, _u, _p, t)
            else
                error()
            end
        end
        nothing
    end
end
function _batch_affect(ccw::ContinuousCallbackWrapper)
    usymidxs = collect_c_or_a_indices(ccw, :affect, :sym)
    psymidxs = collect_c_or_a_indices(ccw, :affect, :psym)

    uidxs = SII.variable_index.(Ref(ccw.nw), usymidxs)
    pidxs = SII.parameter_index.(Ref(ccw.nw), psymidxs)

    if any(isnothing, uidxs) || any(isnothing, pidxs)
        missing_u = []
        if any(isnothing, uidxs)
            nidxs = findall(isnothing, uidxs)
            append!(missing_u, usymidxs[nidxs])
        end
        missing_p = []
        if any(isnothing, pidxs)
            nidxs = findall(isnothing, pidxs)
            append!(missing_p, psymidxs[nidxs])
        end
        throw(ArgumentError(
            "Cannot build callback as it contains refrences to undefined symbols:\n"*
            (isempty(missing_u) ? "" : "Missing state symbols: $(missing_u)\n")*
            (isempty(missing_p) ? "" : "Missing param symbols: $(missing_p)\n")
        ))
    end

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
        if cbtype(ccw) <: ContinuousComponentCallback
            ccw.callbacks[i].affect.f(_u, _p, ctx)
        elseif cbtype(ccw) <: VectorContinuousComponentCallback
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
function to_callback(dcw::DiscreteCallbackWrapper)
    kwargs = first(dcw.callbacks).kwargs
    cond = _batch_condition(dcw)
    affect = _batch_affect(dcw)
    DiscreteCallback(cond, affect; kwargs...)
end
function _batch_condition(dcw::DiscreteCallbackWrapper)
    usymidxs = collect_c_or_a_indices(dcw, :condition, :sym)
    psymidxs = collect_c_or_a_indices(dcw, :condition, :psym)
    ucache = DiffCache(zeros(length(usymidxs)), 12)

    obsf = SII.observed(dcw.nw, usymidxs)
    pidxs = SII.parameter_index.(Ref(dcw.nw), psymidxs)

    if any(isnothing, pidxs)
        nidxs = findall(isnothing, pidxs)
        missing_p = psymidxs[nidxs]
        throw(ArgumentError("Cannot build callback as it contains refrences to undefined parameters $(missing_p)"))
    end

    (u, t, integrator) -> begin
        us = PreallocationTools.get_tmp(ucache, u)
        obsf(u, integrator.p, t, us) # fills us inplace

        # OR logic: return true if ANY component condition is true
        for i in 1:length(dcw)
            # symbolic view into u
            uv = view(us, condition_urange(dcw, i))
            _u = SymbolicView(uv, dcw.callbacks[i].condition.sym)

            # symbolic view into p
            pidxsv = view(pidxs, condition_prange(dcw, i))
            pv = view(integrator.p, pidxsv)
            _p = SymbolicView(pv, dcw.callbacks[i].condition.psym)

            # If any condition is true, trigger the callback
            if dcw.condition(_u, _p, t)
                return true
            end
        end
        return false
    end
end
function _batch_affect(dcw::DiscreteCallbackWrapper)
    # Setup for condition re-evaluation
    cusymidxs = collect_c_or_a_indices(dcw, :condition, :sym)
    cpsymidxs = collect_c_or_a_indices(dcw, :condition, :psym)
    cucache = DiffCache(zeros(length(cusymidxs)), 12)
    cobsf = SII.observed(dcw.nw, cusymidxs)
    cpidxs = SII.parameter_index.(Ref(dcw.nw), cpsymidxs)

    # Setup for affect execution
    ausymidxs = collect_c_or_a_indices(dcw, :affect, :sym)
    apsymidxs = collect_c_or_a_indices(dcw, :affect, :psym)

    auidxs = SII.variable_index.(Ref(dcw.nw), ausymidxs)
    apidxs = SII.parameter_index.(Ref(dcw.nw), apsymidxs)

    if any(isnothing, auidxs) || any(isnothing, apidxs)
        missing_u = []
        if any(isnothing, auidxs)
            nidxs = findall(isnothing, auidxs)
            append!(missing_u, ausymidxs[nidxs])
        end
        missing_p = []
        if any(isnothing, apidxs)
            nidxs = findall(isnothing, apidxs)
            append!(missing_p, apsymidxs[nidxs])
        end
        throw(ArgumentError(
            "Cannot build callback as it contains refrences to undefined symbols:\n"*
            (isempty(missing_u) ? "" : "Missing state symbols: $(missing_u)\n")*
            (isempty(missing_p) ? "" : "Missing param symbols: $(missing_p)\n")
        ))
    end

    (integrator) -> begin
        # Re-evaluate all conditions to determine which affects to execute
        # the affects might mutate p, therfor we ceate a copy to evaluate all
        # conditions on the unaltered state!
        cus = PreallocationTools.get_tmp(cucache, integrator.u)
        cobsf(integrator.u, integrator.p, integrator.t, cus)
        cps = copy(integrator.p)

        any_uchanged = false
        any_pchanged = false

        for i in 1:length(dcw)
            # Re-evaluate condition for component i
            cuv = view(cus, condition_urange(dcw, i))
            c_u = SymbolicView(cuv, dcw.callbacks[i].condition.sym)
            cpidxsv = view(cpidxs, condition_prange(dcw, i))
            cpv = view(cps, cpidxsv)
            c_p = SymbolicView(cpv, dcw.callbacks[i].condition.psym)

            # Only execute affect if condition is true
            if dcw.condition(c_u, c_p, integrator.t)
                # Execute affect for component i
                auidxsv = view(auidxs, affect_urange(dcw, i))
                auv = view(integrator.u, auidxsv)
                a_u = SymbolicView(auv, dcw.callbacks[i].affect.sym)

                apidxsv = view(apidxs, affect_prange(dcw, i))
                apv = view(integrator.p, apidxsv)
                a_p = SymbolicView(apv, dcw.callbacks[i].affect.psym)

                ctx = get_ctx(integrator, dcw.components[i])

                uhash = hash(auv)
                phash = hash(apv)
                dcw.callbacks[i].affect.f(a_u, a_p, ctx)
                pchanged = hash(apv) != phash
                uchanged = hash(auv) != uhash

                any_uchanged = any_uchanged || uchanged
                any_pchanged = any_pchanged || pchanged
            end
        end

        (any_uchanged || any_pchanged) && SciMLBase.auto_dt_reset!(integrator)
        any_pchanged && save_parameters!(integrator)
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

    # Create affect function for the single component
    uidxtype = component isa EIndex ? EIndex : VIndex
    pidxtype = component isa EIndex ? EPIndex : VPIndex
    usymidxs = uidxtype(component.compidx, callback.affect.sym)
    psymidxs = pidxtype(component.compidx, callback.affect.psym)

    uidxs = SII.variable_index.(Ref(ptcw.nw), usymidxs)
    pidxs = SII.parameter_index.(Ref(ptcw.nw), psymidxs)

    affect = (integrator) -> begin
        uv = view(integrator.u, uidxs)
        _u = SymbolicView(uv, callback.affect.sym)
        pv = view(integrator.p, pidxs)
        _p = SymbolicView(pv, callback.affect.psym)
        ctx = get_ctx(integrator, component)

        uhash = hash(uv)
        phash = hash(pv)
        callback.affect.f(_u, _p, ctx)
        pchanged = hash(pv) != phash
        uchanged = hash(uv) != uhash

        (pchanged || uchanged) && SciMLBase.auto_dt_reset!(integrator)
        pchanged && save_parameters!(integrator)
    end

    DiffEqCallbacks.PresetTimeCallback(ts, affect; kwargs...)
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

####
#### Backward compatibility aliases for old misspellings (deprecated)
####
const ContinousComponentCallback = ContinuousComponentCallback
const VectorContinousComponentCallback = VectorContinuousComponentCallback
