function (nw::Network)(du::dT, u::T, p, t; perturb=nothing, perturb_maps=nothing, RET=Val(:du)) where {dT,T}
    if dT isa AbstractVector && !(eachindex(du) == eachindex(u) == 1:nw.im.lastidx_dynamic)
        throw(ArgumentError("du or u does not have expected size $(nw.im.lastidx_dynamic)"))
    end
    if pdim(nw) > 0 && !(eachindex(p) == 1:nw.im.lastidx_p)
        throw(ArgumentError("p does not has expecte size $(nw.im.lastidx_p)"))
    end
    cacheT = _cachetype(du, u, p, t, perturb, RET)

    # For an untyped core this call is the function barrier. A dynamic call boxes every argument
    # which is not a heap object, so the time goes through `tbuf` if it fits.
    if _stash_time!(nw.tbuf, t)
        nw.core(Val(cacheT), du, u, p, nw.tbuf, Val(typeof(t)), perturb, perturb_maps, RET)
    else
        nw.core(Val(cacheT), du, u, p, t, perturb, perturb_maps, RET)
    end
    return nothing
end
function _cachetype(du, u, p, t, perturb, RET)
    if RET == Val(:du)
        eltype(du) # if du is required du needs to be preallocated with the correct type!
    else
        prelimCacheT = cachetype(du, u, p, perturb)
        # don't use t to widen cache type if its just a plain number or nothing;
        # plain Float64 t would otherwise upcast Float32 GPU arrays, and nothing t (from SII) corrupts cacheT
        if t isa AbstractFloat || isnothing(t)
            prelimCacheT
        else
            promote_type(prelimCacheT, typeof(t))
        end
    end
end
_stash_time!(tbuf::Vector{FT}, t::FT) where {FT} = (tbuf[1] = t; true)
function _stash_time!(tbuf::Vector{FT}, t::ForwardDiff.Dual{<:Any,FT,N}) where {FT,N}
    N + 1 ≤ length(tbuf) || return false
    reinterpret(typeof(t), view(tbuf, 1:N+1))[1] = t
    true
end
_stash_time!(_, _) = false

# entry from an untyped network: read the time back and run the coreloop
function (core::NetworkCore)(cacheT::Val, du, u, p, tbuf, ::Val{TT}, perturb, perturb_maps, RET) where {TT}
    t = _unstash_time(tbuf, TT)
    core(cacheT, du, u, p, t, perturb, perturb_maps, RET)
end
_unstash_time(tbuf::Vector{FT}, ::Type{FT}) where {FT} = tbuf[1]
function _unstash_time(tbuf::Vector{FT}, ::Type{TT}) where {FT,N,TT<:ForwardDiff.Dual{<:Any,FT,N}}
    reinterpret(TT, view(tbuf, 1:N+1))[1]
end

# the actual coreloop
function (nwc::NetworkCore)(::Val{cacheT}, du, u, p, t, perturb, perturb_maps, RET) where {cacheT}
    ex = executionstyle(nwc)
    isnothing(du) || fill!(du, zero(eltype(du)))
    o = get_output_cache(nwc, cacheT)
    extbuf = has_external_input(nwc) ? get_extinput_cache(nwc, cacheT) : nothing

    duopt = (du, u, o, p, t)
    aggbuf = get_aggregation_cache(nwc, cacheT)
    fill!(aggbuf, _appropriate_zero(aggbuf))
    gbuf = get_gbuf(nwc.gbufprovider, o)

    # EARLY Return if just buffers are needed, they are read with `get_buffers`
    RET isa Val{:buf_uninit} && return nothing

    # vg without ff
    process_batches!(ex, Val{:g}(), !hasff, nwc.vertexbatches, (nothing, nothing), duopt)
    # eg without ff
    process_batches!(ex, Val{:g}(), !hasff, nwc.layer.edgebatches, (nothing, nothing), duopt)

    # process batches might be async so sync before next step
    ex isa KAExecution && KernelAbstractions.synchronize(get_backend(du))

    # if loopback edges are present, we direclty copy cluster output to satelite input
    !isnothing(nwc.loopbackmap) && apply_loopback!(aggbuf, o, nwc.loopbackmap)

    if !isnothing(perturb)
        @assert aggfun(nwc.layer.aggregator) == (+) "currently only + aggregation is supported for vertex input perturbation, got $(nwc.layer.aggregator)"
        apply_perturb!(aggbuf, perturb, perturb_maps.vi_map)
    end

    # process vg WITH ff (only allowed on loopback edges)
    process_batches!(ex, Val{:g}(), hasff, nwc.vertexbatches, (aggbuf, nothing), duopt)

    # process batches might be async so sync before next step
    ex isa KAExecution && KernelAbstractions.synchronize(get_backend(du))

    # gather the external inputs
    has_external_input(nwc) && collect_externals!(nwc.extmap, extbuf, u, o)

    if !isnothing(perturb)
        apply_perturb!(o, perturb, perturb_maps.vo_map)
    end
    # gather the vertex results for edges with ff
    gather!(nwc.gbufprovider, gbuf, o)

    if !isnothing(perturb)
        @assert nwc.gbufprovider isa EagerGBufProvider "edge input perturbation is only supported with buffered execution schemes!"
        apply_perturb!(gbuf, perturb, perturb_maps.ei_map)
    end

    if RET isa Val{:du} # normal coreloop
        # execute f for the edges without ff
        process_batches!(ex, Val{:f}(), !hasff, nwc.layer.edgebatches, (gbuf, extbuf), duopt)
        # execute f&g for edges with ff
        process_batches!(ex, Val{:fg}(), hasff, nwc.layer.edgebatches, (gbuf, extbuf), duopt)
    else # This is the pass if we only fill the buffers, :f does not need to be executed
        process_batches!(ex, Val{:g}(), hasff, nwc.layer.edgebatches, (gbuf, extbuf), duopt)
    end

    # process batches might be async so sync before next step
    ex isa KAExecution && KernelAbstractions.synchronize(get_backend(du))

    if !isnothing(perturb)
        apply_perturb!(o, perturb, perturb_maps.eo_map)
    end
    # aggegrate the results
    aggregate!(nwc.layer.aggregator, aggbuf, o)

    RET isa Val{:buf_init} && return nothing

    # vf for all vertices (including FF/injector type)
    process_batches!(ex, Val{:f}(), nofilt, nwc.vertexbatches, (aggbuf, extbuf), duopt)

    # process batches might be async so sync before next step
    ex isa KAExecution && KernelAbstractions.synchronize(get_backend(du))
    return nothing
end
# Output, aggregation and external input buffers, filled for the given state if `initbufs`. The
# network call returns nothing (a returned tuple would be boxed behind the barrier), so the buffers
# are read from the caches afterwards. Only type stable if the caches are passed in.
get_buffers(nw::Network, u, p, t; kwargs...) = get_buffers(nw, buffer_caches(nw), u, p, t; kwargs...)
function get_buffers(nw::Network, bc, u, p, t; initbufs=true, perturb=nothing, kwargs...)
    if initbufs
        nw(nothing, u, p, t; RET=Val(:buf_init), perturb, kwargs...)
    else
        nw(nothing, u, p, t; RET=Val(:buf_uninit), perturb, kwargs...)
    end
    T = _cachetype(nothing, u, p, t, perturb, Val(:buf_init))
    extbuf = isnothing(bc.extmap) ? nothing : get_tmp(bc.caches.external, T)
    return get_tmp(bc.caches.output, T), get_tmp(bc.caches.aggregation, T), extbuf
end

@inline function process_batches!(::SequentialExecution, fg, filt::F, batches, inbufs, duopt) where {F}
    (du, u, o, p, t) = duopt
    unrolled_foreach(filt, batches, fg, inbufs, du, u, o, p, t) do batch, fg, inbufs, du, u, o, p, t
        for i in 1:length(batch)
            _type = dispatchT(batch)
            apply_comp!(_type, fg, batch, i, du, u, o, inbufs, p, t)
        end
    end
end

@inline function process_batches!(::ThreadedExecution, fg, filt::F, batches, inbufs, duopt) where {F}
    (du, u, o, p, t) = duopt
    unrolled_foreach(filt, batches, fg, inbufs, du, u, o, p, t) do batch, fg, inbufs, du, u, o, p, t
        Threads.@threads for i in 1:length(batch)
            _type = dispatchT(batch)
            apply_comp!(_type, fg, batch, i, du, u, o, inbufs, p, t)
        end
    end
end

@inline function process_batches!(::PolyesterExecution, fg, filt::F, batches, inbufs, duopt) where {F}
    (du, u, o, p, t) = duopt
    unrolled_foreach(filt, batches, fg, inbufs, du, u, o, p, t) do batch, fg, inbufs, du, u, o, p, t
        Polyester.@batch for i in 1:length(batch)
            _type = dispatchT(batch)
            apply_comp!(_type, fg, batch, i, du, u, o, inbufs, p, t)
        end
    end
end

@inline function process_batches!(::KAExecution, fg, filt::F, batches, inbufs, duopt) where {F}
    _backend = get_backend(duopt[1])
    unrolled_foreach(filt, batches) do batch
        (du, u, o, p, t) = duopt
        _type = dispatchT(batch)
        kernel = if evalf(fg, batch) && evalg(fg, batch)
            compkernel_fg!(_backend)
        elseif evalf(fg, batch)
            compkernel_f!(_backend)
        elseif evalg(fg, batch)
            compkernel_g!(_backend)
        end
        isnothing(kernel) || kernel(_type, fg, batch, du, u, o, inbufs, p, t; ndrange=length(batch))
    end
end
@kernel function compkernel_f!(::Type{T}, @Const(fg), @Const(batch),
                               du, @Const(u), @Const(o), @Const(inbufs), @Const(p), @Const(t)) where {T}
    I = @index(Global)
    apply_comp!(T, fg, batch, I, du, u, o, inbufs, p, t)
    nothing
end
@kernel function compkernel_g!(::Type{T}, @Const(fg), @Const(batch),
                               @Const(du), @Const(u), o, @Const(inbufs), @Const(p), @Const(t)) where {T}
    I = @index(Global)
    apply_comp!(T, fg, batch, I, du, u, o, inbufs, p, t)
    nothing
end
@kernel function compkernel_fg!(::Type{T}, @Const(fg), @Const(batch),
                                du, @Const(u), o, @Const(inbuf), @Const(p), @Const(t)) where {T}
    I = @index(Global)
    apply_comp!(T, fg, batch, I, du, u, o, inbufs, p, t)
    nothing
end


@inline function apply_comp!(::Type{<:VertexModel}, fg, batch, i, du, u, o, inbufs, p, t)
    @inbounds begin
        aggbuf, extbuf = inbufs
        _o   = _needs_out(fg, batch) ? view(o, out_range(batch, i))         : nothing
        _du  = _needs_du(fg, batch)  ? view(du, state_range(batch, i))      : nothing
        _u   = _needs_u(fg, batch)   ? view(u,  state_range(batch, i))      : nothing
        _ins = _needs_in(fg, batch)  ? (view(aggbuf, in_range(batch, i)),)  : nothing
        _p   = _needs_p(fg, batch)   ? view(p,  parameter_range(batch, i))  : nothing
        if has_external_input(batch) && _needs_in(fg, batch)
            _ext = view(extbuf, extbuf_range(batch, i))
            _ins = (_ins..., _ext)
        end
        evalf(fg, batch) && apply_compf(compf(batch), _du, _u, _ins, _p, t)
        evalg(fg, batch) && apply_compg(fftype(batch), compg(batch), (_o,), _u, _ins, _p, t)
    end
    nothing
end

@inline function apply_comp!(::Type{<:EdgeModel}, fg, batch, i, du, u, o, inbufs, p, t)
    @inbounds begin
        gbuf, extbuf = inbufs
        _osrc = _needs_out(fg, batch) ? view(o, out_range(batch, i, :src))   : nothing
        _odst = _needs_out(fg, batch) ? view(o, out_range(batch, i, :dst))   : nothing
        _du   = _needs_du(fg, batch)  ? view(du, state_range(batch, i))      : nothing
        _u    = _needs_u(fg, batch)   ? view(u,  state_range(batch, i))      : nothing
        _ins  = _needs_in(fg, batch)  ? get_src_dst(gbuf, batch, i)          : nothing
        _p    = _needs_p(fg, batch)   ? view(p,  parameter_range(batch, i))  : nothing
        if has_external_input(batch) && _needs_in(fg, batch)
            _ext = view(extbuf, extbuf_range(batch, i))
            _ins = (_ins..., _ext)
        end
        evalf(fg, batch) && apply_compf(compf(batch), _du, _u, _ins, _p, t)
        evalg(fg, batch) && apply_compg(fftype(batch), compg(batch), (_osrc, _odst), _u, _ins, _p, t)
    end
    nothing
end

@propagate_inbounds function apply_compf(f::F, du, u, ins, p, t) where {F}
    f(du, u, ins..., p, t)
    nothing
end

@propagate_inbounds function apply_compg(::PureFeedForward, g::G, outs, u, ins, p, t) where {G}
    g(outs..., ins..., p, t)
    nothing
end
@propagate_inbounds function apply_compg(::FeedForward, g::G, outs, u, ins, p, t) where {G}
    g(outs..., u, ins..., p, t)
    nothing
end
@propagate_inbounds function apply_compg(::NoFeedForward, g::G, outs, u, ins, p, t) where {G}
    g(outs..., u, p, t)
    nothing
end
@propagate_inbounds function apply_compg(::PureStateMap, g::G, outs, u, ins, p, t) where {G}
    g(outs..., u)
    nothing
end

# check if the function arguments are actually used
_needs_du(fg, batch)  = evalf(fg, batch)
_needs_u(fg, batch)   = evalf(fg, batch) || fftype(batch) != PureFeedForward()
_needs_out(fg, batch) = evalg(fg, batch)
_needs_in(fg, batch)  = evalf(fg, batch) || hasff(batch)
_needs_p(fg, batch)   = !iszero(pdim(batch)) && (evalf(fg, batch) || fftype(batch) != PureStateMap())

# check if eval of f or g is necessary
evalf(::Val{:f}, batch) = !isnothing(compf(batch))
evalf(::Val{:g}, batch) = false
evalf(::Val{:fg}, batch) = !isnothing(compf(batch))
evalg(::Val{:f}, _) = false
evalg(::Val{:g}, _) = true
evalg(::Val{:fg}, _) = true

function _appropriate_zero(x)
    if isconcretetype(eltype(x))
        zero(eltype(x))
    else
        0.0 # hopefully that casts to what is needed
    end
end

function apply_perturb!(buf, perturb, map)
    for (bufidx, perturbidx) in map
        buf[bufidx] += perturb[perturbidx]
    end
end
