function (nw::Network{A,B,C,D,E})(du::dT, u::T, p, t) where {A,B,C,D,E,dT,T}
    if !(eachindex(du) == eachindex(u) == 1:nw.im.lastidx_dynamic)
        throw(ArgumentError("du or u does not have expected size $(nw.im.lastidx_dynamic)"))
    end
    if pdim(nw) > 0 && !(eachindex(p) == 1:nw.im.lastidx_p)
        throw(ArgumentError("p does not has expecte size $(nw.im.lastidx_p)"))
    end

    ex = executionstyle(nw)
    fill!(du, zero(eltype(du)))
    o = get_output_cache(nw, du)
    extbuf = has_external_input(nw) ? get_extinput_cache(nw, du) : nothing

    duopt = (du, u, o, p, t)
    aggbuf = get_aggregation_cache(nw, du)
    gbuf = get_gbuf(nw.gbufprovider, o)

    # vg without ff
    process_batches!(ex, Val{:g}(), !hasff, nw.vertexbatches, (nothing, nothing), duopt)
    # eg without ff
    process_batches!(ex, Val{:g}(), !hasff, nw.layer.edgebatches, (nothing, nothing), duopt)

    # process batches might be async so sync before next step
    ex isa KAExecution && KernelAbstractions.synchronize(get_backend(du))

    # gather the external inputs
    has_external_input(nw) && collect_externals!(nw.extmap, extbuf, u, o)

    # gather the vertex results for edges with ff
    gather!(nw.gbufprovider, gbuf, o)

    # execute f for the edges without ff
    process_batches!(ex, Val{:f}(), !hasff, nw.layer.edgebatches, (gbuf, extbuf), duopt)
    # execute f&g for edges with ff
    process_batches!(ex, Val{:fg}(), hasff, nw.layer.edgebatches, (gbuf, extbuf), duopt)

    # process batches might be async so sync before next step
    ex isa KAExecution && KernelAbstractions.synchronize(get_backend(du))

    # aggegrate the results
    aggregate!(nw.layer.aggregator, aggbuf, o)

    # vf for vertices without ff
    process_batches!(ex, Val{:f}(), !hasff, nw.vertexbatches, (aggbuf, extbuf), duopt)

    # process batches might be async so sync before next step
    ex isa KAExecution && KernelAbstractions.synchronize(get_backend(du))
    return nothing
end
function get_buffers(nw, u, p, t; initbufs)
    o = get_output_cache(nw, u)
    aggbuf = get_aggregation_cache(nw, u)
    extbuf = has_external_input(nw) ? get_extinput_cache(nw, u) : nothing
    if initbufs
        du = nothing
        duopt = (du, u, o, p, t)
        ex = executionstyle(nw)
        fill!(o, convert(eltype(o), NaN))
        gbuf = get_gbuf(nw.gbufprovider, o)

        # vg without ff
        process_batches!(ex, Val{:g}(), !hasff, nw.vertexbatches, (nothing, nothing), duopt)
        # eg without ff
        process_batches!(ex, Val{:g}(), !hasff, nw.layer.edgebatches, (nothing, nothing), duopt)
        # process batches might be async so sync before next step
        ex isa KAExecution && KernelAbstractions.synchronize(get_backend(u))
        # gather the external inputs
        has_external_input(nw) && collect_externals!(nw.extmap, extbuf, u, o)
        # gather the vertex results for edges with ff
        gather!(nw.gbufprovider, gbuf, o) # 2.6ms
        # execute g for edges with ff
        process_batches!(ex, Val{:g}(), hasff, nw.layer.edgebatches, (gbuf, nothing), duopt)
        # process batches might be async so sync before next step
        ex isa KAExecution && KernelAbstractions.synchronize(get_backend(u))
        # aggegrate the results
        aggregate!(nw.layer.aggregator, aggbuf, o)
    end
    o, aggbuf, extbuf
end

@inline function process_batches!(::SequentialExecution, fg, filt::F, batches, inbufs, duopt) where {F}
    unrolled_foreach(filt, batches) do batch
        (du, u, o, p, t) = duopt
        for i in 1:length(batch)
            _type = dispatchT(batch)
            apply_comp!(_type, fg, batch, i, du, u, o, inbufs, p, t)
        end
    end
end

@inline function process_batches!(ex::ThreadedExecution, fg, filt::F, batches, inbufs, duopt) where {F}
    # unrolled_foreach(filt, batches) do batch
    #     (du, u, o, p, t) = duopt
    #     Threads.@threads for i in 1:length(batch)
    #         _type = dispatchT(batch)
    #         apply_comp!(_type, fg, batch, i, du, u, o, inbufs, p, t)
    #     end
    # end
    # return

    Nchunks = Threads.nthreads()
    # Nchunks = 4
    # chunking is kinda expensive, so we cache it
    key = hash((Base.objectid(batches), filt, fg, Nchunks))
    chunks = get!(ex.chunk_cache, key) do
        _chunk_batches(batches, filt, fg, Nchunks)
    end

    _progress_in_batch = function(batch, ci, processed, N)
        (du, u, o, p, t) = duopt
        _type = dispatchT(batch)
        while ci ≤ length(batch) && processed < N
            apply_comp!(_type, fg, batch, ci, du, u, o, inbufs, p, t)
            ci += 1
            processed += 1
        end
        processed, ci
    end

    Threads.@sync for chunk in chunks
        chunk.N == 0 && continue
        Threads.@spawn begin
            local N = chunk.N
            local bi = chunk.batch_i
            local ci = chunk.comp_i
            local processed = 0
            while processed < N
                batch = batches[bi]
                filt(batch) || continue
                processed, ci = @noinline _progress_in_batch(batch, ci, processed, N)
                bi += 1
                ci = 1
            end
        end
    end
end
function _chunk_batches(batches, filt, fg, workers)
    Ncomp = 0
    total_eqs = 0
    unrolled_foreach(filt, batches) do batch
        Ncomp += length(batch)::Int
        total_eqs += length(batch)::Int * _N_eqs(fg, batch)::Int
    end
    chunks = Vector{@NamedTuple{batch_i::Int, comp_i::Int, N::Int}}(undef, workers)

    eqs_per_worker = total_eqs / workers
    bi = 1
    ci = 1
    assigned = 0
    eqs_assigned = 0
    for w in 1:workers
        ci_start = ci
        bi_start = bi
        eqs_in_worker = 0
        assigned_in_worker = 0
        while assigned < Ncomp
            batch = batches[bi]
            filt(batch) || continue

            Neqs = _N_eqs(fg, batch)
            stop_collecting = false
            while true
                if ci > length(batch)
                    break
                end

                # compare, whether adding the new component helps to come closer to eqs_per_worker
                diff_now  = abs(eqs_in_worker - eqs_per_worker)
                diff_next = abs(eqs_in_worker + Neqs - eqs_per_worker)
                stop_collecting = assigned == Ncomp || diff_now ≤ diff_next
                if stop_collecting
                    break
                end

                # add component to worker
                eqs_assigned += Neqs
                eqs_in_worker += Neqs
                assigned_in_worker += 1
                assigned += 1
                ci += 1
            end
            # if the hard stop collection is reached, break, otherwise jump to next batch and continue
            stop_collecting && break

            bi += 1
            ci = 1
        end
        chunk = (; batch_i=bi_start, comp_i=ci_start, N=assigned_in_worker)
        chunks[w] = chunk

        # update eqs per worker estimate for the other workders
        eqs_per_worker = (total_eqs - eqs_assigned) / (workers - w)
    end
    @assert assigned == Ncomp
    return chunks
end
_N_eqs(::Val{:f}, batch) = Int(dim(batch))
_N_eqs(::Val{:g}, batch) = Int(outdim(batch))
_N_eqs(::Val{:fg}, batch) = Int(dim(batch)) + Int(outdim(batch))


@inline function process_batches!(::PolyesterExecution, fg, filt::F, batches, inbufs, duopt) where {F}
    unrolled_foreach(filt, batches) do batch
        (du, u, o, p, t) = duopt
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
