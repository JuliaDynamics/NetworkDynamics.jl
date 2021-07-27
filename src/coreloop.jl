function (nw::Network)(du, u::T, p, t) where {T}
    # for now just one layer
    layer = nw.nl
    dupt = (du, u, p, t)

    # get cache array for accumulator of appropriate type
    _acc = accumulator_cache(layer, T)
    fill!(_acc, zero(eltype(T)))

    # go through all colobatches separatly to avoid writing conflicts to accumulator
    unroll_colorbatches!(nw, layer, layer.colorbatches, dupt, _acc)

    # can be run parallel
    if nw.parallel
        parallell_unroll_batches!(nw, layer, nw.vertexbatches, dupt, _acc)
    else
        unroll_batches!(nw, layer, nw.vertexbatches, dupt, _acc)
    end
end

@unroll function unroll_colorbatches!(nw, layer, colorbatches, dupt, _acc)
    @unroll for cbatch in colorbatches
        if nw.parallel
            parallell_unroll_batches!(nw, layer, cbatch.edgebatches, dupt, _acc)
        else
            unroll_batches!(nw, layer, cbatch.edgebatches, dupt, _acc)
        end
    end
end

@unroll function unroll_batches!(nw, layer, batches, dupt, _acc)
    @unroll for batch in batches
        process_batch!(nw, layer, batch, dupt, _acc)
    end
end

@unroll function parallell_unroll_batches!(nw, layer, batches, dupt, _acc)
    ch = Channel(Inf)
    @unroll for batch in batches
        t = async_process_batch!(nw, layer, batch, dupt, _acc)
        put!(ch, t)
    end
    Base.sync_end(ch)
end

async_process_batch!(nw, layer, batch, dupt, _acc) = Threads.@spawn process_batch!(nw, layer, batch, dupt, _acc)

function process_batch!(nw, layer, batch::VertexBatch{F}, dupt, _acc) where {F}
    @cond_threads nw.parallel for i in 1:length(batch)
        du, u, p, t = dupt
        (v_r, acc_r) = vertex_ranges(layer, batch, i)
        pv_r = parameter_range(batch, i)

        vdu = view(du, v_r)
        vu  = view(u, v_r)
        acc = view(_acc, acc_r)
        pv = p===nothing ? nothing : view(p, pv_r)

        batch.fun.f(vdu, vu, acc, pv, t)
    end
end

function process_batch!(nw, layer, batch::EdgeBatch{F}, dupt::T, _acc) where {F<:StaticEdge, T}
    # a cache of size dim per thread
    dim = batch.dim
    _, u, p, t = dupt
    cachesize = Threads.nthreads() * dim
    # XXX: WTF why allocations for this line below?
    # cachesize = nw.parallel ? Threads.nthreads() * dim : dim
    _cache = getcache(layer.cachepool, typeof(u), cachesize)

    @cond_threads nw.parallel for i in 1:length(batch)
        # each thread should get it's portion of the cache
        cidx = 1 + (Threads.threadid() - 1) * dim
        _c = view(_cache, cidx:cidx+dim-1)

        # collect all the ranges to index into the data arrays
        (src_r, dst_r, src_acc_r, dst_acc_r, e_r) = src_dst_ranges(layer, batch, i)
        pe_r = parameter_range(batch, i)

        # create the views into the data & parameters
        vs = view(u, src_r)
        vd = view(u, dst_r)
        pe = p===nothing ? nothing : view(p, pe_r)

        # apply the edge function
        batch.fun.f(_c, vs, vd, pe, t)

        apply_accumulation!(coupling(F), layer.accumulator, _acc, src_acc_r, dst_acc_r, _c)
    end
end

function process_batch!(nw, layer, batch::EdgeBatch{F}, dupt, _acc) where {F<:ODEEdge}
    du, u, p, t = dupt

    @cond_threads nw.parallel for i in 1:length(batch)
        # collect all the ranges to index into the data arrays
        (src_r, dst_r, src_acc_r, dst_acc_r, e_r) = src_dst_ranges(layer, batch, i)
        pe_r = parameter_range(batch, i)

        # create the views into the data & parameters
        de = view(du, e_r)
        e  = view(u, e_r)
        vs = view(u, src_r)
        vd = view(u, dst_r)
        pe = p===nothing ? nothing : view(p, pe_r)

        # apply the edge function
        batch.fun.f(de, e, vs, vd, pe, t)

        apply_accumulation!(coupling(F), layer.accumulator, _acc, src_acc_r, dst_acc_r, e)
    end
end

Base.@propagate_inbounds function apply_accumulation!(::AntiSymmetric, f, _acc, src_acc_r, dst_acc_r, edge)
    _acc_dst = view(_acc, dst_acc_r)
    _acc_src = view(_acc, src_acc_r)

    for i in 1:length(dst_acc_r)
        _acc_dst[i] = f(_acc_dst[i],  edge[i])
        _acc_src[i] = f(_acc_src[i], -edge[i])
    end
end

Base.@propagate_inbounds function apply_accumulation!(::Symmetric, f, _acc, src_acc_r, dst_acc_r, edge)
    _acc_dst = view(_acc, dst_acc_r)
    _acc_src = view(_acc, src_acc_r)

    for i in 1:length(dst_acc_r)
        _acc_dst[i] = f(_acc_dst[i], edge[i])
        _acc_src[i] = f(_acc_src[i], edge[i])
    end
end

Base.@propagate_inbounds function apply_accumulation!(::Directed, f, _acc, _, dst_acc_r, edge)
    _acc_dst = view(_acc, dst_acc_r)

    for i in 1:length(dst_acc_r)
        _acc_dst[i] = f(_acc_dst[i], edge[i])
    end
end

Base.@propagate_inbounds function apply_accumulation!(::Fiducial, f, _acc, src_acc_r, dst_acc_r, edge)
    _acc_dst = view(_acc, dst_acc_r)
    _acc_src = view(_acc, src_acc_r)
    offset = Int(length(edge) / 2)

    for i in 1:length(dst_acc_r)
        _acc_dst[i] = f(_acc_dst[i], edge[i])
        _acc_src[i] = f(_acc_src[i], edge[offset + i])
    end
end

