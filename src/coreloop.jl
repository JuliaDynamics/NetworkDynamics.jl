function (nw::NetworkDynamic)(du, u::T, p, t) where {T}
    # for now just one layer
    layer = nw.nl

    # get cache array for accumulator of appropriate type
    _acc = accumulator_cache(layer, T)
    fill!(_acc, zero(T))

    # go through all colobatches separatly to avoid writing conflicts to accumulator
    for colorbatch in layer.colorbatches
        # inner loop can be parallel
        process_colorbatch!(du, u, p, t, _acc, layer, colorbatch)
    end

    # can be run parallel
    for vertexbatch in nw.vertexbatches
        process_vertexbatch!(du, u, p, t, _acc, layer, vertexbatch)
    end

    # calculate vertex functions (needs  fill!(du, 0) ? )
    @inbounds for vidx in 1:nv(nw)
        vdu = view(du, :, vidx)
        vu  = view(u, :, vidx)
        agg = view(_aggregation, :, vidx)
        nw.vfun(vdu, vu, agg, p, t)
    end
end

function process_colorbatch!(du, u, p, t, _acc, layer, colorbatch::ColorBatch)
    for batch in colorbatch.edgebatches
        process_edgebatch!(_acc, du, u, p, t, batch, layer)
    end
end

function process_edgebatch!(_, u::T, p, t, _acc, layer, batch::EdgeBatch{F}) where {T,F<:StaticEdge}
    # a cache of size dim per thread
    dim = batch.dim
    _cache =  getcache(layer.cachepool, T, Threads.nthreads() * dim)

    @inbounds for i in 1:length(batch)
        # each thread should get it's portion of the cache
        cidx = (Threads.threadid() - 1) * dim
        _c = view(_cache, cidx:cidx+dim)

        # collect all the ranges to index into the data arrays
        (src_r, dst_r, src_acc_r, dst_acc_r) = src_dst_ranges(layer, batch, i)
        pe_r = parameter_range(batch, i)

        # create the views into the data & parameters
        vs = view(u, src_r)
        vd = view(u, dst_r)
        pe = view(p, pe_r)

        # apply the edge function
        batch.fun(_c, vs, vd, pe, t)

        # apply the accumulator
        # XXX: move this to function and dispatch based on couplingtype of batch?
        _c_part = dim == layer.accdim ? _c : view(_c, 1:layer.accdim)
        _acc[dst_acc_r] .=   layer.accumulator.(_acc[dst_acc_r], _c_part)
        _acc[src_acc_r] .= - layer.accumulator.(_acc[src_acc_r], _c_part)
    end
end

function process_vertexbatch!(du, u, p, t, _acc, layer, batch::VertexBatch{F}) where {F}
    @inbounds for i in 1:length(batch)
        (v_r, acc_r) = vertex_ranges(layer, batch, i)
        pv_r = parameter_range(batch, i)

        vdu = view(du, v_r)
        vu  = view(u, v_r)
        acc = view(_acc, acc_r)
        pv  = view(p, pv_r)

        batch.fun(vdu, vu, acc, pv, t)
    end
end
