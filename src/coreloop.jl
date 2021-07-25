function (nw::NetworkDynamic)(du, u::T, p, t) where {T}
    # for now just one layer
    layer = nw.nl

    # get cache array for accumulator of appropriate type
    _acc = accumulator_cache(layer, T)
    fill!(_acc, zero(eltype(T)))

    # go through all colobatches separatly to avoid writing conflicts to accumulator
    for colorbatch in layer.colorbatches
        # inner loop can be parallel
        process_colorbatch!(du, u, p, t, _acc, layer, colorbatch)
    end

    # can be run parallel
    for vertexbatch in nw.vertexbatches
        process_vertexbatch!(du, u, p, t, _acc, layer, vertexbatch)
    end

end

function process_colorbatch!(du, u, p, t, _acc, layer, colorbatch::ColorBatch)
    for batch in colorbatch.edgebatches
        process_edgebatch!(du, u, p, t, _acc, layer, batch)
    end
end

function process_edgebatch!(_, u::T, p, t, _acc, layer, batch::EdgeBatch{F}) where {T,F<:StaticEdge}
    # a cache of size dim per thread
    dim = batch.dim
    _cache =  getcache(layer.cachepool, T, Threads.nthreads() * dim)

    @inbounds for i in 1:length(batch)
        # each thread should get it's portion of the cache
        cidx = (Threads.threadid() - 1) * dim
        _c = view(_cache, cidx:cidx+dim-1)

        # collect all the ranges to index into the data arrays
        (src_r, dst_r, src_acc_r, dst_acc_r, e_r) = src_dst_ranges(layer, batch, i)
        pe_r = parameter_range(batch, i)

        # create the views into the data & parameters
        vs = view(u, src_r)
        vd = view(u, dst_r)
        pe = view(p, pe_r)

        # apply the edge function
        batch.fun.f(_c, vs, vd, pe, t)

        # apply the accumulator
        # XXX: move this to function and dispatch based on couplingtype of batch?
        _c_part = dim == layer.accdim ? _c : view(_c, 1:layer.accdim)
        _acc_src = view(_acc, src_acc_r)
        _acc_dst = view(_acc, dst_acc_r)

        for i in 1:layer.accdim
            _acc_dst[i] = layer.accumulator(_acc_src[i],  _c_part[i])
            _acc_src[i] = layer.accumulator(_acc_src[i], -_c_part[i])
        end
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

        batch.fun.f(vdu, vu, acc, pv, t)
    end
end
