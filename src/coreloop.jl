function (nw::Network)(du, u::T, p, t) where {T}
    # for now just one layer
    layer = nw.nl

    # get cache array for accumulator of appropriate type
    _acc = accumulator_cache(layer, T)
    fill!(_acc, zero(eltype(T)))

    # go through all colobatches separatly to avoid writing conflicts to accumulator
    unroll_colorbatches!(du, u, p, t, _acc, layer)

    # can be run parallel
    unroll_vertexbatches!(du, u, p, t, _acc, layer, nw.vertexbatches)
end

@unroll function unroll_colorbatches!(du, u, p, t, _acc, layer, cbs=layer.colorbatches)
    @unroll for colorbatch in cbs
        unroll_edgebatches!(du, u, p, t, _acc, layer, colorbatch)
    end
end

@unroll function unroll_edgebatches!(du, u, p, t, _acc, layer, cb::ColorBatch, ebs=cb.edgebatches)
    @unroll for batch in ebs
        process_edgebatch!(du, u, p, t, _acc, layer, batch)
    end
end

function process_edgebatch!(_, u::T, p, t, _acc, layer, batch::EdgeBatch{F}) where {T,F<:StaticEdge}
    # a cache of size dim per thread
    dim = batch.dim
    _cache =  getcache(layer.cachepool, T, Threads.nthreads() * dim)

    @inbounds for i in 1:length(batch)
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

@unroll function unroll_vertexbatches!(du, u, p, t, _acc, layer, vertexbatches)
    @unroll for batch in vertexbatches
        process_vertexbatch!(du, u, p, t, _acc, layer, batch)
    end
end

function process_vertexbatch!(du, u, p, t, _acc, layer, batch::VertexBatch{F}) where {F}
    @inbounds for i in 1:length(batch)
        (v_r, acc_r) = vertex_ranges(layer, batch, i)
        pv_r = parameter_range(batch, i)

        vdu = view(du, v_r)
        vu  = view(u, v_r)
        acc = view(_acc, acc_r)
        pv = p===nothing ? nothing : view(p, pv_r)

        batch.fun.f(vdu, vu, acc, pv, t)
    end
end
