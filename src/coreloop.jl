using Base.Threads: @threads

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
    process_batches!(nw, layer, nw.vertexbatches, dupt, _acc)
    return nothing
end

@unroll function unroll_colorbatches!(nw, layer, colorbatches, dupt, _acc)
    @unroll for cbatch in colorbatches
        process_batches!(nw, layer, cbatch.edgebatches, dupt, _acc)
    end
end

@unroll function process_batches!(nw, layer, batches, dupt, acc)
    @unroll for batch in batches
        process_batch!(nw, layer, batch, dupt, acc)
    end
end

function process_batch!(nw::Network{SequentialExecution}, layer, batch, dupt, acc)
    for idx in 1:length(batch)
        element_kernel!(layer, batch, dupt, acc, idx)
    end
end

function process_batch!(nw::Network{ThreadedExecution}, layer, batch, dupt, acc)
    @threads for idx in 1:length(batch)
        element_kernel!(layer, batch, dupt, acc, idx)
    end
end

@inline function element_kernel!(layer, batch::VertexBatch{F}, dupt, acc, i) where {F}
    @inbounds begin
    du, u, p, t = dupt
    (v_r, acc_r) = vertex_ranges(layer, batch, i)
    pv_r = parameter_range(batch, i)

    vdu = view(du, v_r)
    vu  = view(u, v_r)
    vacc = view(acc, acc_r)
    pv = p===nothing ? nothing : view(p, pv_r)

    batch.fun.f(vdu, vu, vacc, pv, t)
    end
end

@inline function element_kernel!(layer, batch::EdgeBatch{F}, dupt::T, acc, i) where {F<:StaticEdge, T}
    @inbounds begin
    # a cache of size dim per thread
    _, u, p, t = dupt

    # collect all the ranges to index into the data arrays
    (src_r, dst_r, src_acc_r, dst_acc_r, _) = src_dst_ranges(layer, batch, i)
    pe_r = parameter_range(batch, i)

    # create the views into the data & parameters
    vs = view(u, src_r)
    vd = view(u, dst_r)
    pe = p===nothing ? nothing : view(p, pe_r)

    # apply the edge function
    e = batch.fun.f(vs, vd, pe, t)

    apply_accumulation!(coupling(F), layer.accumulator, acc, src_acc_r, dst_acc_r, e)
    end
end

@inline function element_kernel!(layer, batch::EdgeBatch{F}, dupt, acc, i) where {F<:ODEEdge}
    @inbounds begin
    du, u, p, t = dupt
    # for i in 1:length(batch)
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

    apply_accumulation!(coupling(F), layer.accumulator, acc, src_acc_r, dst_acc_r, e)
    end
end

Base.@propagate_inbounds function apply_accumulation!(::AntiSymmetric, f, acc, src_acc_r, dst_acc_r, edge)
    _acc_dst = view(acc, dst_acc_r)
    _acc_src = view(acc, src_acc_r)

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
