function (nw::Network)(du, u::T, p, t) where {T}
    # resize!(u, length(u) + im.size_static)
    # s = nw.cachepool[u, im.size_dynamic + im.size_static]
    # s[1:im.size_dynamic] .= u
    s = nw.cachepool[u, im.size_static]

    dupt = (du, u, s, p, t)


    # NOTE: first all static batches, than all dynamic batches
    # maybe disallow static vertices entierly, otherwise the order gets complicated

    # execute all edgebatches
    # for now just one layer
    layer = nw.nl
    process_layer!(nw, layer)

    # can be run parallel
    process_batches!(nw, layer, nw.vertexbatches, dupt)
    return nothing
end

function process_layer!(nw::Network{}, layer, dupt)
    process_batches!(nw, layer, layer.edgebatchesbatch, dupt)
end

@unroll function process_batches!(nw, layer, batches, dupt)
    @unroll for batch in batches
        process_batch!(nw, layer, batch, dupt)
    end
end

function process_batch!(nw::Network{SequentialExecution}, layer, batch, dupt)
    for idx in 1:length(batch)
        element_kernel!(layer, batch, dupt, idx)
    end
end

function process_batch!(nw::Network{ThreadedExecution}, layer, batch, dupt)
    @threads for idx in 1:length(batch)
        element_kernel!(layer, batch, dupt, idx)
    end
end


@inline function element_kernel!(layer, batch::EdgeBatch{F}, dupt::T, acc, i) where {F<:StaticEdge, T}
    @inbounds begin
    # a cache of size dim per thread
    _, u, s, p, t = dupt

    # collect all the ranges to index into the data arrays
    (src_r, dst_r, src_acc_r, dst_acc_r, _) = src_dst_ranges(layer, batch, i)
    pe_r = parameter_range(batch, i)

    # create the views into the data & parameters
    vs = view(u, src_r)
    vd = view(u, dst_r)
    pe = p==SciMLBase.NullParameters() ? p : view(p, pe_r)

    # apply the edge function
    e = batch.fun.f(vs, vd, pe, t)

    apply_accumulation!(coupling(F), layer.accumulator, acc, src_acc_r, dst_acc_r, e)
    end
end

@inline function element_kernel!(layer, batch::EdgeBatch{F}, dupt, acc, i) where {F<:ODEEdge}
    @inbounds begin
    du, u, s, p, t = dupt
    # for i in 1:length(batch)
    # collect all the ranges to index into the data arrays
    (src_r, dst_r, src_acc_r, dst_acc_r, e_r) = src_dst_ranges(layer, batch, i)
    pe_r = parameter_range(batch, i)

    # create the views into the data & parameters
    de = view(du, e_r)
    e  = view(u, e_r)
    vs = view(u, src_r)
    vd = view(u, dst_r)
    pe = p==SciMLBase.NullParameters() ? p : view(p, pe_r)

    # apply the edge function
    batch.fun.f(de, e, vs, vd, pe, t)

    apply_accumulation!(coupling(F), layer.accumulator, acc, src_acc_r, dst_acc_r, e)
    end
end

@inline function element_kernel!(layer, batch::VertexBatch{F}, dupt, acc, i) where {F}
    @inbounds begin
    du, u, s, p, t = dupt

    s_r = state_range(batch, i)
    p_r = parameter_range(batch, i)

    vdu = view(du, s_r)
    vu  = view(u, s_r)

    vacc = view(acc, acc_r)

    pv = p==SciMLBase.NullParameters() ? p : view(p, pv_r)

    batch.fun.f(vdu, vu, vacc, pv, t)
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
