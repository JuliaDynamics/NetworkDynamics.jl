function (nw::Network)(du, u::T, p, t) where {T}
    @timeit_debug "coreloop" begin
        @timeit_debug "fill zeros" begin
            fill!(du, zero(eltype(du)))
        end
        @timeit_debug "create _u" begin
            _u = nw.cachepool[u, nw.im.lastidx_static]
            _u[1:nw.im.lastidx_dynamic] .= u
        end

        dupt = (du, _u, p, t)

        # NOTE: first all static vertices, than all edges (regarless) then dyn vertices
        # maybe disallow static vertices entierly, otherwise the order gets complicated

        @timeit_debug "process layer" process_layer!(nw, nw.layer, dupt)

        @timeit_debug "aggregate" begin
            aggbuf = nw.cachepool[_u, nw.im.lastidx_aggr]
            aggregate!(nw.layer.aggregator, aggbuf, _u)
        end

        @timeit_debug "process vertices" process_vertices!(nw, aggbuf, dupt)
    end
    return nothing
end


####
#### Vertex Execution
####
@inline function process_vertices!(nw::Network{<:SequentialExecution}, aggbuf, dupt)
    unrolled_foreach(nw.vertexbatches) do batch
        (du, u, p, t) = dupt
        for i in 1:length(batch)
            apply_vertex!(batch, i, du, u, aggbuf, p, t)
        end
    end
end

@inline function process_vertices!(nw::Network{<:ThreadedExecution}, aggbuf, dupt)
    _backend = get_backend(dupt[2])
    unrolled_foreach(nw.vertexbatches) do batch
        (du, u, p, t) = dupt
        kernel = vkernel!(_backend)
        kernel(batch, du, u, aggbuf, p, t; ndrange=length(batch))
        # kernel = vkernel!(_backend, 32, length(batch))
        # kernel(batch, du, u, aggbuf, p, t)
    end
    KernelAbstractions.synchronize(_backend)
end
@kernel function vkernel!(@Const(batch::VertexBatch{<:ODEVertex}),
                          du, @Const(u),
                          @Const(aggbuf), @Const(p), @Const(t))
    I = @index(Global)
    @inline apply_vertex!(batch, I, du, u, aggbuf, p, t)
    nothing
end

@inline function apply_vertex!(batch::VertexBatch{<:ODEVertex}, i, du, u, aggbuf, p, t)
    @inbounds begin
        _du  = @views du[state_range(batch, i)]
        _u   = @views u[state_range(batch, i)]
        _p   = isnothing(p) ? p : @views p[parameter_range(batch, i)]
        _agg = @views aggbuf[aggbuf_range(batch, i)]
        batch.fun.f(_du, _u, _agg, _p, t)
    end
    nothing
end


####
#### Edge Layer Execution unbuffered
####
@inline function process_layer!(nw::Network{<:SequentialExecution{false}}, layer, dupt)
    unrolled_foreach(layer.edgebatches) do batch
        (du, u, p, t) = dupt
        for i in 1:length(batch)
            apply_edge_unbuffered!(batch, i, du, u, nw.im.e_src, nw.im.e_dst, p, t)
        end
    end
end

@inline function process_layer!(nw::Network{<:ThreadedExecution{false}}, layer, dupt)
    _backend = get_backend(dupt[2])
    unrolled_foreach(layer.edgebatches) do batch
        (du, u, p, t) = dupt
        kernel = ekernel!(_backend)
        kernel(batch, du, u, nw.im.e_src, nw.im.e_dst, p, t; ndrange=length(batch))
        # kernel = ekernel!(_backend, 32, length(batch))
        # kernel(batch, du, u, nw.im.e_src, nw.im.e_dst, p, t)
    end
    KernelAbstractions.synchronize(_backend)
end
@kernel function ekernel!(@Const(batch::EdgeBatch{<:StaticEdge}),
                          @Const(du), u,
                          @Const(srcrange), @Const(dstrange),
                          @Const(p), @Const(t))
    I = @index(Global)
    @inline apply_edge_unbuffered!(batch, I, du, u, srcrange, dstrange, p, t)
end
@kernel function ekernel!(@Const(batch::EdgeBatch{<:ODEEdge}),
                          du, @Const(u),
                          @Const(srcrange), @Const(dstrange),
                          @Const(p), @Const(t))
    I = @index(Global)
    @inline apply_edge_unbuffered!(batch, I, du, u, srcrange, dstrange, p, t)
end

@inline function apply_edge_unbuffered!(batch::EdgeBatch{<:StaticEdge}, i,
                                        du, u, srcrange, dstrange, p, t)
    @inbounds begin
        _u   = @views u[state_range(batch, i)]
        _p   = isnothing(p) ? p : @views p[parameter_range(batch, i)]
        eidx = @views batch.indices[i]
        _src = @views u[srcrange[eidx]]
        _dst = @views u[dstrange[eidx]]
        batch.comp.f(_u, _src, _dst, _p, t)
    end
    nothing
end

@inline function apply_edge_unbuffered!(batch::EdgeBatch{<:ODEEdge}, i,
                                        du, u, srcrange, dstrange, p, t)
    @inbounds begin
        _du  = @views du[state_range(batch, i)]
        _u   = @views u[state_range(batch, i)]
        _p   = isnothing(p) ? p : @views p[parameter_range(batch, i)]
        eidx = @views batch.indices[i]
        _src = @views u[srcrange[eidx]]
        _dst = @views u[dstrange[eidx]]
        batch.comp.f(_du, _u, _src, _dst, _p, t)
    end
    nothing
end

####
#### Edge Layer Execution buffered
####
@inline function process_layer!(nw::Network{<:SequentialExecution{true}}, layer, dupt)
    u = dupt[2]
    gbuf = nw.cachepool[u, size(layer.gather_map)]
    NNlib.gather!(gbuf, u, layer.gather_map)

    unrolled_foreach(layer.edgebatches) do batch
        (_du, _u, _p, _t) = dupt
        for i in 1:length(batch)
            apply_edge_buffered!(batch, i, _du, _u, gbuf, _p, _t)
        end
    end
end

@inline function process_layer!(nw::Network{<:ThreadedExecution{true}}, layer, dupt)
    # buffered/gathered
    u = dupt[2]
    gbuf = nw.cachepool[u, size(layer.gather_map)]
    NNlib.gather!(gbuf, u, layer.gather_map)

    backend = get_backend(u)
    unrolled_foreach(layer.edgebatches) do batch
        (_du, _u, _p, _t) = dupt
        kernel = ekernel_buffered!(backend)
        kernel(batch, _du, _u, gbuf, _p, _t; ndrange=length(batch))
    end
    KernelAbstractions.synchronize(backend)
end
@kernel function ekernel_buffered!(@Const(batch::EdgeBatch{<:StaticEdge}),
                                   @Const(du), u,
                                   @Const(gbuf), @Const(p), @Const(t))
    I = @index(Global)
    apply_edge_buffered!(batch, I, du, u, gbuf, p, t)
end
@kernel function ekernel_buffered!(@Const(batch::EdgeBatch{<:ODEEdge}),
                                   du, @Const(u),
                                   @Const(gbuf), @Const(p), @Const(t))
    I = @index(Global)
    apply_edge_buffered!(batch, I, du, u, gbuf, p, t)
end

@inline function apply_edge_buffered!(batch::EdgeBatch{<:StaticEdge}, i,
                                      du, u, gbuf, p, t)
    @inbounds begin
        _u   = @views u[state_range(batch, i)]
        _p   = isnothing(p) ? p : @views p[parameter_range(batch, i)]
        bufr = @views gbuf_range(batch, i)
        _src = @views gbuf[bufr, 1]
        _dst = @views gbuf[bufr, 2]
        batch.comp.f(_u, _src, _dst, _p, t)
    end
    nothing
end

@inline function apply_edge_buffered!(batch::EdgeBatch{<:ODEEdge}, i,
                                      du, u, gbuf, p, t)
    @inbounds begin
        _du  = @views du[state_range(batch, i)]
        _u   = @views u[state_range(batch, i)]
        _p   = isnothing(p) ? p : @views p[parameter_range(batch, i)]
        bufr = @views gbuf_range(batch, i)
        _src = @views gbuf[bufr, 1]
        _dst = @views gbuf[bufr, 2]
        batch.comp.f(_du, _u, _src, _dst, _p, t)
    end
    nothing
end
