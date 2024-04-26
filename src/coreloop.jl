function (nw::Network)(du, u::T, p, t) where {T}
    # @timeit_debug "coreloop" begin
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

    # end
    return nothing
end

function process_vertices!(nw::Network{<:ThreadedExecution}, aggbuf, dupt)
    (du, u, p, t) = dupt
    @timeit_debug "launch kernels" begin
        _backend = get_backend(du)
        unrolled_foreach(nw.vertexbatches) do batch
            kernel = vkernel!(_backend)
            kernel(batch, du, u, aggbuf, p, t; ndrange=length(batch))
        end
        KernelAbstractions.synchronize(_backend)
    end
end
@kernel function vkernel!(@Const(batch::VertexBatch{<:ODEVertex}), du, @Const(u),
                          @Const(aggbuf), @Const(p), @Const(t))
    I = @index(Global)
    @views begin
        _du  = du[state_range(batch, I)]
        _u   = u[state_range(batch, I)]
        _p   = p[parameter_range(batch, I)]
        _agg = aggbuf[aggbuf_range(batch, I)]
    end
    batch.fun.f(_du, _u, _agg, _p, t)
end

function process_layer!(nw::Network{<:ThreadedExecution{false}}, layer, dupt)
    @timeit_debug "launch kernels" begin
        (du, u, p, t) = dupt
        _backend = get_backend(du)
        unrolled_foreach(layer.edgebatches) do batch
            kernel = ekernel!(_backend)
            kernel(batch, du, u, nw.im.e_src, nw.im.e_dst, p, t; ndrange=length(batch))
        end
        KernelAbstractions.synchronize(_backend)
    end
end
@kernel function ekernel!(@Const(batch::EdgeBatch{<:StaticEdge}),
                          @Const(du), u,
                          @Const(srcrange), @Const(dstrange), @Const(p), @Const(t))
    I = @index(Global)
    @views begin
        _u   = u[state_range(batch, I)]
        _p   = p[parameter_range(batch, I)]
        eidx = batch.indices[I]
        _src = u[srcrange[eidx]]
        _dst = u[dstrange[eidx]]
    end
    batch.fun.f(_u, _src, _dst, _p, t)
end
@kernel function ekernel!(@Const(batch::EdgeBatch{<:ODEEdge}),
                          du, @Const(u),
                          @Const(srcrange), @Const(dstrange), @Const(p), @Const(t))
    I = @index(Global)
    @views begin
        _du  = du[state_range(batch, I)]
        _u   = u[state_range(batch, I)]
        _p   = p[parameter_range(batch, I)]
        eidx = batch.indices[I]
        _src = u[srcrange[eidx]]
        _dst = u[dstrange[eidx]]
    end
    batch.fun.f(_du, _u, _src, _dst, _p, t)
end

function process_layer!(nw::Network{<:ThreadedExecution{true}}, layer, dupt)
    # buffered/gathered
    (du, u, p, t) = dupt
    @timeit_debug "gather" begin
        gbuf = nw.cachepool[u, size(layer.gather_map)]
        NNlib.gather!(gbuf, u, layer.gather_map)
    end

    @timeit_debug "launch kernels" begin
        _backend = get_backend(du)
        unrolled_foreach(layer.edgebatches) do batch
            kernel = ekernel_buffered!(_backend)
            kernel(batch, du, u, gbuf, p, t; ndrange=length(batch))
        end
        KernelAbstractions.synchronize(_backend)
    end
end
@kernel function ekernel_buffered!(@Const(batch::EdgeBatch{<:StaticEdge}), @Const(du), u,
                                   @Const(gbuf), @Const(p), @Const(t))
    I = @index(Global)
    @views begin
        _u   = u[state_range(batch, I)]
        _p   = p[parameter_range(batch, I)]
        bufr = gbuf_range(batch, I)
        _src = gbuf[bufr, 1]
        _dst = gbuf[bufr, 2]
    end
    batch.fun.f(_u, _src, _dst, _p, t)
end
@kernel function ekernel_buffered!(@Const(batch::EdgeBatch{<:ODEEdge}), du, @Const(u),
                                   @Const(gbuf), @Const(p), @Const(t))
    I = @index(Global)
    @views begin
        _du  = du[state_range(batch, I)]
        _u   = u[state_range(batch, I)]
        _p   = p[parameter_range(batch, I)]
        bufr = gbuf_range(batch, I)
        _src = gbuf[bufr, 1]
        _dst = gbuf[bufr, 2]
    end
    batch.fun.f(_du, _u, _src, _dst, _p, t)
end
