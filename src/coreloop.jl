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
function process_vertices!(nw::Network{<:SequentialExecution}, aggbuf, dupt)
    @timeit_debug "execute vertex batches" begin
        unrolled_foreach(nw.vertexbatches) do batch
            (du, u, p, t) = dupt
            for i in 1:length(batch)
                apply_vertex!(batch, i, du, u, aggbuf, p, t)
            end
        end
    end
end

function process_vertices!(nw::Network{<:ThreadedExecution}, aggbuf, dupt)
    @timeit_debug "launch vertex kernels" begin
        _backend = get_backend(dupt[2])
        unrolled_foreach(nw.vertexbatches) do batch
            (du, u, p, t) = dupt
            kernel = vkernel!(_backend)
            kernel(batch, du, u, aggbuf, p, t; ndrange=length(batch))
        end
        KernelAbstractions.synchronize(_backend)
    end
end
@kernel function vkernel!(@Const(batch::VertexBatch{<:ODEVertex}),
                          du, @Const(u),
                          @Const(aggbuf), @Const(p), @Const(t))
    I = @index(Global)
    @inline apply_vertex!(batch, I, du, u, aggbuf, p, t)
    nothing
end

@inline function apply_vertex!(batch::VertexBatch{<:ODEVertex}, i, du, u, aggbuf, p, t)
    _du  = @views du[state_range(batch, i)]
    _u   = @views u[state_range(batch, i)]
    _p   = @views p[parameter_range(batch, i)]
    _agg = @views aggbuf[aggbuf_range(batch, i)]
    batch.fun.f(_du, _u, _agg, _p, t)
    nothing
end


####
#### Edge Layer Execution unbuffered
####
function process_layer!(nw::Network{<:SequentialExecution{false}}, layer, dupt)
    @timeit_debug "execute edge batches" begin
        unrolled_foreach(layer.edgebatches) do batch
            (du, u, p, t) = dupt
            for i in 1:length(batch)
                apply_edge_unbuffered!(batch, i, du, u, nw.im.e_src, nw.im.e_dst, p, t)
            end
        end
    end
end

function process_layer!(nw::Network{<:ThreadedExecution{false}}, layer, dupt)
    @timeit_debug "launch kernels" begin
        _backend = get_backend(dupt[2])
        unrolled_foreach(layer.edgebatches) do batch
            (du, u, p, t) = dupt
            kernel = ekernel!(_backend)
            kernel(batch, du, u, nw.im.e_src, nw.im.e_dst, p, t; ndrange=length(batch))
        end
        KernelAbstractions.synchronize(_backend)
    end
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
    _u   = @views u[state_range(batch, i)]
    _p   = @views p[parameter_range(batch, i)]
    eidx = @views batch.indices[i]
    _src = @views u[srcrange[eidx]]
    _dst = @views u[dstrange[eidx]]
    batch.fun.f(_u, _src, _dst, _p, t)
end

@inline function apply_edge_unbuffered!(batch::EdgeBatch{<:ODEEdge}, i,
                                        du, u, srcrange, dstrange, p, t)
    _du  = @views du[state_range(batch, i)]
    _u   = @views u[state_range(batch, i)]
    _p   = @views p[parameter_range(batch, i)]
    eidx = @views batch.indices[i]
    _src = @views u[srcrange[eidx]]
    _dst = @views u[dstrange[eidx]]
    batch.fun.f(_du, _u, _src, _dst, _p, t)
end

####
#### Edge Layer Execution buffered
####
function process_layer!(nw::Network{<:ThreadedExecution{true}}, layer, dupt)
    # buffered/gathered
    u = dupt[2]
    @timeit_debug "gather" begin
        gbuf = nw.cachepool[u, size(layer.gather_map)]
        NNlib.gather!(gbuf, u, layer.gather_map)
    end

    @timeit_debug "launch kernels" begin
        _backend = get_backend(u)
        unrolled_foreach(layer.edgebatches) do batch
            (du, u, p, t) = dupt
            kernel = ekernel_buffered!(_backend)
            kernel(batch, du, u, gbuf, p, t; ndrange=length(batch))
        end
        KernelAbstractions.synchronize(_backend)
    end
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
    _u   = @views u[state_range(batch, i)]
    _p   = @views p[parameter_range(batch, i)]
    bufr = @views gbuf_range(batch, i)
    _src = @views gbuf[bufr, 1]
    _dst = @views gbuf[bufr, 2]
    batch.fun.f(_u, _src, _dst, _p, t)
    nothing
end

@inline function apply_edge_buffered!(batch::EdgeBatch{<:ODEEdge}, i,
                                      du, u, gbuf, p, t)
    _du  = @views du[state_range(batch, i)]
    _u   = @views u[state_range(batch, i)]
    _p   = @views p[parameter_range(batch, i)]
    bufr = @views gbuf_range(batch, i)
    _src = @views gbuf[bufr, 1]
    _dst = @views gbuf[bufr, 2]
    batch.fun.f(_du, _u, _src, _dst, _p, t)
    nothing
end
