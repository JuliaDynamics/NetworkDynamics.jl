function (nw::Network{A,B,C,D,E})(du::dT, u::T, p, t) where {A,B,C,D,E,dT,T}
    if !(eachindex(du) == eachindex(u) == 1:nw.im.lastidx_dynamic)
        throw(ArgumentError("du or u does not have expected size $(nw.im.lastidx_dynamic)"))
    end
    if nw.im.lastidx_p > 0 && _indexable(p) && !(eachindex(p) == 1:nw.im.lastidx_p)
        throw(ArgumentError("p does not has expecte size $(nw.im.lastidx_p)"))
    end
    ex = executionstyle(nw)
    @timeit_debug "coreloop" begin
        @timeit_debug "fill zeros" begin
            fill!(du, zero(eltype(du)))
        end
        @timeit_debug "create _u" begin
            _u = nw.cachepool[du, nw.im.lastidx_static]::dT
            _u[1:nw.im.lastidx_dynamic] .= u
        end

        s = _u # currently defaults to stacked vector
        dupt = (du, _u, s, p, t)

        # NOTE: first all static vertices, than all edges (regarless) then dyn vertices
        # maybe disallow static vertices entierly, otherwise the order gets complicated

        @timeit_debug "process layer" process_layer!(ex, nw, nw.layer, dupt)

        @timeit_debug "aggregate" begin
            if nw.im.lastidx_aggr == nw.im.lastidx_static
                error("Aggbuf and _u buf cannot be the same size! This is a known bug.")
            end
            aggbuf = nw.cachepool[_u, nw.im.lastidx_aggr]::T
            aggregate!(nw.layer.aggregator, aggbuf, _u)
        end

        @timeit_debug "process vertices" process_vertices!(ex, nw, aggbuf, dupt)
    end
    return nothing
end

# current compat, no diff between static and dynamic range
dynamic_range(batch, i) = state_range(batch, i)
static_range(batch, i) = state_range(batch, i)

####
#### Vertex Execution
####
@inline function process_vertices!(::SequentialExecution, nw, aggbuf, dupt)
    unrolled_foreach(nw.vertexbatches) do batch
        (du, u, s, p, t) = dupt
        for i in 1:length(batch)
            _type = comptype(batch)
            _batch = essence(batch)
            apply_vertex!(_type, _batch, i, du, u, s, aggbuf, p, t)
        end
    end
end

@inline function process_vertices!(::ThreadedExecution, nw, aggbuf, dupt)
    unrolled_foreach(nw.vertexbatches) do batch
        (du, u, s, p, t) = dupt
        Threads.@threads for i in 1:length(batch)
            _type = comptype(batch)
            _batch = essence(batch)
            apply_vertex!(_type, _batch, i, du, u, s, aggbuf, p, t)
        end
    end
end

@inline function process_vertices!(::PolyesterExecution, nw, aggbuf, dupt)
    unrolled_foreach(nw.vertexbatches) do batch
        (du, u, s, p, t) = dupt
        Polyester.@batch for i in 1:length(batch)
            _type = comptype(batch)
            _batch = essence(batch)
            apply_vertex!(_type, _batch, i, du, u, s, aggbuf, p, t)
        end
    end
end

@inline function process_vertices!(::KAExecution, nw, aggbuf, dupt)
    _backend = get_backend(dupt[2])
    unrolled_foreach(nw.vertexbatches) do batch
        (du, u, s, p, t) = dupt
        kernel = vkernel!(_backend)
        kernel(comptype(batch), essence(batch),
               du, u, s, aggbuf, p, t; ndrange=length(batch))
    end
    KernelAbstractions.synchronize(_backend)
end
@kernel function vkernel!(::Type{T}, @Const(batch),
                          du, @Const(u), s,
                          @Const(aggbuf), @Const(p), @Const(t)) where {T}
    I = @index(Global)
    apply_vertex!(T, batch, I, du, u, s, aggbuf, p, t)
    nothing
end

@inline function apply_vertex!(::Type{T}, batch, i, du, u, s, aggbuf, p, t) where {T}
    @inbounds begin
        _du  = _has_dynamic(T) ? view(du, dynamic_range(batch, i))   : nothing
        _u   = _has_dynamic(T) ? view(u,  dynamic_range(batch, i))   : nothing
        _s   = _has_static(T)  ? view(s,  static_range(batch, i))    : nothing
        _p   = _indexable(p)   ? view(p,  parameter_range(batch, i)) : p
        _agg = view(aggbuf, aggbuf_range(batch, i))
        apply_compf(T, compf(batch), _du, _u, _s, _agg, _p, t)
    end
    nothing
end

@propagate_inbounds function apply_compf(::Type{<:ODEVertex}, f::F, du, u, s, agg, p, t) where {F}
    f(du, u, agg, p, t)
end
@propagate_inbounds function apply_compf(::Type{<:StaticVertex}, f::F, du, u, s, agg, p, t) where {F}
    f(s, agg, p, t)
end

####
#### Edge Layer Execution unbuffered
####
@inline function process_layer!(::SequentialExecution{false}, nw, layer, dupt)
    unrolled_foreach(layer.edgebatches) do batch
        (du, u, s, p, t) = dupt
        for i in 1:length(batch)
            _type = comptype(batch)
            _batch = essence(batch)
            apply_edge_unbuffered!(_type, _batch, i, du, u, s, nw.im.e_src, nw.im.e_dst, p, t)
        end
    end
end

@inline function process_layer!(::ThreadedExecution{false}, nw, layer, dupt)
    unrolled_foreach(layer.edgebatches) do batch
        (du, u, s, p, t) = dupt
        Threads.@threads for i in 1:length(batch)
            _type = comptype(batch)
            _batch = essence(batch)
            apply_edge_unbuffered!(_type, _batch, i, du, u, s, nw.im.e_src, nw.im.e_dst, p, t)
        end
    end
end

@inline function process_layer!(::PolyesterExecution{false}, nw, layer, dupt)
    unrolled_foreach(layer.edgebatches) do batch
        (du, u, s, p, t) = dupt
        Polyester.@batch for i in 1:length(batch)
            _type = comptype(batch)
            _batch = essence(batch)
            apply_edge_unbuffered!(_type, _batch, i, du, u, s, nw.im.e_src, nw.im.e_dst, p, t)
        end
    end
end

@inline function process_layer!(::KAExecution{false}, nw, layer, dupt)
    _backend = get_backend(dupt[2])
    unrolled_foreach(layer.edgebatches) do batch
        (du, u, s, p, t) = dupt
        kernel = ekernel!(_backend)
        kernel(comptype(batch), essence(batch),
               du, u, s, nw.im.e_src, nw.im.e_dst, p, t; ndrange=length(batch))
    end
    KernelAbstractions.synchronize(_backend)
end
@kernel function ekernel!(::Type{T}, @Const(batch),
                          du, @Const(u), s,
                          @Const(srcrange), @Const(dstrange),
                          @Const(p), @Const(t)) where {T}
    I = @index(Global)
    apply_edge_unbuffered!(T, batch, I, du, u, s, srcrange, dstrange, p, t)
end

@inline function apply_edge_unbuffered!(::Type{T}, batch, i,
                                        du, u, s, srcrange, dstrange, p, t) where {T}
    @inbounds begin
        _du  = _has_dynamic(T) ? view(du, dynamic_range(batch, i))   : nothing
        _u   = _has_dynamic(T) ? view(u,  dynamic_range(batch, i))   : nothing
        _s   = _has_static(T)  ? view(s,  static_range(batch, i))    : nothing
        _p   = _indexable(p)   ? view(p,  parameter_range(batch, i)) : p
        eidx = @views batch.indices[i]
        _src = @views u[srcrange[eidx]]
        _dst = @views u[dstrange[eidx]]
        apply_compf(T, compf(batch), _du, _u, _s, _src, _dst, _p, t)
    end
    nothing
end

####
#### Edge Layer Execution buffered
####
@inline function process_layer!(::SequentialExecution{true}, nw, layer, dupt)
    u = dupt[2]
    gbuf = nw.cachepool[u, size(layer.gather_map)]
    NNlib.gather!(gbuf, u, layer.gather_map)

    unrolled_foreach(layer.edgebatches) do batch
        (_du, _u, _s, _p, _t) = dupt
        for i in 1:length(batch)
            _type = comptype(batch)
            _batch = essence(batch)
            apply_edge_buffered!(_type, _batch, i, _du, _u, _s, gbuf, _p, _t)
        end
    end
end

@inline function process_layer!(::ThreadedExecution{true}, nw, layer, dupt)
    u = dupt[2]
    gbuf = nw.cachepool[u, size(layer.gather_map)]
    NNlib.gather!(gbuf, u, layer.gather_map)

    unrolled_foreach(layer.edgebatches) do batch
        (_du, _u, _s, _p, _t) = dupt
        Threads.@threads for i in 1:length(batch)
            _type = comptype(batch)
            _batch = essence(batch)
            apply_edge_buffered!(_type, _batch, i, _du, _u, _s, gbuf, _p, _t)
        end
    end
end

@inline function process_layer!(::PolyesterExecution{true}, nw, layer, dupt)
    u = dupt[2]
    gbuf = nw.cachepool[u, size(layer.gather_map)]
    NNlib.gather!(gbuf, u, layer.gather_map)

    unrolled_foreach(layer.edgebatches) do batch
        (_du, _u, _s, _p, _t) = dupt
        Polyester.@batch for i in 1:length(batch)
            _type = comptype(batch)
            _batch = essence(batch)
            apply_edge_buffered!(_type, _batch, i, _du, _u, _s, gbuf, _p, _t)
        end
    end
end

@inline function process_layer!(::KAExecution{true}, nw, layer, dupt)
    # buffered/gathered
    u = dupt[2]
    gbuf = nw.cachepool[u, size(layer.gather_map)]
    NNlib.gather!(gbuf, u, layer.gather_map)

    backend = get_backend(u)
    unrolled_foreach(layer.edgebatches) do batch
        (_du, _u, _s, _p, _t) = dupt
        kernel = ekernel_buffered!(backend)
        kernel(comptype(batch), essence(batch),
               _du, _u, _s, gbuf, _p, _t; ndrange=length(batch))
    end
    KernelAbstractions.synchronize(backend)
end

@kernel function ekernel_buffered!(::Type{T}, @Const(batch),
                                   du, @Const(u), s,
                                   @Const(gbuf), @Const(p), @Const(t)) where {T}
    I = @index(Global)
    apply_edge_buffered!(T, batch, I, du, u, s, gbuf, p, t)
end

@inline function apply_edge_buffered!(::Type{T}, batch, i,
                                      du, u, s, gbuf, p, t) where {T}
    @inbounds begin
        _du  = _has_dynamic(T) ? view(du, dynamic_range(batch, i))   : nothing
        _u   = _has_dynamic(T) ? view(u,  dynamic_range(batch, i))   : nothing
        _s   = _has_static(T)  ? view(s,  static_range(batch, i))    : nothing
        _p   = _indexable(p)   ? view(p,  parameter_range(batch, i)) : p
        bufr = @views gbuf_range(batch, i)
        _src = @views gbuf[bufr, 1]
        _dst = @views gbuf[bufr, 2]
        apply_compf(T, compf(batch), _du, _u, _s, _src, _dst, _p, t)
    end
    nothing
end

@propagate_inbounds function apply_compf(::Type{<:ODEEdge}, f::F, du, u, s, src, dst, p, t) where {F}
    f(du, u, src, dst, p, t)
end
@propagate_inbounds function apply_compf(::Type{<:StaticEdge}, f::F, du, u, s, src, dst, p, t) where {F}
    f(s, src, dst, p, t)
end

_has_dynamic(T) = _has_dynamic(statetype(T))
_has_static(T) = _has_static(statetype(T))
_has_dynamic(::Dynamic) = true
_has_dynamic(::Static) = false
_has_static(::Dynamic) = false
_has_static(::Static) = true
# check if indexing into p is necessary
_indexable(::Nothing) = false
_indexable(::SciMLBase.NullParameters) = false
_indexable(::AbstractVector) = true
