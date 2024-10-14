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
            _u = get_state_cache(nw, du)
            _u[1:nw.im.lastidx_dynamic] .= u
        end

        # We thought about splitting u and s instead of stacking
        # but this leads to problems with unbuffered execution (additional lookup necessary)
        # and possible gathering for buffered execution
        dupt = (du, _u, p, t)

        # NOTE: first all static vertices, than all edges (regarless) then dyn vertices
        # maybe disallow static vertices entierly, otherwise the order gets complicated
        gbuf =  get_gbuf(nw.gbufprovider, _u)
        @timeit_debug "process layer" process_layer!(ex, nw, nw.layer, gbuf, dupt)

        @timeit_debug "aggregate" begin
            aggbuf = get_aggregation_cache(nw, _u)
            aggregate!(nw.layer.aggregator, aggbuf, _u)
        end

        @timeit_debug "process vertices" process_vertices!(ex, nw, aggbuf, dupt)
    end
    return nothing
end

get_ustacked_buf(s) = get_ustacked_buf(s.nw, uflat(s), pflat(s), s.t)
function get_ustacked_buf(nw, u, p, t)
    _u = get_state_cache(nw, u)
    _u[1:nw.im.lastidx_dynamic] .= u
    dupt = (nothing, _u, p, t)
    ex = executionstyle(nw)
    gbuf = get_gbuf(nw.gbufprovider, _u)
    process_layer!(ex, nw, nw.layer, gbuf, dupt; filt=isstatic)
    aggbuf = get_aggregation_cache(nw, _u)
    aggregate!(nw.layer.aggregator, aggbuf, _u)
    process_vertices!(ex, nw, aggbuf, dupt; filt=isstatic)
    _u, aggbuf
end

####
#### Vertex Execution
####
@inline function process_vertices!(::SequentialExecution, nw, aggbuf, dupt; filt=nofilt)
    unrolled_foreach(filt, nw.vertexbatches) do batch
        (du, u, p, t) = dupt
        for i in 1:length(batch)
            _type = dispatchT(batch)
            apply_vertex!(_type, batch, i, du, u, aggbuf, p, t)
        end
    end
end

@inline function process_vertices!(::ThreadedExecution, nw, aggbuf, dupt; filt=nofilt)
    unrolled_foreach(filt, nw.vertexbatches) do batch
        (du, u, p, t) = dupt
        Threads.@threads for i in 1:length(batch)
            _type = dispatchT(batch)
            apply_vertex!(_type, batch, i, du, u, aggbuf, p, t)
        end
    end
end

@inline function process_vertices!(::PolyesterExecution, nw, aggbuf, dupt; filt=nofilt)
    unrolled_foreach(filt, nw.vertexbatches) do batch
        (du, u, p, t) = dupt
        Polyester.@batch for i in 1:length(batch)
            _type = dispatchT(batch)
            apply_vertex!(_type, batch, i, du, u, aggbuf, p, t)
        end
    end
end

@inline function process_vertices!(::KAExecution, nw, aggbuf, dupt; filt=nofilt)
    _backend = get_backend(dupt[2])
    unrolled_foreach(filt, nw.vertexbatches) do batch
        (du, u, p, t) = dupt
        kernel = vkernel!(_backend)
        kernel(dispatchT(batch), batch,
               du, u, aggbuf, p, t; ndrange=length(batch))
    end
    KernelAbstractions.synchronize(_backend)
end
@kernel function vkernel!(::Type{T}, @Const(batch),
                          du, @Const(u),
                          @Const(aggbuf), @Const(p), @Const(t)) where {T<:ODEVertex}
    I = @index(Global)
    apply_vertex!(T, batch, I, du, u, aggbuf, p, t)
    nothing
end
@kernel function vkernel!(::Type{T}, @Const(batch),
                          @Const(du), u,
                          @Const(aggbuf), @Const(p), @Const(t)) where {T<:StaticVertex}
    I = @index(Global)
    apply_vertex!(T, batch, I, du, u, aggbuf, p, t)
    nothing
end

@inline function apply_vertex!(::Type{T}, batch, i, du, u, aggbuf, p, t) where {T}
    @inbounds begin
        _du  = _has_dynamic(T) ? view(du, state_range(batch, i))   : nothing
        _u   = _has_dynamic(T) ? view(u,  state_range(batch, i))   : nothing
        _s   = _has_static(T)  ? view(u,  state_range(batch, i))    : nothing
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
#### Edge Layer Execution
####
@inline function process_layer!(::SequentialExecution, nw, layer, gbuf, dupt; filt=nofilt)
    unrolled_foreach(filt, layer.edgebatches) do batch
        (du, u, p, t) = dupt
        for i in 1:length(batch)
            _type = dispatchT(batch)
            apply_edge!(_type, batch, i, du, u, gbuf, p, t)
        end
    end
end

@inline function process_layer!(::ThreadedExecution, nw, layer, gbuf, dupt; filt=nofilt)
    unrolled_foreach(filt, layer.edgebatches) do batch
        (du, u, p, t) = dupt
        Threads.@threads for i in 1:length(batch)
            _type = dispatchT(batch)
            apply_edge!(_type, batch, i, du, u, gbuf, p, t)
        end
    end
end

@inline function process_layer!(::PolyesterExecution, nw, layer, gbuf, dupt; filt=nofilt)
    unrolled_foreach(filt, layer.edgebatches) do batch
        (du, u, p, t) = dupt
        Polyester.@batch for i in 1:length(batch)
            _type = dispatchT(batch)
            apply_edge!(_type, batch, i, du, u, gbuf, p, t)
        end
    end
end

@inline function process_layer!(::KAExecution, nw, layer, gbuf, dupt; filt=nofilt)
    _backend = get_backend(dupt[2])
    unrolled_foreach(filt, layer.edgebatches) do batch
        (du, u, p, t) = dupt
        kernel = ekernel!(_backend)
        kernel(dispatchT(batch), batch,
               du, u, gbuf, p, t; ndrange=length(batch))
    end
    KernelAbstractions.synchronize(_backend)
end
@kernel function ekernel!(::Type{T}, @Const(batch),
                          du, @Const(u), @Const(gbuf),
                          @Const(p), @Const(t)) where {T<:ODEEdge}
    I = @index(Global)
    apply_edge!(T, batch, I, du, u, gbuf, p, t)
end
@kernel function ekernel!(::Type{T}, @Const(batch),
                          @Const(du), u, @Const(gbuf),
                          @Const(p), @Const(t)) where {T<:StaticEdge}
    I = @index(Global)
    apply_edge!(T, batch, I, du, u, gbuf, p, t)
end

@inline function apply_edge!(::Type{T}, batch, i, du, u, gbuf, p, t) where {T}
    @inbounds begin
        _du  = _has_dynamic(T) ? view(du, state_range(batch, i))     : nothing
        _u   = _has_dynamic(T) ? view(u,  state_range(batch, i))     : nothing
        _s   = _has_static(T)  ? view(u,  state_range(batch, i))     : nothing
        _p   = _indexable(p)   ? view(p,  parameter_range(batch, i)) : p
        _src, _dst = get_src_dst(gbuf, batch, i)
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
