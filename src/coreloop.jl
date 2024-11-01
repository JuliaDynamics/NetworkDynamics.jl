function (nw::Network{A,B,C,D,E})(du::dT, u::T, p, t) where {A,B,C,D,E,dT,T}
    # if !(eachindex(du) == eachindex(u) == 1:nw.im.lastidx_dynamic)
    #     throw(ArgumentError("du or u does not have expected size $(nw.im.lastidx_dynamic)"))
    # end
    # if nw.im.lastidx_p > 0 && _indexable(p) && !(eachindex(p) == 1:nw.im.lastidx_p)
    #     throw(ArgumentError("p does not has expecte size $(nw.im.lastidx_p)"))
    # end

    ex = executionstyle(nw)
    fill!(du, zero(eltype(du)))
    o = get_output_cache(nw, du)
    fill!(o, convert(eltype(o), NaN))


    duopt = (du, u, o, p, t)
    aggbuf = get_aggregation_cache(nw, du)
    gbuf = get_gbuf(nw.gbufprovider, o)

    # vg without ff
    process_batches!(ex, Val{:g}(), !hasff, nw.vertexbatches, nothing, duopt)
    # eg without ff
    process_batches!(ex, Val{:g}(), !hasff, nw.layer.edgebatches, nothing, duopt)

    # gather the vertex results for edges with ff
    gather!(nw.gbufprovider, gbuf, o) # 2.6ms

    # execute f for the edges without ff
    process_batches!(ex, Val{:f}(), !hasff, nw.layer.edgebatches, gbuf, duopt)
    # execute f&g for edges with ff
    process_batches!(ex, Val{:fg}(), hasff, nw.layer.edgebatches, gbuf, duopt)

    # aggegrate the results
    aggregate!(nw.layer.aggregator, aggbuf, o)

    # vf for vertices without ff
    process_batches!(ex, Val{:f}(), !hasff, nw.vertexbatches, aggbuf, duopt)
    return nothing
end

@inline function process_batches!(::SequentialExecution, fg, filt::F, batches, inbuf, duopt) where {F}
    unrolled_foreach(filt, batches) do batch
        (du, u, o, p, t) = duopt
        for i in 1:length(batch)
            _type = dispatchT(batch)
            apply_comp!(_type, fg, batch, i, du, u, o, inbuf, p, t)
        end
    end
end

@inline function apply_comp!(::Type{<:UnifiedVertex}, fg, batch, i, du, u, o, aggbuf, p, t)
    @inbounds begin
        _o   = _needs_out(fg, batch) ? view(o, out_range(batch, i))         : nothing
        _du  = _needs_du(fg, batch)  ? view(du, state_range(batch, i))      : nothing
        _u   = _needs_u(fg, batch)   ? view(u,  state_range(batch, i))      : nothing
        _agg = _needs_in(fg, batch)  ? view(aggbuf, aggbuf_range(batch, i)) : nothing
        _p   = _needs_p(fg, batch) && _indexable(p) ? view(p,  parameter_range(batch, i))  : p
        evalf(fg, batch) && apply_compf(compf(batch), _du, _u, (_agg,), _p, t)
        evalg(fg, batch) && apply_compg(fftype(batch), compg(batch), (_o,), _u, (_agg,), _p, t)
    end
    nothing
end

@inline function apply_comp!(::Type{<:UnifiedEdge}, fg, batch, i, du, u, o, gbuf, p, t)
    @inbounds begin
        _odst = _needs_out(fg, batch) ? view(o, out_range(batch, i, 1))     : nothing
        _osrc = _needs_out(fg, batch) ? view(o, out_range(batch, i, 2))     : nothing
        _du   = _needs_du(fg, batch)  ? view(du, state_range(batch, i))     : nothing
        _u    = _needs_u(fg, batch)   ? view(u,  state_range(batch, i))     : nothing
        _ins  = _needs_in(fg, batch) ? get_src_dst(gbuf, batch, i)          : nothing
        _p    = _needs_p(fg, batch) && _indexable(p) ? view(p,  parameter_range(batch, i)) : p
        evalf(fg, batch) && apply_compf(compf(batch), _du, _u, _ins, _p, t)
        evalg(fg, batch) && apply_compg(fftype(batch), compg(batch), (_osrc, _odst), _u, _ins, _p, t)
    end
    nothing
end

@propagate_inbounds function apply_compf(f::F, du, u, ins, p, t) where {F}
    f(du, u, ins..., p, t)
    nothing
end

@propagate_inbounds function apply_compg(::PureFeedForward, g::G, outs, u, ins, p, t) where {G}
    g(outs..., ins..., p, t)
    nothing
end
@propagate_inbounds function apply_compg(::FeedForward, g::G, outs, u, ins, p, t) where {G}
    g(outs..., u, ins..., p, t)
    nothing
end
@propagate_inbounds function apply_compg(::NoFeedForward, g::G, outs, u, ins, p, t) where {G}
    g(outs..., u, p, t)
    nothing
end
@propagate_inbounds function apply_compg(::PureStateMap, g::G, outs, u, ins, p, t) where {G}
    g(outs..., u)
    nothing
end

_needs_du(fg, batch) = evalf(fg, batch)
_needs_u(fg, batch) = evalf(fg, batch) || fftype(batch) != PureFeedForward()
_needs_out(fg, batch) = evalg(fg, batch)
_needs_in(fg, batch) = evalf(fg, batch) || hasff(batch)
_needs_p(fg, batch) = evalf(fg, batch) || fftype(batch) != PureStateMap()

evalf(::Val{:f}, batch) = !isnothing(compf(batch))
evalf(::Val{:g}, batch) = false
evalf(::Val{:fg}, batch) = !isnothing(compf(batch))
evalg(::Val{:f}, _) = false
evalg(::Val{:g}, _) = true
evalg(::Val{:fg}, _) = true


# get_ustacked_buf(s) = get_ustacked_buf(s.nw, uflat(s), pflat(s), s.t)
# function get_ustacked_buf(nw, u, p, t)
#     _u = get_state_cache(nw, u)
#     _u[1:nw.im.lastidx_dynamic] .= u
#     dupt = (nothing, _u, p, t)
#     ex = executionstyle(nw)
#     gbuf = get_gbuf(nw.gbufprovider, _u)
#     process_layer!(ex, nw, nw.layer, gbuf, dupt; filt=isstatic)
#     aggbuf = get_aggregation_cache(nw, _u)
#     aggregate!(nw.layer.aggregator, aggbuf, _u)
#     process_vertices!(ex, nw, aggbuf, dupt; filt=isstatic)
#     _u, aggbuf
# end

####
#### Vertex Execution
####
# @inline function process_vertices!(::SequentialExecution, vbatches, aggbuf, dupt; filt=nofilt)
#     unrolled_foreach(filt, vbatches) do batch
#         (du, u, p, t) = dupt
#         for i in 1:length(batch)
#             _type = dispatchT(batch)
#             apply_vertex!(_type, batch, i, du, u, aggbuf, p, t)
#         end
#     end
# end

# @inline function process_vertices!(::ThreadedExecution, nw, aggbuf, dupt; filt=nofilt)
#     unrolled_foreach(filt, nw.vertexbatches) do batch
#         (du, u, p, t) = dupt
#         Threads.@threads for i in 1:length(batch)
#             _type = dispatchT(batch)
#             apply_vertex!(_type, batch, i, du, u, aggbuf, p, t)
#         end
#     end
# end

# @inline function process_vertices!(::PolyesterExecution, nw, aggbuf, dupt; filt=nofilt)
#     unrolled_foreach(filt, nw.vertexbatches) do batch
#         (du, u, p, t) = dupt
#         Polyester.@batch for i in 1:length(batch)
#             _type = dispatchT(batch)
#             apply_vertex!(_type, batch, i, du, u, aggbuf, p, t)
#         end
#     end
# end

# @inline function process_vertices!(::KAExecution, nw, aggbuf, dupt; filt=nofilt)
#     _backend = get_backend(dupt[2])
#     unrolled_foreach(filt, nw.vertexbatches) do batch
#         (du, u, p, t) = dupt
#         kernel = vkernel!(_backend)
#         kernel(dispatchT(batch), batch,
#                du, u, aggbuf, p, t; ndrange=length(batch))
#     end
#     KernelAbstractions.synchronize(_backend)
# end
# @kernel function vkernel!(::Type{T}, @Const(batch),
#                           du, @Const(u),
#                           @Const(aggbuf), @Const(p), @Const(fg, batch)) where {T<:ODEVertex}
#     I = @index(Global)
#     apply_vertex!(T, batch, I, du, u, aggbuf, p, t)
#     nothing
# end
# @kernel function vkernel!(::Type{T}, @Const(batch),
#                           @Const(du), u,
#                           @Const(aggbuf), @Const(p), @Const(t)) where {T<:StaticVertex}
#     I = @index(Global)
#     apply_vertex!(T, batch, I, du, u, aggbuf, p, t)
#     nothing
# end

# @inline function apply_vertex!(::Type{T}, batch, i, du, u, aggbuf, p, t) where {T}
#     @inbounds begin
#         _du  = _has_dynamic(T) ? view(du, state_range(batch, i))   : nothing
#         _u   = _has_dynamic(T) ? view(u,  state_range(batch, i))   : nothing
#         _s   = _has_static(T)  ? view(u,  state_range(batch, i))    : nothing
#         _p   = _indexable(p)   ? view(p,  parameter_range(batch, i)) : p
#         _agg = view(aggbuf, aggbuf_range(batch, i))
#         apply_compf(T, compf(batch), _du, _u, _s, _agg, _p, t)
#     end
#     nothing
# end

# @propagate_inbounds function apply_compf(::Type{<:ODEVertex}, f::F, du, u, s, agg, p, t) where {F}
#     f(du, u, agg, p, t)
# end
# @propagate_inbounds function apply_compf(::Type{<:StaticVertex}, f::F, du, u, s, agg, p, t) where {F}
#     f(s, agg, p, t)
# end

# ####
# #### Edge Layer Execution
# ####
# @inline function process_layer!(::SequentialExecution, nw, layer, gbuf, dupt; filt=nofilt)
#     unrolled_foreach(filt, layer.edgebatches) do batch
#         (du, u, p, t) = dupt
#         for i in 1:length(batch)
#             _type = dispatchT(batch)
#             apply_edge!(_type, batch, i, du, u, gbuf, p, t)
#         end
#     end
# end

# @inline function process_layer!(::ThreadedExecution, nw, layer, gbuf, dupt; filt=nofilt)
#     unrolled_foreach(filt, layer.edgebatches) do batch
#         (du, u, p, t) = dupt
#         Threads.@threads for i in 1:length(batch)
#             _type = dispatchT(batch)
#             apply_edge!(_type, batch, i, du, u, gbuf, p, t)
#         end
#     end
# end

# @inline function process_layer!(::PolyesterExecution, nw, layer, gbuf, dupt; filt=nofilt)
#     unrolled_foreach(filt, layer.edgebatches) do batch
#         (du, u, p, t) = dupt
#         Polyester.@batch for i in 1:length(batch)
#             _type = dispatchT(batch)
#             apply_edge!(_type, batch, i, du, u, gbuf, p, t)
#         end
#     end
# end

# @inline function process_layer!(::KAExecution, nw, layer, gbuf, dupt; filt=nofilt)
#     _backend = get_backend(dupt[2])
#     unrolled_foreach(filt, layer.edgebatches) do batch
#         (du, u, p, t) = dupt
#         kernel = ekernel!(_backend)
#         kernel(dispatchT(batch), batch,
#                du, u, gbuf, p, t; ndrange=length(batch))
#     end
#     KernelAbstractions.synchronize(_backend)
# end
# @kernel function ekernel!(::Type{T}, @Const(batch),
#                           du, @Const(u), @Const(gbuf),
#                           @Const(p), @Const(t)) where {T<:ODEEdge}
#     I = @index(Global)
#     apply_edge!(T, batch, I, du, u, gbuf, p, t)
# end
# @kernel function ekernel!(::Type{T}, @Const(batch),
#                           @Const(du), u, @Const(gbuf),
#                           @Const(p), @Const(t)) where {T<:StaticEdge}
#     I = @index(Global)
#     apply_edge!(T, batch, I, du, u, gbuf, p, t)
# end

# @inline function apply_edge!(::Type{T}, batch, i, du, u, gbuf, p, t) where {T}
#     @inbounds begin
#         _du  = _has_dynamic(T) ? view(du, state_range(batch, i))     : nothing
#         _u   = _has_dynamic(T) ? view(u,  state_range(batch, i))     : nothing
#         _s   = _has_static(T)  ? view(u,  state_range(batch, i))     : nothing
#         _p   = _indexable(p)   ? view(p,  parameter_range(batch, i)) : p
#         _src, _dst = get_src_dst(gbuf, batch, i)
#         apply_compf(T, compf(batch), _du, _u, _s, _src, _dst, _p, t)
#     end
#     nothing
# end

# @propagate_inbounds function apply_compf(::Type{<:ODEEdge}, f::F, du, u, s, src, dst, p, t) where {F}
#     f(du, u, src, dst, p, t)
# end
# @propagate_inbounds function apply_compf(::Type{<:StaticEdge}, f::F, du, u, s, src, dst, p, t) where {F}
#     f(s, src, dst, p, t)
# end

# _has_dynamic(T) = _has_dynamic(statetype(T))
# _has_static(T) = _has_static(statetype(T))
# _has_dynamic(::Dynamic) = true
# _has_dynamic(::Static) = false
# _has_static(::Dynamic) = false
# _has_static(::Static) = true
# check if indexing into p is necessary
_indexable(::Nothing) = false
_indexable(::SciMLBase.NullParameters) = false
_indexable(::AbstractVector) = true
