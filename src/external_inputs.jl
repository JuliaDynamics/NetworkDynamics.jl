struct StateBufIdx
    idx::Int
end
struct OutBufIdx
    idx::Int
end
struct ExtMap{M<:AbstractVector{<:Union{StateBufIdx, OutBufIdx}}}
    map::M
end

function ExtMap(im::IndexManager)
    map = Vector{Union{StateBufIdx, OutBufIdx}}(undef, im.lastidx_extbuf)
    isempty(map) && error("There are no external inputs.")

    for (vi, vm) in pairs(im.vertexm)
        has_external_input(vm) || continue
        for (i, si) in pairs(vm.extin)
            map[im.v_ext[vi][i]] = _symidx_to_extidx(im, si)
        end
    end
    for (ei, em) in pairs(im.edgem)
        has_external_input(em) || continue
        for (i, si) in pairs(em.extin)
            map[im.e_ext[ei][i]] = _symidx_to_extidx(im, si)
        end
    end

    # narrow down type if all are output or state
    if all(e -> typeof(e) == typeof(first(map)), map)
        map = Vector{typeof(first(map))}(map)
    end
    ExtMap(map)
end

function _symidx_to_extidx(im, sni)
    cm = getcomp(im, sni)
    if subsym_has_idx(sni.subidx, sym(cm))
        range = getcomprange(im, sni)
        return StateBufIdx(range[subsym_to_idx(sni.subidx, sym(cm))])
    elseif subsym_has_idx(sni.subidx, outsym_flat(cm))
        range = getcompoutrange(im, sni)
        if hasff(cm)
            throw(ArgumentError("Cannot resolve external input $sni! Outputs of feed-forward components are not allowed as external inputs."))
        end
        return OutBufIdx(range[subsym_to_idx(sni.subidx, sym(cm))])
    else
        throw(ArgumentError("Cannot resolve external input $sni! External inputs musst be states or outputs of non-feed-forward components."))
    end
end

# CPU version
function collect_externals!(map::ExtMap, extbuf::Vector, u::Vector, o::Vector)
    # the performance of this function really depends on the fact
    # that julia somehow inferes that if src isn't a statebuf it musst be a outbuf
    # wasn't able to recreate something this fast in a more general way
    @inbounds for (dst, src) in pairs(map.map)
        if src isa StateBufIdx
            extbuf[dst] = u[src.idx]
        else
            extbuf[dst] = o[src.idx]
        end
    end
    extbuf
end

# GPU Version
function collect_externals!(map::ExtMap, extbuf, u, o)
    _backend = get_backend(extbuf)
    kernel = collect_ext_kernel!(_backend)
    kernel(map.map, extbuf, u, o; ndrange=length(extbuf))
end
@kernel function collect_ext_kernel!(map, extbuf, @Const(u), @Const(o))
    dst = @index(Global)
    src = map[dst]
    if src isa StateBufIdx
        extbuf[dst] = u[src.idx]
    else
        extbuf[dst] = o[src.idx]
    end
    nothing
end

has_external_input(c::ComponentModel) = !isnothing(c.extin)
has_external_input(cb::ComponentBatch) = !iszero(cb.extbufstride.strides)
has_external_input(nw::Network) = !isnothing(nw.extmap)
has_external_input(im::IndexManager) = !iszero(im.lastidx_extbuf)
