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
    isempty(map) && return ExtMap(map)

    for vm in im.vertexm
        for (i, si) in pairs(vm.extsym)
            map[im.v_ext[i]] = _symidx_to_extidx(im, si)
        end
    end
    for em in im.edgem
        for (i, si) in pairs(em.extsym)
            map[im.e_ext[i]] = _symidx_to_extidx(im, si)
        end
    end

    # narrow down type if all are output or state
    if all(e -> typeof(e) == typeof(first(map)), map)
        map = Vector{typeof(first(map))}(map)
    end
    ExtMap(map)
end

function _symidex_to_extidx(im, sni)
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

function collect_externals!(map::ExtMap, extbuf, u, o)
    @inbounds for (dst, src) in pairs(map.map)
        if src isa StateBufIdx
            extbuf[dst] = u[src.idx]
        else
            extbuf[dst] = o[src.idx]
        end
    end
end

has_external_inputs(c::ComponentModel) = !iszero(extdim(c))
has_external_inputs(cb::ComponentBatch) = !iszero(stridesT(cb.extbufstride))
