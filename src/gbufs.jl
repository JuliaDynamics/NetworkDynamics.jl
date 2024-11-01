abstract type GBufProvider end

####
#### Eager version
####

struct EagerGBufProvider{MT,C} <: GBufProvider
    "mapping e_idx -> [v_src_idx_in_fullflat; v_dst_idx_in_fullflat]"
    map::MT # input_map[:, e_idx] = [v_src_idx, v_dst_idx]
    diffcache::C
end

function EagerGBufProvider(im::IndexManager, batches)
    map = zeros(Int, ne(im.g) * im.vdepth, 2)
    for (i, e) in pairs(im.edgevec)
        map[im.e_gbufr[i], 1] .= im.v_out[e.src]
        map[im.e_gbufr[i], 2] .= im.v_out[e.dst]
    end

    N = ForwardDiff.pickchunksize(max(im.lastidx_dynamic, im.lastidx_p))
    EagerGBufProvider(map, DiffCache(Float64.(map), N))
end

get_gbuf(bufp::EagerGBufProvider, o) = get_tmp(bufp.diffcache, o)
gather!(bufp::EagerGBufProvider, gbuf, o) = NNlib.gather!(gbuf, o, bufp.map)

Base.@propagate_inbounds function get_src_dst(gbuf::AbstractArray, batch, i)
    bufr = @views gbuf_range(batch, i)
    src = @views gbuf[bufr, 1]
    dst = @views gbuf[bufr, 2]
    src, dst
end

####
#### Lazy version
####

struct LazyGBufProvider{SM,DM} <: GBufProvider
    e_src::SM
    e_dst::DM
end
struct LazyGBuf{LBP,UT}
    lbp::LBP
    u::UT
end

function LazyGBufProvider(im::IndexManager, _)
    e_src = Vector{UnitRange{Int}}(undef, ne(im.g))
    e_dst = Vector{UnitRange{Int}}(undef, ne(im.g))
    for (i, e) in pairs(im.edgevec)
        e_src[i] = im.v_out[e.src]
        e_dst[i] = im.v_out[e.dst]
    end
    LazyGBufProvider(e_src, e_dst)
end

get_gbuf(bufp::LazyGBufProvider, o) = LazyGBuf(bufp, o)
gather!(bufp::LazyGBufProvider, gbuf, o) = gbuf


Base.@propagate_inbounds function get_src_dst(gbuf::LazyGBuf, batch, i)
    eidx =  batch.indices[i]
    lbp = gbuf.lbp
    src = view(gbuf.u, lbp.e_src[eidx])
    dst = view(gbuf.u, lbp.e_dst[eidx])
    src, dst
end
