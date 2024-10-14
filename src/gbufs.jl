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
    for i in 1:ne(im.g)
        map[im.e_gbufr[i], 1] .= im.e_src[i]
        map[im.e_gbufr[i], 2] .= im.e_dst[i]
    end

    N = ForwardDiff.pickchunksize(max(im.lastidx_dynamic, im.lastidx_p))
    EagerGBufProvider(map, DiffCache(Float64.(map), N))
end

function get_gbuf(bufp::EagerGBufProvider, u)
    gbuf = get_tmp(bufp.diffcache, u)
    NNlib.gather!(gbuf, u, bufp.map)
    gbuf
end

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

LazyGBufProvider(im::IndexManager, _) = LazyGBufProvider(copy(im.e_src), copy(im.e_dst))

function get_gbuf(bufp::LazyGBufProvider, u)
    LazyGBuf(bufp, u)
end

Base.@propagate_inbounds function get_src_dst(gbuf::LazyGBuf, batch, i)
    eidx =  batch.indices[i]
    lbp = gbuf.lbp
    src = view(gbuf.u, lbp.e_src[eidx])
    dst = view(gbuf.u, lbp.e_dst[eidx])
    src, dst
end
