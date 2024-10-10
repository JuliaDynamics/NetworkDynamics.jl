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
    for batch in batches
        for i in 1:length(batch)
            eidx = batch.indices[i]
            e = im.edgevec[eidx]
            dst_range = im.v_data[e.src][1:im.vdepth]
            src_range = im.v_data[e.dst][1:im.vdepth]
            range = gbuf_range(batch, i)
            map[range, 1] .= dst_range
            map[range, 2] .= src_range
        end
    end

    map2 = zeros(Int, ne(im.g) * im.vdepth, 2)
    for i in 1:ne(im.g)
        map2[im.e_gbufr[i], 1] .= im.e_src[i]
        map2[im.e_gbufr[i], 2] .= im.e_dst[i]
    end

    if map != map2
        error("nope not the same")
    end

    N = ForwardDiff.pickchunksize(max(im.lastidx_dynamic, im.lastidx_p))
    EagerGBufProvider(map, DiffCache(Float64.(map), N))
end

function get_gbuf(bufp::EagerGBufProvider, u)
    gbuf = get_tmp(bufp.diffcache, u)
    NNlib.gather!(gbuf, u, bufp.map)
    gbuf
end

@inline function get_src_dst(gbuf::AbstractArray, batch, i)
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

@inline function get_src_dst(gbuf::LazyGBuf, batch, i)
    eidx = @views batch.indices[i]
    src = @views gbuf.u[gbuf.lbp.e_src[eidx]]
    dst = @views gbuf.u[gbuf.lbp.e_dst[eidx]]
    src, dst
end
