function Adapt.adapt_structure(to, n::Network)
    vb = adapt(to, n.vertexbatches)
    layer = adapt(to, n.layer)
    mm = adapt(to, n.mass_matrix)
    gbp = adapt(to, n.gbufprovider)
    caches = (;state = adapt_diffcache(to, n.caches.state),
              aggregation = adapt_diffcache(to, n.caches.aggregation))
    exT = typeof(executionstyle(n))
    gT = typeof(n.im.g)

    Network{exT,gT,typeof(layer),typeof(vb),typeof(mm),eltype(caches),typeof(gbp)}(
        vb, layer, n.im, caches, mm, gbp)
end

Adapt.@adapt_structure NetworkLayer
Adapt.@adapt_structure KAAggregator
Adapt.@adapt_structure AggregationMap
Adapt.@adapt_structure SparseAggregator
Adapt.@adapt_structure LazyGBufProvider

function Adapt.adapt_structure(to, gbp::EagerGBufProvider)
    map = adapt(to, gbp.map)
    cache = adapt_diffcache(to, gbp.diffcache)
    EagerGBufProvider(map, cache)
end

function adapt_diffcache(to, c::DiffCache)
    du = adapt(to, c.du)
    dual_du = adapt(to, c.dual_du)
    DiffCache(du, dual_du, c.any_du)
end

function Adapt.adapt_structure(to, b::VertexBatch)
    idxs = adapt(to, b.indices)
    VertexBatch{dispatchT(b), typeof(b.compf), typeof(idxs)}(
        idxs, b.compf, b.statestride, b.pstride, b.aggbufstride)
end
function Adapt.adapt_structure(to, b::EdgeBatch)
    idxs = adapt(to, b.indices)
    EdgeBatch{dispatchT(b), typeof(b.compf), typeof(idxs)}(
        idxs, b.compf, b.statestride, b.pstride, b.gbufstride)
end
