function Adapt.adapt_structure(to, l::NetworkLayer)
    # ebatches = (adapt(to, eb) for eb in l.edgebatches)
    aggr = adapt(to, l.aggregator)
    gmap = adapt(to, l.gather_map)
    NetworkLayer(l.g, l.edgebatches, aggr, l.edepth, l.vdepth, gmap)
end

function Adapt.adapt_structure(to, n::Network)
    layer = adapt(to, n.layer)
    Network{executionstyle(n),typeof(n.im.g),typeof(layer),typeof(n.vertexbatches)}(n.vertexbatches, layer, n.im, n.cachepool)
end

function Adapt.adapt_structure(to, a::NNlibScatter)
    dstmaps = [adapt(to, map) for map in a.dstmaps]
    srcmaps = [adapt(to, map) for map in a.srcmaps]
    NNlibScatter(a.f, a.batchranges, a.couplings, dstmaps, srcmaps, a.aggrsize)
end

function Adapt.adapt_structure(to, a::KAAggregator)
    KAAggregator(a.f, adapt(to, a.m))
end
function Adapt.adapt_structure(to, m::AggregationMap)
    AggregationMap(m.range, adapt(to, m.map),
                   m.symrange, adapt(to, m.symmap))
end


# XXX: get rid of essence() hack in favor for adapt, since Metal is not an important backand for now
function essence(b::VertexBatch)
    (;f=compf(b),
     statestride=b.statestride,
     pstride=b.pstride,
     aggbufstride=b.aggbufstride)
end
function essence(b::EdgeBatch)
    (;f=compf(b),
     statestride=b.statestride,
     pstride=b.pstride,
     gbufstride=b.gbufstride,
     indices=b.indices)
end
