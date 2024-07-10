function Adapt.adapt_structure(to, n::Network)
    vb = adapt(to, n.vertexbatches)
    layer = adapt(to, n.layer)
    mm = adapt(to, n.mass_matrix)
    exT = typeof(executionstyle(n))
    gT = typeof(n.im.g)
    Network{exT,gT,typeof(layer),typeof(vb),typeof(mm)}(
        vb, layer, n.im, n.cachepool, mm)
end

Adapt.@adapt_structure NetworkLayer
Adapt.@adapt_structure KAAggregator
Adapt.@adapt_structure AggregationMap

function Adapt.adapt_structure(to, b::VertexBatch)
    idxs = adapt(to, b.indices)
    VertexBatch{compT(b), typeof(b.compf), typeof(idxs)}(
        idxs, b.compf, b.statestride, b.pstride, b.aggbufstride)
end
function Adapt.adapt_structure(to, b::EdgeBatch)
    idxs = adapt(to, b.indices)
    EdgeBatch{compT(b), typeof(b.compf), typeof(idxs)}(
        idxs, b.compf, b.statestride, b.pstride, b.gbufstride)
end
