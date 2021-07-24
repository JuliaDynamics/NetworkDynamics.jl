using LightGraphs

struct DynamicNetwork{NL,VTup}
    "vertex batches of same function"
    vertexbatches::VTup
    "network layer"
    nl::NL
end

struct NetworkLayer{GT,CTup,AF}
    "graph/toplogy of layer"
    g::GT
    "batches of edges with same color"
    colorbatches::CTup
    "accumulation function"
    accumulator::AF
    "depth/dim of accumulation"
    accdim::Int
    "lazy cache pool"
    cachepool::CachePool
end

accumulator_cache(nl::NetworkLayer, T) = getcache(nl.network.cachepool, T, nl.accdim, nv(nl))

struct ColorBatch{ETup}
    "edge indices (as in LG edge iterator) contained in batch"
    edges::Vector{Int}
    "edge batches of same function"
    edgebatches::ETup
end

abstract type EdgeBatch end

struct StaticEdgeBatch{F} <: EdgeBatch
    "edge indices (as in LG edge iterator) contained in batch"
    edges::Vector{Int}
    "edge function"
    fun::F
    "edge dimension"
    dim::Int
    "parameter dimensions"
    pdim::Int
    "first index of first batch element in flat parameter Vector"
    pfirstidx::Int
    "first index of first batch element in acc cache"
    accfirstidx::Int
    "src idx in data, src dim, dst idx in data, dst dim, dst index in graph"
    vertex_indices::Vector{Int}
end

struct ODEEdgeBatch{F} <: EdgeBatch
    "edge indices (as in LG edge iterator) contained in batch"
    edges::Vector{Int}
    "edge function"
    fun::F
    "edge dimension"
    dim::Int
    "first index of first batch element in flat state Vector"
    firstidx::Int
    "parameter dimensions"
    pdim::Int
    "first index of first batch element in flat parameter Vector"
    pfirstidx::Int
    "first index of first batch element in acc cache"
    accfirstidx::Int
    "src idx in data, src dim, dst idx in data, dst dim, dst index in graph"
    vertex_indices::Vector{Int}
end

@inline function vertex_indices(batch::EdgeBatch, i)
    v = batch.vertex_indices
    idx = (i-1) * 5
    @inbounds begin
        (v[idx + 1], v[idx + 2], v[idx + 3], v[idx + 4], v[idx + 5])
    end
end

@inline function vertex_ranges(layer::NetworkLayer, batch::StaticEdgeBatch, i)
    (src_dat, src_dim, dst_dat, dst_dim, dst_idx) = vertex_indices(batch, i)
    # ranges of src and data idx in data array
    src = src_dat : src_dat + src_dim - 1
    dst = src_dat : src_dat + src_dim - 1

    # range of dst vertex in accumolator array
    astart = batch.accfirstidx + (dst_idx-1) * layer.accdim
    acc = astart : astart + layer.accdim - 1

    return (src, dst, acc)
end

struct VertexBatch{F}
    "vertex indices contained in batch"
    vertices::Vector{Int}
    "vertex function"
    fun::F
    "vertex dimension"
    dim::Int
    "first index of first batch element in flat state Vector"
    firstidx::Int
    "parameter dimensions"
    pdim::Int
    "first index of first batch element in flat parameter Vector"
    pfirstidx::Int
end

@inline function vertex_ranges(layer::NetworkLayer, batch::VertexBatch, i)
    # range of the vertex in data array
    vstart = batch.firstidx + (i-1) * batch
    v = vstart : vstart + batch.dim - 1

    # range in accumulator array
    v_idx = batch.vertices[i]
    astart = (v_idx-1) * layer.accdim
    acc = astart : astart + layer.accdim - 1

    return (v, acc)
end

@inline function parameter_range(batch::Union{EdgeBatch,VertexBatch}, i)
    pstart = batch.pfirstidx + (i-1) * batch.pdim
    p = pstart : pstart + batch.pdim - 1
end

Base.length(cb::ColorBatch) = length(cb.edges)
Base.length(eb::EdgeBatch) = length(eb.edges)
Base.length(vb::VertexBatch) = length(vb.vertices)

LightGraphs.nv(nw::DynamicNetwork) = sum(nv.(nw.vertexbatches))
LightGraphs.nv(vb::VertexBatch) = length(vb.vertices)

LightGraphs.ne(nw::DynamicNetwork) = size(ne.(nl))
LightGraphs.ne(nl::NetworkLayer) = ne(nl.g)
