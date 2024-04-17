using Graphs

export dim, pdim

struct IndexManager
    v_data::OrderedDict{Int, UnitRange{Int}}
    e_data::OrderedDict{Int, UnitRange{Int}}
    v_para::OrderedDict{Int, UnitRange{Int}}
    e_para::OrderedDict{Int, UnitRange{Int}}
    function IndexManager()
        dict() = OrderedDict{Int, UnitRange{Int}}()
        new((dict() for i in 1:4)...)
    end
end

const EXECUTION_STYLES = (:seq, :threaded)

struct Network{EX,NL,VTup}
    "vertex batches of same function"
    vertexbatches::VTup
    "network layer"
    nl::NL
    "index manager"
    im::IndexManager
end

dim(nw::Network) = full_data_range(nw.im)[end]
pdim(nw::Network) = full_para_range(nw.im)[end]

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

accumulator_cache(nl::NetworkLayer, T) = getcache(nl.cachepool, T, nl.accdim * nv(nl.g))

struct ColorBatch{ETup}
    "edge indices (as in LG edge iterator) contained in batch"
    edges::Vector{Int}
    "edge batches of same function"
    edgebatches::ETup
end

struct EdgeBatch{F}
    "edge indices (as in LG edge iterator) contained in batch"
    edges::Vector{Int}
    "edge function"
    fun::F
    "edge dimension"
    dim::Int
    "first index of first batch element in flat state vector (-1 for static edge)"
    firstidx::Int
    "parameter dimensions"
    pdim::Int
    "first index of first batch element in flat parameter Vector"
    pfirstidx::Int
    "src: idx in g, idx in data, dim ; dst: idx in g, idx in data, dim"
    vertex_indices::Vector{Int}
end

component_idxs(eb::EdgeBatch) = eb.edges

@inline function vertex_indices(batch::EdgeBatch, i)
    v = batch.vertex_indices
    idx = (i-1)*6
    @inbounds begin
        (v[idx + 1], v[idx + 2], v[idx + 3], v[idx + 4], v[idx + 5], v[idx + 6])
    end
end

@inline function src_dst_ranges(layer::NetworkLayer, batch::EdgeBatch, i)
    (src, src_dat, src_dim, dst, dst_dat, dst_dim) = vertex_indices(batch, i)
    # ranges of src and data idx in data array
    src_r = src_dat : src_dat + src_dim - 1
    dst_r = dst_dat : dst_dat + dst_dim - 1

    # range of src vertex in accumulator array
    src_acc_first = 1 + (src-1) * layer.accdim
    src_acc = src_acc_first : src_acc_first + layer.accdim - 1

    # range of dst vertex in accumulator array
    dst_acc_first = 1 + (dst-1) * layer.accdim
    dst_acc = dst_acc_first : dst_acc_first + layer.accdim - 1

    # range of edge in flat state array
    edge_first = batch.firstidx + (i-1) * batch.dim
    edge_r = edge_first : edge_first + batch.dim - 1

    return (src_r, dst_r, src_acc, dst_acc, edge_r)
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

component_idxs(vb::VertexBatch) = vb.vertices

@inline function vertex_ranges(layer::NetworkLayer, batch::VertexBatch, i)
    # range of the vertex in data array
    vstart = batch.firstidx + (i-1) * batch.dim
    v = vstart : vstart + batch.dim - 1

    # range in accumulator array
    v_idx = batch.vertices[i]
    astart = 1 + (v_idx-1) * layer.accdim
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

Graphs.nv(nw::Network) = sum(nv.(nw.vertexbatches))
Graphs.nv(vb::VertexBatch) = length(vb.vertices)

Graphs.ne(nw::Network) = size(ne.(nl))
Graphs.ne(nl::NetworkLayer) = ne(nl.g)

_lastinrange(d::OrderedDict{Int, UnitRange{Int}}) = isempty(d) ? 0 : d[d.keys[end]][end]

_lastindex_data(im::IndexManager) = max(_lastinrange(im.v_data), _lastinrange(im.e_data))
_lastindex_para(im::IndexManager) = max(_lastinrange(im.v_para), _lastinrange(im.e_para))

function _register_comp!(im::IndexManager, data, para, idxs, dim, pdim)
    firstidx = datapos = _lastindex_data(im) + 1
    pfirstidx = parapos = _lastindex_para(im) + 1
    for i in idxs
        if haskey(data, i) || haskey(para, i)
            error("Index $i allready in $data or $para")
        end
        if dim > 0
            data[i] = datapos : datapos + dim - 1
            datapos += dim
        end
        if pdim > 0
            para[i] = parapos : parapos + pdim - 1
            parapos += pdim
        end
    end
    return firstidx, pfirstidx
end

register_edges!(im::IndexManager, idxs, dim, pdim) = _register_comp!(im, im.e_data, im.e_para, idxs, dim, pdim)
register_vertices!(im::IndexManager, idxs, dim, pdim) = _register_comp!(im, im.v_data, im.v_para, idxs, dim, pdim)

vertex_data_range(im::IndexManager, i) = im.v_data[i]
vertex_para_range(im::IndexManager, i) = im.v_para[i]
edge_data_range(im::IndexManager, i) = im.e_data[i]
edge_para_range(im::IndexManager, i) = im.e_para[i]

function _isdense(last, dicts...)
    a = Vector{Int}(undef, last)
    empty!(a)
    for d in dicts
        for v in values(d)
            append!(a, v)
        end
    end
    sort!(a)
    return a == 1:last
end

isdense(im::IndexManager) = (_isdense(_lastindex_data(im), im.v_data, im.e_data) &&
    _isdense(_lastindex_para(im), im.v_para, im.e_para))

full_data_range(im::IndexManager) = 1:_lastindex_data(im)
full_para_range(im::IndexManager) = 1:_lastindex_para(im)
