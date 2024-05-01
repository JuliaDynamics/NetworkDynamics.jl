using Graphs

export dim, pdim

abstract type StateType end
struct Dynamic <: StateType end
struct Static <: StateType end

mutable struct IndexManager{G}
    g::G
    edgevec::Vector{SimpleEdge{Int64}}
    v_data::Vector{UnitRange{Int}}  # v data in flat states
    v_para::Vector{UnitRange{Int}}  # v para in flat para
    v_aggr::Vector{UnitRange{Int}}  # v input in aggbuf
    e_data::Vector{UnitRange{Int}}  # e data in flat states
    e_para::Vector{UnitRange{Int}}  # e para in flat para
    e_src::Vector{UnitRange{Int}}   # e input 1 in flat states
    e_dst::Vector{UnitRange{Int}}   # e input 2 in flat states
    e_gbufr::Vector{UnitRange{Int}} # e input range in gather buffer
    edepth::Int
    vdepth::Int
    lastidx_dynamic::Int
    lastidx_static::Int
    lastidx_p::Int
    lastidx_aggr::Int
    lastidx_gbuf::Int
    function IndexManager(g, dyn_states, edepth, vdepth)
        new{typeof(g)}(g, collect(edges(g)),
                       (Vector{UnitRange{Int}}(undef, nv(g)) for i in 1:3)...,
                       (Vector{UnitRange{Int}}(undef, ne(g)) for i in 1:5)...,
                       edepth, vdepth,
                       0, dyn_states, 0, 0, 0)
    end
end

abstract type ExecutionStyle{buffered} end
struct SequentialExecution{buffered} <: ExecutionStyle{buffered} end
struct ThreadedExecution{buffered} <: ExecutionStyle{buffered} end
usebuffer(::ExecutionStyle{buffered}) where {buffered} = buffered
usebuffer(::Type{<:ExecutionStyle{buffered}}) where {buffered} = buffered

struct Network{EX<:ExecutionStyle,G,NL,VTup}
    "vertex batches of same function"
    vertexbatches::VTup
    "network layer"
    layer::NL
    "index manager"
    im::IndexManager{G}
    "lazy cache pool"
    cachepool::LazyBufferCache
end
executionstyle(::Network{ex}) where {ex} = ex
@inline nvbatches(::Network) = length(vertexbatches)
dim(nw::Network) = nw.im.lastidx_dynamic
pdim(nw::Network) = nw.im.lastidx_p

struct NetworkLayer{GT,ETup,AF,MT}
    "graph/toplogy of layer"
    g::GT
    "edge batches with same function"
    edgebatches::ETup
    "aggregator object"
    aggregator::AF
    "depth of edge accumulation"
    edepth::Int # potential becomes range for multilayer
    "vertex dimensions visible to edges"
    vdepth::Int # potential becomes range for multilayer
    "mapping e_idx -> [v_src_idx_in_fullflat; v_dst_idx_in_fullflat]"
    gather_map::MT # input_map[:, e_idx] = [v_src_idx, v_dst_idx]
end
@inline nebatches(::NetworkLayer) = length(edgebatches)


abstract type ComponentBatch{F} end

struct VertexBatch{F} <: ComponentBatch{F}
    "vertex indices contained in batch"
    indices::Vector{Int}
    "vertex function"
    comp::F
    "state: dimension and first index"
    statestride::BatchStride
    "parameter: dimension and first index"
    pstride::BatchStride
    "aggregation: dimension and first index"
    aggbufstride::BatchStride
end

struct EdgeBatch{F} <: ComponentBatch{F}
    "edge indices (as in edge iterator) contained in batch"
    indices::Vector{Int}
    "edge function"
    comp::F
    "state: dimension and first index"
    statestride::BatchStride
    "parameter: dimension and first index"
    pstride::BatchStride
    "gathered vector: dimension and first index"
    gbufstride::BatchStride
end

@inline Base.length(cb::ComponentBatch) = Base.length(cb.indices)
@inline statetype(::ComponentBatch{F}) where {F} = statetype(F)
@inline coupling(::EdgeBatch{F}) where {F} = coupling(F)
@inline comptype(::ComponentBatch{F}) where {F} = F
@inline compf(b::ComponentBatch) = b.comp.f
@inline compf(b::NamedTuple) = b.f

@inline state_range(batch) = _fullrange(batch.statestride, length(batch))

@inline state_range(batch, i)     = _range(batch.statestride, i)
@inline parameter_range(batch, i) = _range(batch.pstride, i)
@inline aggbuf_range(batch, i)    = _range(batch.aggbufstride, i)
@inline gbuf_range(batch, i)      = _range(batch.gbufstride, i)

function register_vertices!(im::IndexManager, idxs, comp)
    for i in idxs
        im.v_data[i] = _nextdatarange!(im, statetype(comp), dim(comp))
        im.v_para[i] = _nextprange!(im, pdim(comp))
        im.v_aggr[i] = _nextaggrrange!(im, im.edepth)
    end
    (BatchStride(first(im.v_data[first(idxs)]), dim(comp)),
     BatchStride(first(im.v_para[first(idxs)]), pdim(comp)),
     BatchStride(first(im.v_aggr[first(idxs)]), im.edepth))
end
function register_edges!(im::IndexManager, idxs, comp)
    edgevec = collect(edges(im.g))
    for i in idxs
        e = edgevec[i]
        im.e_data[i] = _nextdatarange!(im, statetype(comp), dim(comp))
        im.e_para[i] = _nextprange!(im, pdim(comp))
        im.e_src[i] = im.v_data[e.src][1:im.vdepth]
        im.e_dst[i] = im.v_data[e.dst][1:im.vdepth]
        im.e_gbufr[i] = _nextgbufrange!(im, im.vdepth)
    end
    (BatchStride(first(im.e_data[first(idxs)]), dim(comp)),
     BatchStride(first(im.e_para[first(idxs)]), pdim(comp)),
     BatchStride(first(im.e_gbufr[first(idxs)]), im.vdepth))
end
function _nextdatarange!(im::IndexManager, ::Dynamic, N)
    newlast, range = _nextrange(im.lastidx_dynamic, N)
    im.lastidx_dynamic = newlast
    return range
end
function _nextdatarange!(im::IndexManager, ::Static, N)
    newlast, range = _nextrange(im.lastidx_static, N)
    im.lastidx_static = newlast
    return range
end
function _nextprange!(im::IndexManager, N)
    newlast, range = _nextrange(im.lastidx_p, N)
    im.lastidx_p = newlast
    return range
end
function _nextaggrrange!(im::IndexManager, N)
    newlast, range = _nextrange(im.lastidx_aggr, N)
    im.lastidx_aggr = newlast
    return range
end
function _nextgbufrange!(im::IndexManager, N)
    newlast, range = _nextrange(im.lastidx_gbuf, N)
    im.lastidx_gbuf = newlast
    return range
end
_nextrange(last, N) = last + N, last+1:last+N

function isdense(im::IndexManager)
    pidxs = Int[]
    stateidxs = Int[]
    for dataranges in (im.v_data, im.e_data)
        for range in dataranges
            append!(stateidxs, range)
        end
    end
    for pararanges in (im.v_para, im.e_para)
        for range in pararanges
            append!(pidxs, range)
        end
    end
    sort!(pidxs)
    sort!(stateidxs)
    @assert pidxs == 1:im.lastidx_p
    @assert stateidxs == 1:im.lastidx_static
    return true
end
