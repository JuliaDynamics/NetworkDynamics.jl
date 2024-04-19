using Graphs

export dim, pdim

@enumx StateType dynamic static

mutable struct IndexManager
    v_data::OrderedDict{Int, Tuple{StateType.T, UnitRange{Int}}}
    e_data::OrderedDict{Int, Tuple{StateType.T, UnitRange{Int}}}
    v_para::OrderedDict{Int, UnitRange{Int}}
    e_para::OrderedDict{Int, UnitRange{Int}}
    size_p::Int
    size_dynamic::Int
    size_static::Int
    function IndexManager()
        new((OrderedDict{Int, Tuple{StateType.T,UnitRange{Int}}}() for i in 1:4)..., 0,0,0)
    end
end

# struct FlatGraphData{D,S,P}
#    dynamic::D
#    static::S
#    p::P
#    im::IndexManager
# end

# vertex_data_range(im::IndexManager, i) = im.v_data[i]
# vertex_para_range(im::IndexManager, i) = im.v_para[i]
# edge_data_range(im::IndexManager, i) = im.e_data[i]
# edge_para_range(im::IndexManager, i) = im.e_para[i]

abstract type ExecutionStyle{buffered} end
struct SequentialExecution{buffered} <: ExecutionStyle{buffered} end
struct ThreadedExecution{buffered} <: ExecutionStyle{buffered} end
usebuffer(::ExecutionStyle{buffered}) where {buffered} = buffered

struct Network{EX<:ExecutionStyle,NL,VTup}
    "vertex batches of same function"
    vertexbatches::VTup
    "network layer"
    nl::NL
    "index manager"
    im::IndexManager
    "lazy cache pool"
    cachepool::LazyBufferCache
end
executionstyle(::Network{ex}) where {ex} = ex

dim(nw::Network) = full_data_range(nw.im)[end]
pdim(nw::Network) = full_para_range(nw.im)[end]

struct NetworkLayer{GT,ETup,AF,MT}
    "graph/toplogy of layer"
    g::GT
    "edge batches with same function"
    edgebatches::ETup
    "accumulation function"
    accumulator::AF
    "depth of edge accumulation"
    edepth::Int # potential becomes range for multilayer
    "vertex dimensions visible to edges"
    vdepth::Int # potential becomes range for multilayer
    "mapping e_idx -> [v_src_idx_in_fullflat; v_dst_idx_in_fullflat]"
    gather_map::MT # input_map[:, e_idx] = [v_src_idx, v_dst_idx]
    function NetworkLayer(g, eb, ac, ed, vdepth, execution)
        if usebuffer(execution)
            mt = zeros(Int, ne(g)*vdepth, 2)
        else
            mt = (Vector{UnitRange{Int}}(undef, ne(g)),
                  Vector{UnitRange{Int}}(undef, ne(g)))
        end
        new{typeof(g),typeof(eb),typeof(ac),typeof(mt)}(g, eb, ac, ed, vdepth, mt)
    end
end


abstract type ComponentBatch{F} end

struct VertexBatch{F} <: ComponentBatch{F}
    "vertex indices contained in batch"
    indices::Vector{Int}
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

struct EdgeBatch{F,SM} <: ComponentBatch{F}
    "edge indices (as in edge iterator) contained in batch"
    indices::Vector{Int}
    "edge function"
    fun::F
    "edge dimension"
    dim::Int
    "first index of first batch element in flat state vector"
    firstidx::Int
    "parameter dimensions"
    pdim::Int
    "first index of first batch element in flat parameter Vector"
    pfirstidx::Int
    "Mapping idx_in_fullflat -> index in agg vector"
    scatter_map::SM
    function EdgeBatch(i,f,args...)
        cp = coupling(f)
        sm = if cp ∈ (AntiSymmetric(), Symmetric)
            (zeros(Int,0), zeros(Int,0))
        else
            zeros(Int,0)
        end
        new{typeof(f), typeof(sm)}(i, f, args..., sm)
    end
end

@inline Base.length(cb::ComponentBatch) = Base.length(cb.indices)
@inline statetype(::ComponentBatch{F}) where {F} = statetype(F)
coupling(::EdgeBatch{F}) where {F} = coupling(F)

@inline function state_range(batch::ComponentBatch, i)
    start = batch.firstidx + (i-1) * batch.dim
    start : start + batch.dim - 1
end

@inline function parameter_range(batch::ComponentBatch, i)
    start = batch.pfirstidx + (i-1) * batch.pdim
    start : start + batch.pdim - 1
end

function register_components!(im::IndexManager, idxs, compf)
    type = statetype(compf)
    ddict, pdict = if compf isa VertexFunction
        im.v_data, im.v_para
    elseif compf isa EdgeFunction
        im.e_data, im.e_para
    end

    d_first = _getsize(im, type)
    p_first = im.size_p
    for i in idxs
        d_size = _getsize(im, type)
        d_end = d_size + dim(compf)
        ddict[i] = (type, (d_size+1):d_end)
        _setsize!(im, type, d_end)

        p_size = im.size_p
        p_end = p_size + pdim(compf)
        pdict[i] = (p_size+1):p_end
        im.size_p = p_end
    end
    d_first, p_first
end
function _getsize(im::IndexManager, type::StateType.T)
    type == StateType.dynamic && return im.size_dynamic
    type == StateType.static  && return im.size_static
    error()
end
function _setsize!(im::IndexManager, type::StateType.T, s)
    if type == StateType.dynamic
        im.size_dynamic = s
    elseif type == StateType.static
        im.size_static = s
    else
        error()
    end
end

function isdense(im::IndexManager)
    pidxs = Int[]
    dynamic_idxs = Int[]
    static_idxs = Int[]
    for d in (im.v_data, im.e_data)
        for (k, v) in d
            if v[1] == StateType.dynamic
                append!(dynamic_idxs, v[2])
            elseif v[1] == StateType.static
                append!(static_idxs, v[2])
            else
                error()
            end
        end
    end
    for d in (im.v_para, im.e_para)
        for (k, v) in d
            append!(pidxs, v)
        end
    end
    sort!(pidxs)
    sort!(dynamic_idxs)
    sort!(static_idxs)
    @assert pidxs == 1:im.size_p
    @assert dynamic_idxs == 1:im.size_dynamic
    @assert static_idxs == 1:im.size_static
    return true
end
