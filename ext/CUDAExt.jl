module CUDAExt
using NetworkDynamics: Network, NetworkLayer, VertexBatch, EdgeBatch,
                       KAAggregator, AggregationMap, SparseAggregator,
                       LazyGBufProvider, EagerGBufProvider,
                       dispatchT, compf, iscudacompatible, executionstyle
using NetworkDynamics.PreallocationTools: DiffCache
using NetworkDynamics: KernelAbstractions as KA

using CUDA: CuArray
using Adapt: Adapt, adapt

# main entry for bringing Network to GPU
function Adapt.adapt_structure(to, n::Network)
    if to isa KA.GPU
        throw(ArgumentError("Looks like to passed an KernelAbstractions backend to adapt Network to GPU. \
            this is not supported as the internal cache types cannot be infered without known the eltype. \
            Please adapt using `CuArray{Float32}` or `CuArray{Float64}`!"))
    end
    if !(to isa Type{<:CuArray})
        throw(ArgumentError("Can't handle Adaptor $to. \
            Please adapt using `CuArray{Float32}` or `CuArray{Float64}`!"))
    end
    if eltype(to) ∉ (Float32, Float64)
        throw(ArgumentError("Use adapt on Network with either `CuArray{Float32}` or `CuArray{Float64}` \
            such that internal caches can be created with the correct type!"))
    end
    if !iscudacompatible(n)
        throw(ArgumentError("The provided network has non-cuda compatible aggregator or exectuion strategies."))
    end
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

# overload to retain int types for aggregation map
Adapt.@adapt_structure AggregationMap
function Adapt.adapt_structure(to::Type{<:CuArray{<:AbstractFloat}}, am::AggregationMap)
    map = adapt(CuArray, am.map)
    symmap = adapt(CuArray, am.symmap)
    AggregationMap(am.range, map, am.symrange, symmap)
end

Adapt.@adapt_structure SparseAggregator
Adapt.@adapt_structure LazyGBufProvider

# overload to retain int types for eager gbuf map
function Adapt.adapt_structure(to::Type{<:CuArray{<:AbstractFloat}}, gbp::EagerGBufProvider)
    _adapt_eager_gbufp(CuArray, to, gbp)
end
function Adapt.adapt_structure(to, gbp::EagerGBufProvider)
    _adapt_eager_gbufp(to, to, gbp)
end
function _adapt_eager_gbufp(mapto, cacheto, gbp)
    map = adapt(mapto, gbp.map)
    cache = adapt_diffcache(cacheto, gbp.diffcache)
    EagerGBufProvider(map, cache)
end

# define similar to adapt_structure for DiffCache without type piracy
function adapt_diffcache(to, c::DiffCache)
    du = adapt(to, c.du)
    dual_du = adapt(to, c.dual_du)
    DiffCache(du, dual_du, c.any_du)
    # N = length(c.dual_du) ÷ length(c.du) - 1
    # DiffCache(du, N)
end

# keep Int types in indices for VertexBatch/EdgeBatch
function Adapt.adapt_structure(to::Type{<:CuArray{<:AbstractFloat}}, b::VertexBatch)
    Adapt.adapt_structure(CuArray, b)
end
function Adapt.adapt_structure(to::Type{<:CuArray{<:AbstractFloat}}, b::EdgeBatch)
    Adapt.adapt_structure(CuArray, b)
end
function Adapt.adapt_structure(to, b::VertexBatch)
    idxs = adapt(to, b.indices)
    VertexBatch{dispatchT(b), typeof(compf(b)), typeof(idxs)}(
        idxs, compf(b), b.statestride, b.pstride, b.aggbufstride)
end
function Adapt.adapt_structure(to, b::EdgeBatch)
    idxs = adapt(to, b.indices)
    EdgeBatch{dispatchT(b), typeof(compf(b)), typeof(idxs)}(
        idxs, compf(b), b.statestride, b.pstride, b.gbufstride)
end

end
