module CUDAExt
using NetworkDynamics: Network, NetworkLayer, VertexBatch, EdgeBatch,
                       KAAggregator, AggregationMap, SparseAggregator,
                       LazyGBufProvider, EagerGBufProvider,
                       dispatchT, compf, iscudacompatible, executionstyle
using NetworkDynamics.PreallocationTools: DiffCache
using NetworkDynamics: KernelAbstractions as KA

using CUDA: CuArray
using Adapt: Adapt, adapt

function Adapt.adapt_structure(to, n::Network)
    if to isa KA.GPU
        throw(ArgumentError("Looks like to passed an KernelAbstractions backend to adapt Network to GPU. \
            this is not supported as the internal cache types cannot be infered without known the eltype. \
            Please adapt using `CuArray{Float32}` or `CuArray{Float64}`!"))
    end
    if !(to isa Type{<:CuArray})
        throw(ArgumentError("Can't handle Adaptor $to.\
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
Adapt.@adapt_structure AggregationMap
Adapt.@adapt_structure SparseAggregator
Adapt.@adapt_structure LazyGBufProvider

function Adapt.adapt_structure(to::Type{<:CuArray}, gbp::EagerGBufProvider)
    map = adapt(CuArray{Int32}, gbp.map)
    cache = adapt_diffcache(to, gbp.diffcache)
    EagerGBufProvider(map, cache)
end

function adapt_diffcache(to::Type{<:CuArray}, c::DiffCache)
    du = adapt(to, c.du)
    dual_du = adapt(to, c.dual_du)
    DiffCache(du, dual_du, c.any_du)
end

function Adapt.adapt_structure(to::Type{<:CuArray}, b::VertexBatch)
    # idxs = adapt(CuArray{Int32}, b.indices)
    idxs = adapt(CuArray, b.indices)
    VertexBatch{dispatchT(b), typeof(compf(b)), typeof(idxs)}(
        idxs, compf(b), b.statestride, b.pstride, b.aggbufstride)
end
function Adapt.adapt_structure(to::Type{<:CuArray}, b::EdgeBatch)
    # idxs = adapt(CuArray{Int32}, b.indices)
    idxs = adapt(CuArray, b.indices)
    EdgeBatch{dispatchT(b), typeof(compf(b)), typeof(idxs)}(
        idxs, compf(b), b.statestride, b.pstride, b.gbufstride)
end

end
