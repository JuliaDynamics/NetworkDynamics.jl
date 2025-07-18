module NetworkDynamicsCUDAExt
using NetworkDynamics: Network, NetworkLayer, ComponentBatch,
                       KAAggregator, AggregationMap, SparseAggregator,
                       LazyGBufProvider, EagerGBufProvider, LazyGBuf,
                       dispatchT, iscudacompatible, executionstyle, ExtMap
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
    caches = (;output = _adapt_diffcache(to, n.caches.output),
              aggregation = _adapt_diffcache(to, n.caches.aggregation),
              external = _adapt_diffcache(to, n.caches.external))
    exT = typeof(executionstyle(n))
    gT = typeof(n.im.g)
    extmap = adapt(to, n.extmap)

    Network{exT,gT,typeof(layer),typeof(vb),typeof(mm),eltype(caches),typeof(gbp),typeof(extmap)}(
        vb, layer, n.im, caches, mm, gbp, extmap, getfield(n, :jac_prototype))
end

Adapt.@adapt_structure NetworkLayer


####
#### Adapt Aggregators
####
Adapt.@adapt_structure KAAggregator

# overload to retain int types for aggregation map
Adapt.@adapt_structure AggregationMap
function Adapt.adapt_structure(to::Type{<:CuArray{<:AbstractFloat}}, am::AggregationMap)
    map = adapt(CuArray, am.map)
    AggregationMap(am.range, map)
end

Adapt.@adapt_structure SparseAggregator


####
#### Adapt GBufProviders
####
Adapt.@adapt_structure LazyGBufProvider
function Adapt.adapt_structure(to::Type{<:CuArray{<:AbstractFloat}}, gbp::LazyGBufProvider)
    adapt(CuArray, gbp) # preserve Vector{UnitRange}
end
Adapt.@adapt_structure LazyGBuf

# overload to retain int types for eager gbuf map
function Adapt.adapt_structure(to::Type{<:CuArray{<:AbstractFloat}}, gbp::EagerGBufProvider)
    _adapt_eager_gbufp(CuArray, to, gbp)
end
function Adapt.adapt_structure(to, gbp::EagerGBufProvider)
    _adapt_eager_gbufp(to, to, gbp)
end
function _adapt_eager_gbufp(mapto, cacheto, gbp)
    map = adapt(mapto, gbp.map)
    cache = _adapt_diffcache(cacheto, gbp.diffcache)
    EagerGBufProvider(map, cache)
end

####
#### Adapt external input map
####
Adapt.@adapt_structure ExtMap
function Adapt.adapt_structure(to::Type{<:CuArray{<:AbstractFloat}}, em::ExtMap)
    adapt(CuArray, em)
end


####
#### Adapt VertexBatch/EdgeBatch
####
function Adapt.adapt_structure(to::Type{<:CuArray{<:AbstractFloat}}, b::ComponentBatch)
    Adapt.adapt_structure(CuArray, b)
end
function Adapt.adapt_structure(to, b::ComponentBatch)
    indices = adapt(to, b.indices)
    ComponentBatch(dispatchT(b), indices, b.compf, b.compg, b.ff,
        b.statestride, b.pstride, b.inbufstride, b.outbufstride, b.extbufstride)
end

####
#### utils
####
# define similar to adapt_structure for DiffCache without type piracy
function _adapt_diffcache(to, c::DiffCache)
    du = adapt(to, c.du)
    dual_du = adapt(to, c.dual_du)
    DiffCache(du, dual_du, c.any_du)
    # N = length(c.dual_du) ÷ length(c.du) - 1
    # DiffCache(du, N)
end

end
