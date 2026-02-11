"""
    abstract type Aggregator end

Abstract sypertype for aggregators. Aggregators operate on the output buffer of
all components and fill the aggregation buffer with the aggregatated edge values
per vertex.

All aggregators have the constructor

    Aggegator(aggfun)

for example

    SequentialAggreator(+)
"""
abstract type Aggregator end

struct NaiveAggregator{F,ETup,G} <: Aggregator
    im::IndexManager{G}
    batches::ETup
    f::F
    function NaiveAggregator(im, batches, f)
        tup = Tuple(batches)
        new{typeof(f),typeof(tup),typeof(im.g)}(im, tup, f)
    end
end
NaiveAggregator(f) = (im, batches) -> NaiveAggregator(im, batches, f)

function aggregate!(a::NaiveAggregator, aggbuf, data)
    _aggregate!(a, a.batches, aggbuf, data)
end
function _aggregate!(a::NaiveAggregator, batches, aggbuf, data)
    unrolled_foreach(batches) do batch
        im = a.im
        for eidx in batch.indices
            edge = im.edgevec[eidx]

            # dst mapping
            source = @views data[im.e_out[eidx].dst]
            target = @views aggbuf[im.v_aggr[edge.dst]]
            target .= a.f.(target, source)

            # src mapping
            source = @views data[im.e_out[eidx].src]
            if !isempty(source)
                target = @views aggbuf[im.v_aggr[edge.src]]
                target .= a.f.(target, source)
            end
        end
    end
end


struct AggregationMap{V}
    range::UnitRange{Int} # range in data where aggregation is necessary
    map::V                # maps data idx to destination, if zero skip
end
function AggregationMap(im, batches)
    _map    = zeros(Int, im.lastidx_out)
    for batch in batches
        for eidx in batch.indices
            edge = im.edgevec[eidx]

            # dst mapping
            source = im.e_out[eidx].dst
            target = im.v_aggr[edge.dst]
            _map[source] .= target

            # src mapping
            source = im.e_out[eidx].src
            if !isempty(source)
                target = im.v_aggr[edge.src]
                _map[source] .= target
            end
        end
    end
    range, map       = _tighten_idxrange(_map)
    AggregationMap(range, map)
end
function _tighten_idxrange(v)
    first = findfirst(!iszero, v)
    isnothing(first) && return (1:0, similar(v, 0))
    last = findlast(!iszero, v)
    range = first:last
    return range, v[range]
end

"""
    KAAggregator(aggfun)

Parallel aggregation using [`KernelAbstractions.jl`](https://github.com/JuliaGPU/KernelAbstractions.jl).
Works with both GPU and CPU arrays.
"""
struct KAAggregator{F,V} <: Aggregator
    f::F
    m::AggregationMap{V}
end
KAAggregator(f) = (im, batches) -> KAAggregator(im, batches, f)
KAAggregator(im, batches, f) = KAAggregator(f, AggregationMap(im, batches))

function aggregate!(a::KAAggregator, aggbuf, data)
    am = a.m
    _backend = get_backend(data)
    # kernel = agg_kernel!(_backend, 1024, length(am.map))
    # kernel(a.f, aggbuf, view(data, am.range), am.map)
    kernel = agg_kernel!(_backend)
    kernel(a.f, aggbuf, view(data, am.range), am.map; ndrange=length(am.map))
    # TODO: synchronize after both aggregation sweeps?
    KernelAbstractions.synchronize(_backend)

    nothing
end

@kernel function agg_kernel!(f::F, aggbuf, @Const(data), @Const(idxs)) where {F}
    I = @index(Global)
    @inbounds if I ≤ length(idxs)
        _dst_i = idxs[I]
        if _dst_i != 0
            _dat = data[I]
            ref = Atomix.IndexableRef(aggbuf, (_dst_i,))
            Atomix.modify!(ref, f, _dat)
        end
    end
    nothing
end


"""
    SequentialAggregator(aggfun)

Sequential aggregation.
"""
struct SequentialAggregator{F} <: Aggregator
    f::F
    m::AggregationMap{Vector{Int}}
end
SequentialAggregator(f) = (im, batches) -> SequentialAggregator(im, batches, f)
SequentialAggregator(im, batches, f) = SequentialAggregator(f, AggregationMap(im, batches))

function aggregate!(a::SequentialAggregator, aggbuf, data)
    am = a.m
    @inbounds begin
        for (dat, dst_idx) in zip(view(data, am.range), am.map)
            if dst_idx != 0
                aggbuf[dst_idx] = a.f(aggbuf[dst_idx], dat)
            end
        end
    end

    nothing
end


"""
    PolyesterAggregator(aggfun)

Parallel aggregation using [`Polyester.jl`](https://github.com/JuliaSIMD/Polyester.jl).
"""
struct PolyesterAggregator{F} <: Aggregator
    f::F
    m::Vector{Tuple{Int64,Vector{Int64}}}
end
PolyesterAggregator(f) = (im, batches) -> PolyesterAggregator(im, batches, f)
PolyesterAggregator(im, batches, f) = PolyesterAggregator(f, _inv_aggregation_map(im, batches))

function aggregate!(a::PolyesterAggregator, aggbuf, data)
    length(a.m) == length(aggbuf) || throw(DimensionMismatch("length of aggbuf and a.m must be equal"))

    maxdepth = mapreduce(x -> length(x[2]), max, a.m)

    Polyester.@batch for (dstidx, srcidxs) in a.m
       for srcidx in srcidxs
           aggbuf[dstidx] = a.f(aggbuf[dstidx], data[srcidx])
       end
    end
    nothing
end


"""
    ThreadedAggregator(aggfun)

Parallel aggregation using Julia threads.
"""
struct ThreadedAggregator{F} <: Aggregator
    f::F
    m::Vector{Tuple{Int64,Vector{Int64}}}
end
ThreadedAggregator(f) = (im, batches) -> ThreadedAggregator(im, batches, f)
ThreadedAggregator(im, batches, f) = ThreadedAggregator(f, _inv_aggregation_map(im, batches))

function aggregate!(a::ThreadedAggregator, aggbuf, data)
    length(a.m) == length(aggbuf) || throw(DimensionMismatch("length of aggbuf and a.m must be equal"))

    Threads.@threads for (dstidx, srcidxs) in a.m
       @inbounds for srcidx in srcidxs
           aggbuf[dstidx] = a.f(aggbuf[dstidx], data[srcidx])
       end
    end
    nothing
end

function _inv_aggregation_map(im, batches)
    srcidxs = map(im.v_aggr, Graphs.degree(im.g)) do dst, ndegr
        [empty!(Vector{Int}(undef, ndegr)) for _ in 1:length(dst)]
    end
    unrolled_foreach(batches) do batch
        for eidx in batch.indices
            edge = im.edgevec[eidx]

            # dst mapping
            edat_idx = im.e_out[eidx].dst
            _pusheach!(srcidxs[edge.dst], edat_idx)

            # src mapping
            edat_idx = im.e_out[eidx].src
            if !isempty(edat_idx)
                _pusheach!(srcidxs[edge.src], edat_idx)
            end
        end
    end
    ret = Vector{Tuple{Int,Vector{Int}}}(undef, im.lastidx_aggr)
    empty!(ret)
    for (dstr, srcs) in zip(im.v_aggr, srcidxs)
        for (dst, src) in zip(dstr, srcs)
            push!(ret, (dst, src))
        end
    end
    ret
end
function _pusheach!(target, src)
    for i in eachindex(src)
        push!(target[i], src[i])
    end
end

"""
    SparseAggregator(+)

Only works with additive aggregation `+`. Aggregates via sparse inplace matrix
multiplication. Works with GPU Arrays.
"""
struct SparseAggregator{M} <: Aggregator
    m::M
    SparseAggregator(m::AbstractMatrix) = new{typeof(m)}(m)
end
function SparseAggregator(f)
    @argcheck f===(+) ArgumentError("Sparse Aggregator only works with + as reducer.")
    SparseAggregator
end
function SparseAggregator(im, batches)
    # sparse multiply is faster with Matrix{Float} . Vector{Float} than int!
    # (both on GPU and CPU)
    I, J, V = Float64[], Float64[], Float64[]
    unrolled_foreach(batches) do batch
        for eidx in batch.indices
            edge = im.edgevec[eidx]
            edgef = im.edgem[eidx]

            # dst mapping
            source = im.e_out[eidx].dst
            target = im.v_aggr[edge.dst]
            append!(I, target)
            append!(J, source)
            append!(V, Iterators.repeated(1, length(target)))

            # src mapping
            source = im.e_out[eidx].src
            if !isempty(source)
                target = im.v_aggr[edge.src]
                append!(I, target)
                append!(J, source)
                append!(V, Iterators.repeated(1, length(target)))
            end
        end
    end

    SparseAggregator(sparse(I,J,V, im.lastidx_aggr, im.lastidx_out))
end

function aggregate!(a::SparseAggregator, aggbuf, out)
   LinearAlgebra.mul!(aggbuf, a.m, out, 1, 1)
   nothing
end

# functions to retrieve the constructor of an aggregator for remake of network
get_aggr_constructor(a::NaiveAggregator) = NaiveAggregator(a.f)
get_aggr_constructor(a::KAAggregator) = KAAggregator(a.f)
get_aggr_constructor(a::SequentialAggregator) = SequentialAggregator(a.f)
get_aggr_constructor(a::PolyesterAggregator) = PolyesterAggregator(a.f)
get_aggr_constructor(a::ThreadedAggregator) = ThreadedAggregator(a.f)
get_aggr_constructor(a::SparseAggregator) = SparseAggregator(+)

iscudacompatible(::Type{<:KAAggregator}) = true
iscudacompatible(::Type{<:SparseAggregator}) = true

aggfun(agg::Aggregator) = agg.f
aggfun(agg::SparseAggregator) = +
