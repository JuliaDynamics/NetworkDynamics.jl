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
    fill!(aggbuf, zero(eltype(aggbuf)))
    _aggregate!(a, a.batches, aggbuf, data)
end
function _aggregate!(a::NaiveAggregator, batches, aggbuf, data)
    unrolled_foreach(batches) do batch
        im = a.im
        for eidx in batch.indices
            edge = im.edgevec[eidx]
            # dst mapping
            target = @views aggbuf[im.v_aggr[edge.dst]]
            source = @views data[im.e_data[eidx][1:im.edepth]]
            target .= a.f.(target, source)

            # src mapping
            cplng = coupling(batch)
            if cplng == Symmetric()
                target = @views aggbuf[im.v_aggr[edge.src]]
                source = @views data[im.e_data[eidx][1:im.edepth]]
                target .= a.f.(target, source)
            elseif cplng == AntiSymmetric()
                target = @views aggbuf[im.v_aggr[edge.src]]
                source = @views data[im.e_data[eidx][1:im.edepth]]
                target .= a.f.(target, -1 .* source)
            elseif cplng == Fiducial()
                target = @views aggbuf[im.v_aggr[edge.src]]
                source = @views data[im.e_data[eidx][im.edepth+1:2*im.edepth]]
                target .= a.f.(target, source)
            end
        end
    end
end


struct AggregationMap{V}
    range::UnitRange{Int} # range in data where aggregation is necessary
    map::V                # maps data idx to destination, if zero skip
    symrange::UnitRange{Int}     # same for symmetric/antisymmetric coupling
    symmap::V
end
function AggregationMap(im, batches)
    _map    = zeros(Int, im.lastidx_static)
    _symmap = zeros(Int, im.lastidx_static)
    for batch in batches
        for eidx in batch.indices
            edge = im.edgevec[eidx]

            target = im.v_aggr[edge.dst]
            source = im.e_data[eidx][1:im.edepth]
            _map[source] .= target

            # src mapping
            cplng = coupling(batch)
            if cplng == Symmetric()
                target = im.v_aggr[edge.src]
                source = im.e_data[eidx][1:im.edepth]
                _symmap[source] .= target
            elseif cplng == AntiSymmetric()
                target = -1 .* im.v_aggr[edge.src]
                source = im.e_data[eidx][1:im.edepth]
                _symmap[source] .= target
            elseif cplng == Fiducial()
                target = im.v_aggr[edge.src]
                source = im.e_data[eidx][im.edepth+1:2*im.edepth]
                _symmap[source] .= target
            end
        end
    end
    range, map       = _tighten_idxrange(_map)
    symrange, symmap = _tighten_idxrange(_symmap)
    AggregationMap(range, map, symrange, symmap)
end
function _tighten_idxrange(v)
    first = findfirst(!iszero, v)
    isnothing(first) && return (1:0, similar(v, 0))
    last = findlast(!iszero, v)
    range = first:last
    return range, v[range]
end

struct KAAggregator{F,V} <: Aggregator
    f::F
    m::AggregationMap{V}
end
KAAggregator(f) = (im, batches) -> KAAggregator(im, batches, f)
KAAggregator(im, batches, f) = KAAggregator(f, AggregationMap(im, batches))

function aggregate!(a::KAAggregator, aggbuf, data)
    am = a.m
    fill!(aggbuf, zero(eltype(aggbuf)))
    _backend = get_backend(data)

    # kernel = agg_kernel!(_backend, 1024, length(am.map))
    # kernel(a.f, aggbuf, view(data, am.range), am.map)
    kernel = agg_kernel!(_backend)
    kernel(a.f, aggbuf, view(data, am.range), am.map; ndrange=length(am.map))
    KernelAbstractions.synchronize(_backend)

    if !isempty(am.symrange)
        # symkernel = agg_kernel_sym!(_backend, 1024, length(am.map))
        # symkernel(a.f, aggbuf, view(data, am.symrange), am.symmap)
        symkernel = agg_kernel_sym!(_backend)
        symkernel(a.f, aggbuf, view(data, am.symrange), am.symmap; ndrange=length(am.symmap))
        KernelAbstractions.synchronize(_backend)
    end
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
@kernel function agg_kernel_sym!(f::F, aggbuf, @Const(data), @Const(idxs)) where {F}
    I = @index(Global)
    @inbounds if I ≤ length(idxs)
        dst_idx = idxs[I] # might be < 1 for antisymmetric coupling
        _dst_i = abs(idxs[I])
        if _dst_i != 0
            _dat = sign(dst_idx) * data[I]
            ref = Atomix.IndexableRef(aggbuf, (_dst_i,))
            Atomix.modify!(ref, f, _dat)
        end
    end
    nothing
end


struct SequentialAggregator{F} <: Aggregator
    f::F
    m::AggregationMap{Vector{Int}}
end
SequentialAggregator(f) = (im, batches) -> SequentialAggregator(im, batches, f)
SequentialAggregator(im, batches, f) = SequentialAggregator(f, AggregationMap(im, batches))

function aggregate!(a::SequentialAggregator, aggbuf, data)
    fill!(aggbuf, zero(eltype(aggbuf)))

    @inbounds begin
        am = a.m
        for (dat, dst_idx) in zip(view(data, am.range), am.map)
            if dst_idx != 0
                aggbuf[dst_idx] = a.f(aggbuf[dst_idx], dat)
            end
        end

        for (dat, dst_idx) in zip(view(data, am.symrange), am.symmap)
            if dst_idx != 0
                _dst_idx = abs(dst_idx)
                aggbuf[_dst_idx] = a.f(aggbuf[_dst_idx], sign(dst_idx) * dat)
            end
        end
    end

    nothing
end


struct PolyesterAggregator{F} <: Aggregator
    f::F
    m::Vector{Tuple{Int64,Vector{Int64}}}
end
PolyesterAggregator(f) = (im, batches) -> PolyesterAggregator(im, batches, f)
PolyesterAggregator(im, batches, f) = PolyesterAggregator(f, _inv_aggregation_map(im, batches))

function aggregate!(a::PolyesterAggregator, aggbuf, data)
    length(a.m) == length(aggbuf) || throw(DimensionMismatch("length of aggbuf and a.m must be equal"))
    fill!(aggbuf, zero(eltype(aggbuf)))

    Polyester.@batch for (dstidx, srcidxs) in a.m
       @inbounds for srcidx in srcidxs
           dat = sign(srcidx) * data[abs(srcidx)]
           aggbuf[dstidx] = a.f(aggbuf[dstidx], dat)
       end
    end
    nothing
end


struct ThreadedAggregator{F} <: Aggregator
    f::F
    m::Vector{Tuple{Int64,Vector{Int64}}}
end
ThreadedAggregator(f) = (im, batches) -> ThreadedAggregator(im, batches, f)
ThreadedAggregator(im, batches, f) = ThreadedAggregator(f, _inv_aggregation_map(im, batches))

function aggregate!(a::ThreadedAggregator, aggbuf, data)
    length(a.m) == length(aggbuf) || throw(DimensionMismatch("length of aggbuf and a.m must be equal"))
    fill!(aggbuf, zero(eltype(aggbuf)))

    Threads.@threads for (dstidx, srcidxs) in a.m
       @inbounds for srcidx in srcidxs
           dat = sign(srcidx) * data[abs(srcidx)]
           aggbuf[dstidx] = a.f(aggbuf[dstidx], dat)
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
            edat_idx = im.e_data[eidx][1:im.edepth]
            _pusheach!(srcidxs[edge.dst], edat_idx)

            # src mapping
            cplng = coupling(batch)
            if cplng == Symmetric()
                edat_idx = im.e_data[eidx][1:im.edepth]
                _pusheach!(srcidxs[edge.src], edat_idx)
            elseif cplng == AntiSymmetric()
                edat_idx = -1 .* im.e_data[eidx][1:im.edepth]
                _pusheach!(srcidxs[edge.src], edat_idx)
            elseif cplng == Fiducial()
                edat_idx = im.e_data[eidx][im.edepth+1:2*im.edepth]
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

struct SparseAggregator{M} <: Aggregator
    m::M
    SparseAggregator(m::AbstractMatrix) = new{typeof(m)}(m)
end
function SparseAggregator(f)
    @argcheck f===(+) ArgumentError("Sparse Aggregator only works with + as reducer.")
    SparseAggregator
end
function SparseAggregator(im, batches)
    I, J, V = Float64[], Float64[], Float64[]
    unrolled_foreach(batches) do batch
        for eidx in batch.indices
            edge = im.edgevec[eidx]

            # dst mapping
            edat_idx = im.e_data[eidx][1:im.edepth]
            dst_idx = im.v_aggr[edge.dst]
            append!(I, dst_idx)
            append!(J, edat_idx)
            append!(V, Iterators.repeated(1, length(dst_idx)))

            # src mapping
            cplng = coupling(batch)
            if cplng == Symmetric()
                edat_idx = im.e_data[eidx][1:im.edepth]
                dst_idx = im.v_aggr[edge.src]
                append!(I, dst_idx)
                append!(J, edat_idx)
                append!(V, Iterators.repeated(1, length(dst_idx)))
            elseif cplng == AntiSymmetric()
                edat_idx = im.e_data[eidx][1:im.edepth]
                dst_idx = im.v_aggr[edge.src]
                append!(I, dst_idx)
                append!(J, edat_idx)
                append!(V, Iterators.repeated(-1, length(dst_idx)))
            elseif cplng == Fiducial()
                edat_idx = im.e_data[eidx][im.edepth+1:2*im.edepth]
                dst_idx = im.v_aggr[edge.src]
                append!(I, dst_idx)
                append!(J, edat_idx)
                append!(V, Iterators.repeated(1, length(dst_idx)))
            end
        end
    end

    SparseAggregator(sparse(I,J,V, im.lastidx_aggr, im.lastidx_static))
end

function aggregate!(a::SparseAggregator, aggbuf, data)
   LinearAlgebra.mul!(aggbuf, a.m, data)
   nothing
end
