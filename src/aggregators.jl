abstract type Aggregator end

struct NNlibScatter{T} <: Aggregator
    f::T
    batchranges::Vector{UnitRange{Int}}
    couplings::Vector{CouplingUnion}
    dstmaps::Vector{Vector{Int}}
    srcmaps::Vector{Vector{Int}}
    aggrsize::Int
end

NNlibScatter(f) = (im, batches) -> NNlibScatter(im, batches, f)
function NNlibScatter(im, batches, f)
    batchranges = Vector{UnitRange{Int}}(undef, length(batches))
    dstmaps     = Vector{Vector{Int}}(undef, length(batches))
    srcmaps     = Vector{Vector{Int}}(undef, length(batches))
    couplings   = Vector{CouplingUnion}(undef, length(batches))
    maxaggindex = 0
    _edgevec    = collect(edges(im.g))
    for (batchi, batch) in enumerate(batches)
        cplng = coupling(batch.fun)
        # generate scatter map
        dst = Vector{Int}(undef, length(batch.indices) * dim(batch.fun))
        fill!(dst, -1)
        src = Vector{Int}(undef, length(batch.indices) * dim(batch.fun))
        fill!(src, -1)

        for i in batch.indices
            datarange = im.e_data[i] .- (batch.firstidx - 1) # range in batch slice
            e = _edgevec[i]
            dst[datarange[1:im.edepth]] .= im.v_aggr[e.dst]
            if cplng == Symmetric() || cplng == AntiSymmetric()
                src[datarange[1:im.edepth]] .= im.v_aggr[e.src]
            elseif cplng == Fiducial()
                src[datarange[im.edepth+1:2*im.edepth]] .= im.v_aggr[e.src]
            end
        end

        # unasigned fields will have a -1, number them up to copy to "ghost elements" without access conflict
        unassigned = findall(isequal(-1), dst)
        _maxaggindex = im.lastidx_aggr + length(unassigned)
        dst[unassigned] .= im.lastidx_aggr+1:_maxaggindex
        maxaggindex = max(maxaggindex, _maxaggindex)
        if cplng != Directed()
            unassigned = findall(isequal(-1), src)
            _maxaggindex = im.lastidx_aggr + length(unassigned)
            src[unassigned] .= im.lastidx_aggr+1:_maxaggindex
            maxaggindex = max(maxaggindex, _maxaggindex)
        end

        batchranges[batchi] = state_range(batch)
        couplings[batchi]   = cplng
        dstmaps[batchi]     = dst
        srcmaps[batchi]     = src
    end
    NNlibScatter(f, batchranges, couplings, dstmaps, srcmaps, maxaggindex)
end

function aggregate!(a::NNlibScatter, aggbuf, data)
    fill!(aggbuf, zero(eltype(aggbuf)))
    originallength = length(aggbuf)
    resize!(aggbuf, a.aggrsize)
    @assert length(aggbuf) == a.aggrsize
    for (range, dstmap, srcmap, coupling) in zip(a.batchranges, a.dstmaps, a.srcmaps,
                                                 a.couplings)
        batchdata = view(data, range)
        # scatter the dst
        NNlib.scatter!(a.f, aggbuf, batchdata, dstmap)

        # scatter to source depending on
        if coupling ∈ (Symmetric(), Fiducial())
            NNlib.scatter!(a.f, aggbuf, batchdata, srcmap)
        elseif coupling == AntiSymmetric()
            batchdata .= -1 .* batchdata
            NNlib.scatter!(a.f, aggbuf, batchdata, srcmap)
            batchdata .= -1 .* batchdata
        end
    end
    resize!(aggbuf, originallength)
    nothing
end


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
            cplng = coupling(batch.fun)
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


struct KAAggregator{F} <: Aggregator
    f::F
    range::UnitRange{Int} # range in data where aggregation is necessary
    map::Vector{Int}      # maps data idx to destination, if zero skip
    symrange::UnitRange{Int}     # same for symmetric/antisymmetric coupling
    symmap::Vector{Int}
end
KAAggregator(f) = (im, batches) -> KAAggregator(im, batches, f)
function KAAggregator(im, batches, f)
    _map    = zeros(Int, im.lastidx_static)
    _symmap = zeros(Int, im.lastidx_static)
    for batch in batches
        for eidx in batch.indices
            edge = im.edgevec[eidx]

            target = im.v_aggr[edge.dst]
            source = im.e_data[eidx][1:im.edepth]
            _map[source] .= target

            # src mapping
            cplng = coupling(batch.fun)
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
    KAAggregator(f, range, map, symrange, symmap)
end
function _tighten_idxrange(v)
    first = findfirst(!iszero, v)
    isnothing(first) && return (1:0, similar(v, 0))
    last = findlast(!iszero, v)
    range = first:last
    return range, v[range]
end

function aggregate!(a::KAAggregator, aggbuf, data)
    fill!(aggbuf, zero(eltype(aggbuf)))
    _backend = get_backend(data)

    # kernel = agg_kernel!(_backend, 1024, length(a.map))
    # kernel(a.f, aggbuf, view(data, a.range), a.map)
    kernel = agg_kernel!(_backend)
    kernel(a.f, aggbuf, view(data, a.range), a.map; ndrange=length(a.map))
    KernelAbstractions.synchronize(_backend)

    if !isempty(a.symrange)
        # symkernel = agg_kernel_sym!(_backend, 1024, length(a.map))
        # symkernel(a.f, aggbuf, view(data, a.symrange), a.symmap)
        symkernel = agg_kernel_sym!(_backend)
        symkernel(a.f, aggbuf, view(data, a.symrange), a.symmap; ndrange=length(a.symmap))
        KernelAbstractions.synchronize(_backend)
    end
    nothing
end

@kernel function agg_kernel!(f::F, aggbuf, data, idxs) where {F}
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
@kernel function agg_kernel_sym!(f::F, aggbuf, data, idxs) where {F}
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
    kaagg::KAAggregator{F}
end
function SequentialAggregator(f)
    (im, batches) -> SequentialAggregator(KAAggregator(im, batches, f))
end

function aggregate!(sa::SequentialAggregator, aggbuf, data)
    fill!(aggbuf, zero(eltype(aggbuf)))

    @inbounds begin
        a = sa.kaagg
        for (dat, dst_idx) in zip(view(data, a.range), a.map)
            aggbuf[dst_idx] = a.f(aggbuf[dst_idx], dat)
        end

        for (dat, dst_idx) in zip(view(data, a.symrange), a.symmap)
            _dst_idx = abs(dst_idx)
            aggbuf[_dst_idx] = a.f(aggbuf[_dst_idx], sign(dst_idx) * dat)
        end
    end

    nothing
end
