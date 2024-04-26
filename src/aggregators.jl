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
    for (batchi, batch) in enumerate(batches)
        cplng = coupling(batch.fun)
        # generate scatter map
        dst = Vector{Int}(undef, length(batch.indices) * dim(batch.fun))
        fill!(dst, -1)
        src = Vector{Int}(undef, length(batch.indices) * dim(batch.fun))
        fill!(src, -1)

        for i in batch.indices
            datarange = im.e_data[i] .- (batch.firstidx - 1) # range in batch slice
            e = edgebyidx(im.g, i)
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
    originallength = length(aggbuf)
    resize!(aggbuf, a.aggrsize)
    @assert length(aggbuf) == a.aggrsize
    fill!(aggbuf, zero(eltype(aggbuf)))
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
@unroll function _aggregate!(a::NaiveAggregator, batches, aggbuf, data)
    im = a.im
    @unroll for batch in batches
        for eidx in batch.indices
            edge = edgebyidx(im.g, eidx)
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
