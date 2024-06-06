abstract type Aggregator end

struct NNlibScatter{T, V} <: Aggregator
    f::T
    batchranges::Vector{UnitRange{Int}}
    couplings::Vector{CouplingUnion}
    dstmaps::Vector{V}
    srcmaps::Vector{V}
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
        cplng = coupling(batch)
        # generate scatter map
        dst = Vector{Int}(undef, length(state_range(batch)))
        fill!(dst, -1)
        src = Vector{Int}(undef, length(state_range(batch)))
        fill!(src, -1)

        for i in batch.indices
            datarange = im.e_data[i] .- (batch.statestride.first - 1) # range in batch slice
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
        @assert get_backend(aggbuf) == get_backend(dstmap)
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
    m::AggregationMap{Vector{Int}}
end
PolyesterAggregator(f) = (im, batches) -> PolyesterAggregator(im, batches, f)
PolyesterAggregator(im, batches, f) = PolyesterAggregator(f, AggregationMap(im, batches))

function aggregate!(a::PolyesterAggregator, aggbuf, data)
    fill!(aggbuf, zero(eltype(aggbuf)))

    let _data=view(data, a.m.range), idxs=a.m.map, f=a.f
        Polyester.@batch for I in 1:length(idxs)
            @inbounds begin
                _dst_i = idxs[I]
                if _dst_i != 0
                    _dat = _data[I]
                    ref = Atomix.IndexableRef(aggbuf, (_dst_i,))
                    Atomix.modify!(ref, f, _dat)
                end
            end
        end
    end

    let _data=view(data, a.m.symrange), idxs=a.m.symmap, f=a.f
        Polyester.@batch for I in 1:length(idxs)
            @inbounds begin
                dst_idx = idxs[I] # might be < 1 for antisymmetric coupling
                _dst_i = abs(idxs[I])
                if _dst_i != 0
                    _dat = sign(dst_idx) * _data[I]
                    ref = Atomix.IndexableRef(aggbuf, (_dst_i,))
                    Atomix.modify!(ref, f, _dat)
                end
            end
        end
    end

    nothing
end
