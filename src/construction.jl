export Network

function Network(g::AbstractGraph,
                 vertexf::Union{VertexFunction,Vector{<:VertexFunction}},
                 edgef::Union{EdgeFunction,Vector{<:EdgeFunction}};
                 accumulator=NaiveAggregator(+),
                 edepth=:auto,
                 vdepth=:auto,
                 execution=SequentialExecution{true}(),
                 verbose=false)
    verbose && println("Create dynamic network with $(nv(g)) vertices and $(ne(g)) edges:")
    @argcheck execution isa ExecutionStyle "Exectuion type $execution not supportet (choose from $(subtypes(ExecutionStyle)))"

    _maxedepth = maxedepth(edgef)
    if edepth === :auto
        edepth = _maxedepth
        verbose && println(" - auto accumulation depth = $edepth")
    end
    @argcheck edepth<=_maxedepth "For this system accumulation depth is limited to $edepth by the edgefunctions"

    _maxvdepth = maxvdepth(vertexf)
    if vdepth === :auto
        vdepth = _maxvdepth
        verbose && println(" - auto edge input depth = $vdepth")
    end
    @argcheck vdepth<=_maxvdepth "For this system edge input depth is limited to $edepth by the vertex dimensions"

    # batch identical edge and vertex functions
    vtypes, vidxs = _batch_identical(vertexf, collect(1:nv(g)))
    etypes, eidxs = _batch_identical(edgef, collect(1:ne(g)))

    # count dynamic states
    dynstates = 0
    for (t, idxs) in zip(vtypes, vidxs)
        if statetype(t) == Dynamic()
            dynstates += length(idxs) * dim(t)
        end
    end
    for (t, idxs) in zip(etypes, eidxs)
        if statetype(t) == Dynamic()
            dynstates += length(idxs) * dim(t)
        end
    end

    # create index manager
    im = IndexManager(g, dynstates, edepth, vdepth)

    # create vertex batches and initialize with index manager
    vertexbatches = collect(VertexBatch(im, i, t; verbose) for (i, t) in zip(vidxs, vtypes))

    # create edge batches and initialize with index manager
    edgebatches = collect(EdgeBatch(im, i, t; verbose) for (i, t) in zip(eidxs, etypes))

    aggregator = accumulator(im, edgebatches)

    nl = NetworkLayer(im, edgebatches, aggregator)

    @assert isdense(im)
    nw = Network{typeof(execution),typeof(nl),typeof(vertexbatches)}(vertexbatches, nl, im,
                                                                     LazyBufferCache())

    return nw
end

function VertexBatch(im::IndexManager,
                     idxs::Vector{Int},
                     vertexf::VertexFunction;
                     verbose)
    (firstidx, pfirstidx, aggrfirstidx) = register_vertices!(im, idxs, vertexf)

    verbose &&
        println(" - VertexBatch: dim=$(vertexf.dim), pdim=$(vertexf.pdim), length=$(length(idxs))")

    VertexBatch(idxs, vertexf,
                vertexf.dim, firstidx,
                vertexf.pdim, pfirstidx,
                im.edepth, aggrfirstidx)
end

function EdgeBatch(im::IndexManager,
                   idxs::Vector{Int},
                   edgef::EdgeFunction;
                   verbose)
    (firstidx, pfirstidx) = register_edges!(im, idxs, edgef)

    verbose &&
        println(" - EdgeBatch: dim=$(edgef.dim), pdim=$(edgef.pdim), length=$(length(idxs))")

    EdgeBatch(idxs, edgef,
              edgef.dim, firstidx,
              edgef.pdim, pfirstidx)
end

function NetworkLayer(im::IndexManager, eb, agg)
    map = zeros(Int, ne(im.g) * im.vdepth, 2)
    for (i, e) in enumerate(edges(im.g))
        startidx = (i - 1) * im.vdepth + 1
        range = startidx:(startidx+im.vdepth-1)
        dst_range = im.v_data[e.src][1:im.vdepth]
        src_range = im.v_data[e.dst][1:im.vdepth]
        map[range, 1] .= dst_range
        map[range, 2] .= src_range
    end
    NetworkLayer(im.g, eb, agg, im.edepth, im.vdepth, map)
end

maxedepth(e::EdgeFunction) = aggrdepth(e)
maxedepth(e::AbstractVector{<:EdgeFunction}) = minimum(aggrdepth.(e))
maxvdepth(v::VertexFunction) = dim(v)
maxvdepth(v::AbstractVector{<:VertexFunction}) = minimum(dim.(v))

batch_by_idxs(v, idxs::Vector{Vector{Int}}) = [v for batch in idxs]
function batch_by_idxs(v::AbstractVector, batches::Vector{Vector{Int}})
    @assert length(v) == sum(length.(batches))
    [v[batch] for batch in batches]
end

_batch_identical(el, idxs::Vector{Int}) = [el], [idxs]
function _batch_identical(v::Vector{T}, indices::Vector{Int}) where {T}
    @assert length(v) == length(indices)
    idxs_per_type = Vector{Int}[]
    types = T[]
    for i in eachindex(v)
        found = false
        for j in eachindex(types)
            if v[i] === types[j]
                found = true
                push!(idxs_per_type[j], indices[i])
                break
            end
        end
        if !found
            push!(types, v[i])
            push!(idxs_per_type, [indices[i]])
        end
    end
    @assert length(types) == length(idxs_per_type)
    return types, idxs_per_type
end
