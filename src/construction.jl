export Network

function Network(g::SimpleGraph,
                 vertexf::Union{VertexFunction,Vector{<:VertexFunction}},
                 edgef::Union{EdgeFunction,Vector{<:EdgeFunction}};
                 accumulator=+,
                 edepth=:auto,
                 vdepth=:auto,
                 execution=SequentialExecution{true}(),
                 verbose=false)
    verbose && println("Create dynamic network with $(nv(g)) vertices and $(ne(g)) edges:")
    @argcheck execution isa ExecutionStyle "Exectuion type $execution not supportet (choose from $(subtypes(ExecutionStyle)))"

    im = IndexManager()

    vtypes, idxs = _batch_identical(vertexf, collect(1:nv(g)))
    vertexbatches = collect(VertexBatch(im, i, t; verbose) for (i, t) in zip(idxs, vtypes))

    nl = NetworkLayer(im, g, vertexf, edgef, accumulator, edepth, vdepth, execution;
                      verbose)

    @assert isdense(im)
    nw = Network{typeof(execution),typeof(nl),typeof(vertexbatches)}(vertexbatches, nl, im,
                                                                     LazyBufferCache())

    init_gathermap!(nw)
    init_scattermap!(nw)

    return nw
end

function VertexBatch(im::IndexManager,
                     vertex_idxs::Vector{Int},
                     vertexf::VertexFunction;
                     verbose)
    (firstidx, pfirstidx) = register_components!(im, vertex_idxs, vertexf)

    dim = vertexf.dim
    pdim = vertexf.pdim
    verbose &&
        println(" - VertexBatch: dim=$(dim), pdim=$(pdim), length=$(length(vertex_idxs))")

    VertexBatch(vertex_idxs, vertexf, dim, firstidx, pdim, pfirstidx)
end

function NetworkLayer(im::IndexManager,
                      g::SimpleGraph,
                      vertexf::Union{VertexFunction,Vector{<:VertexFunction}},
                      edgef::Union{EdgeFunction,Vector{<:EdgeFunction}},
                      accumulator,
                      edepth, vdepth,
                      execution;
                      verbose)
    @argcheck hasmethod(accumulator, NTuple{2,AbstractFloat}) "Accumulator needs `acc(AbstractFloat, AbstractFloat) method`"

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

    etypes, idxs = _batch_identical(edgef, collect(1:ne(g)))
    edgebatches = collect(EdgeBatch(im, i, t; verbose) for (i, t) in zip(idxs, etypes))

    # edgecolors = color_edges_greedy(g)
    # verbose && println(" - found $(length(unique(edgecolors))) edgecolors (optimum would be $(maximum(degree(g))))")

    # colors = unique(edgecolors)

    # idx_per_color = [findall(isequal(c), edgecolors) for c in colors]

    # edgef_per_color = batch_by_idxs(edgef, idx_per_color)

    # colorbatches = collect(ColorBatch(im, g, idx, ef; verbose)
    #                      for (idx, ef) in zip(idx_per_color, edgef_per_color))

    NetworkLayer(g, edgebatches, accumulator, edepth, vdepth, execution)
end

maxedepth(e::EdgeFunction) = accdepth(e)
maxedepth(e::AbstractVector{<:EdgeFunction}) = minimum(accdepth.(e))
maxvdepth(v::VertexFunction) = dim(v)
maxvdepth(v::AbstractVector{<:EdgeFunction}) = minimum(dim.(v))

function EdgeBatch(im::IndexManager,
                   edge_idxs::Vector{Int},
                   edgef::EdgeFunction;
                   verbose)
    (firstidx, pfirstidx) = register_components!(im, edge_idxs, edgef)

    dim = edgef.dim
    pdim = edgef.pdim
    verbose &&
        println(" - EdgeBatch: dim=$(dim), pdim=$(pdim), length=$(length(edge_idxs))")

    EdgeBatch(edge_idxs, edgef, dim, firstidx, pdim, pfirstidx)
end

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

function init_scattermap!(nw::Network)
    # just one layer currently
    for batch in nw.nl.edgebatches
        init_scattermap!(coupling(batch), nw, nw.nl, batch)
    end
end

function init_scattermap!(::Union{AntiSymmetric,Symmetric}, nw, layer, batch)
    edepth = layer.edepth
    nullelement = nv(layer.g) * edepth + 1
    resize!(batch.scatter_map[1], length(batch.indices) * dim(batch.fun))
    resize!(batch.scatter_map[2], length(batch.indices) * dim(batch.fun))
    fill!(batch.scatter_map[1], nullelement)
    fill!(batch.scatter_map[2], nullelement)
    for (i, eidx) in enumerate(batch.indices)
        (; src, dst) = edgebyidx(layer.g, eidx)

        src_start = (src - 1) * edepth + 1
        dst_start = (src - 1) * edepth + 1
        src_range = src_start:(src_start+edepth-1)
        dst_range = dst_start:(dst_start+edepth-1)

        startidx = (i - 1) * dim(batch.fun) + 1
        idxrange = startidx:(startidx+edepth-1)
        batch.scatter_map[1][idxrange] .= dst_range
        batch.scatter_map[2][idxrange] .= src_range
    end
end
function init_scattermap!(::Fiducial, nw, layer, batch)
    edepth = layer.edepth
    nullelement = nv(layer.g) * edepth + 1
    resize!(batch.scatter_map, length(batch.indices) * dim(batch.fun))
    fill!(batch.scatter_map, nullelement)
    for (i, eidx) in enumerate(batch.indices)
        (; src, dst) = edgebyidx(layer.g, eidx)

        src_start = (src - 1) * edepth + 1
        dst_start = (src - 1) * edepth + 1
        src_range = src_start:(src_start+edepth-1)
        dst_range = dst_start:(dst_start+edepth-1)

        startidx = (i - 1) * dim(batch.fun) + 1
        idxrange = startidx:(startidx+2*edepth-1)
        batch.scatter_map[idxrange] .= vcat(dst_range, src_range)
    end
end
function init_scattermap!(::Directed, nw, layer, batch)
    edepth = layer.edepth
    nullelement = nv(layer.g) * edepth + 1
    resize!(batch.scatter_map, length(batch.indices) * dim(batch.fun))
    fill!(batch.scatter_map, nullelement)
    for (i, eidx) in enumerate(batch.indices)
        (; dst) = edgebyidx(layer.g, eidx)

        dst_start = (src - 1) * edepth + 1
        dst_range = dst_start:(dst_start+edepth-1)

        startidx = (i - 1) * dim(batch.fun) + 1
        idxrange = startidx:(startidx+edepth-1)
        batch.scatter_map[idxrange] .= dst_range
    end
end

# buffered version
function init_gathermap!(nw::Network{<:ExecutionStyle{false}})
    layer = nw.nl
    for (i, e) in enumerate(edges(layer.g))
        layer.gather_map[1][i] = nw.im.v_data[e.dst][2]
        layer.gather_map[2][i] = nw.im.v_data[e.src][2]
    end
end

# unbuffered version
function init_gathermap!(nw::Network{<:ExecutionStyle{true}})
    layer = nw.nl
    vdepth = layer.vdepth
    for (i, e) in enumerate(edges(layer.g))
        startidx = (i - 1) * vdepth + 1
        range = startidx:(startidx+vdepth-1)
        dst_range = nw.im.v_data[e.src][2][1:vdepth]
        src_range = nw.im.v_data[e.dst][2][1:vdepth]
        layer.gather_map[range, 1] .= dst_range
        layer.gather_map[range, 2] .= src_range
    end
end
