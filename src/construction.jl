export Network

function Network(g::SimpleGraph,
                 vertexf::Union{VertexFunction, Vector{<:VertexFunction}},
                 edgef::Union{EdgeFunction, Vector{<:EdgeFunction}};
                 accumulator=+,
                 edepth=:auto,
                 vdepth=:auto,
                 execution=SequentialExecution(),
                 verbose=false)
    verbose && println("Create dynamic network with $(nv(g)) vertices and $(ne(g)) edges:")
    @argcheck execution isa ExecutionStyle "Exectuion type $execution not supportet (choose from $(subtypes(ExecutionStyle)))"

    im = IndexManager()

    vtypes, idxs = _batch_identical(vertexf, collect(1:nv(g)))
    vertexbatches = collect(VertexBatch(im, i, t; verbose) for (i, t) in zip(idxs, vtypes))

    nl = NetworkLayer(im, g, vertexf, edgef, accumulator, edepth, vdepth; verbose)

    @assert isdense(im)
    return Network{typeof(execution), typeof(nl), typeof(vertexbatches)}(vertexbatches, nl, im, LazyBufferCache())
end

function VertexBatch(im::IndexManager,
                     vertex_idxs::Vector{Int},
                     vertexf::VertexFunction;
                     verbose)
    (firstidx, pfirstidx) = register_components!(im, vertex_idxs, vertexf)

    dim = vertexf.dim
    pdim = vertexf.pdim
    verbose && println(" - VertexBatch: dim=$(dim), pdim=$(pdim), length=$(length(vertex_idxs))")

    VertexBatch(vertex_idxs, vertexf, dim, firstidx, pdim, pfirstidx)
end

function NetworkLayer(im::IndexManager, g::SimpleGraph,
                      vertexf::Union{VertexFunction, Vector{<:VertexFunction}},
                      edgef::Union{EdgeFunction, Vector{<:EdgeFunction}},
                      accumulator, edepth, vdepth;
                      verbose)
    @argcheck hasmethod(accumulator, NTuple{2, AbstractFloat}) "Accumulator needs `acc(AbstractFloat, AbstractFloat) method`"

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

    NetworkLayer(g, edgebatches, accumulator, edepth, vdepth)
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
    verbose && println(" - EdgeBatch: dim=$(dim), pdim=$(pdim), length=$(length(edge_idxs))")

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
