export Network

function Network(g::SimpleGraph,
                 vertexf::Union{VertexFunction, Vector{<:VertexFunction}},
                 edgef::Union{EdgeFunction, Vector{<:EdgeFunction}};
                 accumulator=+, accdim=:auto,
                 execution=:seq,
                 verbose=false)
    verbose && println("Create dynamic network with $(nv(g)) vertices and $(ne(g)) edges:")
    @check execution ∈ EXECUTION_STYLES "Exectuion type $execution not supportet (choose from $(EXECUTION_STYLES))"

    im = IndexManager()

    vtypes, idxs = batch_identical(vertexf, collect(1:nv(g)))
    vertexbatches = collect(VertexBatch(im, i, t; verbose) for (i, t) in zip(idxs, vtypes))

    nl = NetworkLayer(im, g, edgef, accumulator, accdim; verbose)

    @assert isdense(im) "The index structure is not dense"
    return Network{execution, typeof(nl), typeof(vertexbatches)}(vertexbatches, nl, im)
end

function VertexBatch(im::IndexManager,
                     vertex_idxs::Vector{Int},
                     vertexf::VertexFunction;
                     verbose)
    dim = vertexf.dim
    pdim = vertexf.pdim
    (firstidx, pfirstidx) = register_vertices!(im, vertex_idxs, dim, pdim)

    verbose && println(" - VertexBatch: dim=$dim, pdim=$pdim, length=$(length(vertex_idxs))")

    VertexBatch(vertex_idxs,
                vertexf,
                dim,
                firstidx,
                pdim,
                pfirstidx)
end

function NetworkLayer(im::IndexManager, g::SimpleGraph,
                      edgef::Union{EdgeFunction, Vector{<:EdgeFunction}},
                      accumulator, accdim;
                      verbose)
    @check hasmethod(accumulator, NTuple{2, AbstractFloat}) "Accumulator needs `acc(AbstractFloat, AbstractFloat) method`"

    accdim_max = maxaccdim(edgef)
    if accdim === :auto
        accdim = accdim_max
        verbose && println(" - auto accumulation dimension = $accdim")
    end
    @check accdim<=accdim_max "For this system acc dim is limited to $accdim_max by the edgefunctions"

    gc = color_edges_greedy(g)
    verbose && println(" - found $(length(uniquecolors(gc))) edgecolors (optimum would be $(maximum(degree(g))))")

    edgecolors = [color(gc, e) for e in edges(g)]
    colors = unique(edgecolors)

    idx_per_color = [findall(isequal(c), edgecolors) for c in colors]

    edgef_per_color = batch_by_idxs(edgef, idx_per_color)

    colorbatches = collect(ColorBatch(im, g, idx, ef; verbose)
                         for (idx, ef) in zip(idx_per_color, edgef_per_color))

    NetworkLayer(g, colorbatches, accumulator, accdim, CachePool())
end

maxaccdim(e::EdgeFunction) = accdim(e)
maxaccdim(e::AbstractVector{<:EdgeFunction}) = minimum(accdim.(e))

function ColorBatch(im::IndexManager, g::SimpleGraph,
                    edge_idxs::Vector{Int},
                    edgef::Union{EdgeFunction, Vector{<:EdgeFunction}};
                    verbose)
    verbose && println(" - ColorBatch with $(length(edge_idxs)) edges")

    eftypes, idxs = batch_identical(edgef, edge_idxs)
    edgebatches = collect(EdgeBatch(im, g, i, ef; verbose) for (i, ef) in zip(idxs, eftypes))

    return ColorBatch(edge_idxs, edgebatches)
end

function EdgeBatch(im::IndexManager, g::SimpleGraph,
                   edge_idxs::Vector{Int},
                   edgef::EdgeFunction;
                   verbose)
    dim = edgef.dim
    pdim = edgef.pdim

    if edgef isa StaticEdge
        (firstidx, pfirstidx) = register_edges!(im, edge_idxs, 0, pdim)
    else
        (firstidx, pfirstidx) = register_edges!(im, edge_idxs, dim, pdim)
    end

    vidx_array = Array{Int}(undef, 6*length(edge_idxs))
    empty!(vidx_array)
    alledges = collect(edges(g))
    for (i, eidx) in enumerate(edge_idxs)
        edge = alledges[eidx]
        src, dst = edge.src, edge.dst
        src_r = vertex_data_range(im, src)
        dst_r = vertex_data_range(im, dst)
        src_idx, dst_idx = src_r[begin], dst_r[begin]
        src_dim, dst_dim = length(src_r), length(dst_r)

        append!(vidx_array, src, src_idx, src_dim, dst, dst_idx, dst_dim)
    end
    @assert length(vidx_array) == 6*length(edge_idxs)

    verbose && println("   - EdgeBatch: dim=$dim, pdim=$pdim, length=$(length(edge_idxs))")

    EdgeBatch(edge_idxs,
              edgef,
              dim,
              firstidx,
              pdim,
              pfirstidx,
              vidx_array)
end

batch_by_idxs(v, idxs::Vector{Vector{Int}}) = [v for batch in idxs]
function batch_by_idxs(v::AbstractVector, batches::Vector{Vector{Int}})
    @assert length(v) == sum(length.(batches))
    [v[batch] for batch in batches]
end

batch_identical(el, idxs::Vector{Int}) = [el], [idxs]
function batch_identical(v::Vector{T}, indices::Vector{Int}) where {T}
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
