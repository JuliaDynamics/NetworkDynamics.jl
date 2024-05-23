function Network(g::AbstractGraph,
                 vertexf::Union{VertexFunction,Vector{<:VertexFunction}},
                 edgef::Union{EdgeFunction,Vector{<:EdgeFunction}};
                 accumulator=NNlibScatter(+),
                 edepth=:auto,
                 vdepth=:auto,
                 execution=SequentialExecution{true}(),
                 verbose=false)
    reset_timer!()
    @timeit_debug "Construct Network" begin
        # collect all vertex/edgf to vector
        _vertexf = vertexf isa Vector ? vertexf : [vertexf for _ in vertices(g)]
        _edgef = edgef isa Vector ? edgef : [edgef for _ in edges(g)]
        @assert _vertexf isa Vector{<:VertexFunction}
        @assert _edgef isa Vector{<:EdgeFunction}

        verbose &&
            println("Create dynamic network with $(nv(g)) vertices and $(ne(g)) edges:")
        @argcheck execution isa ExecutionStyle "Exectuion type $execution not supportet (choose from $(subtypes(ExecutionStyle)))"

        _maxedepth = mapreduce(aggrdepth, min, _edgef)
        if edepth === :auto
            edepth = _maxedepth
            verbose && println(" - auto accumulation depth = $edepth")
        end
        @argcheck edepth<=_maxedepth "For this system accumulation depth is limited to $edepth by the edgefunctions"

        _maxvdepth = mapreduce(dim, min, _vertexf)
        if vdepth === :auto
            vdepth = _maxvdepth
            verbose && println(" - auto edge input depth = $vdepth")
        end
        @argcheck vdepth<=_maxvdepth "For this system edge input depth is limited to $edepth by the vertex dimensions"

        dynstates = mapreduce(+, Iterators.flatten((_vertexf,_edgef))) do t
            statetype(t) == Dynamic() ? dim(t) : 0
        end

        # create index manager
        im = IndexManager(g, dynstates, edepth, vdepth, _vertexf, _edgef)

        # batch identical edge and vertex functions
        vidxs = _find_identical(_vertexf)
        eidxs = _find_identical(_edgef)

        # create vertex batches and initialize with index manager
        @timeit_debug "create vertex batches" begin
            vertexbatches = map(vidxs) do idxs
                VertexBatch(im, idxs; verbose)
            end

            if length(vertexbatches) ≤ 50
                vertexbatches = Tuple(vertexbatches)
            else
                verbose &&
                    println(" $(length(vertexbatches)) > 50 unique indices: don't unroll!")
            end
        end

        # create edge batches and initialize with index manager
        @timeit_debug "create edge batches" begin
            edgebatches = map(eidxs) do idxs
                EdgeBatch(im, idxs; verbose)
            end
            if length(edgebatches) ≤ 50
                edgebatches = Tuple(edgebatches)
            else
                verbose &&
                    println("$(length(edgebatches)) > 50 unique edges: don't unroll!")
            end
        end

        @timeit_debug "initialize aggregator" begin
            aggregator = accumulator(im, edgebatches)
        end

        nl = NetworkLayer(im, edgebatches, aggregator)

        @assert isdense(im)
        nw = Network{typeof(execution),typeof(g),typeof(nl),typeof(vertexbatches)}(vertexbatches,
                                                                                   nl, im,
                                                                                   LazyBufferCache())

    end
    # print_timer()
    return nw
end

function VertexBatch(im::IndexManager, idxs::Vector{Int}; verbose)
    components = @view im.vertexf[idxs]

    _compT = compT(only(unique(compT, components)))
    _compf = compf(only(unique(compf, components)))
    _statetype = statetype(only(unique(statetype, components)))
    _dim = dim(only(unique(dim, components)))
    _pdim = pdim(only(unique(pdim, components)))

    (statestride, pstride, aggbufstride) = register_vertices!(im, _statetype, _dim, _pdim, idxs)

    verbose &&
        println(" - VertexBatch: dim=$(_dim), pdim=$(_pdim), length=$(length(idxs))")

    VertexBatch{_compT, typeof(_compf)}(idxs, _compf, statestride, pstride, aggbufstride)
end

function EdgeBatch(im::IndexManager, idxs::Vector{Int}; verbose)
    components = @view im.edgef[idxs]

    _compT = compT(only(unique(compT, components)))
    _compf = compf(only(unique(compf, components)))
    _statetype = statetype(only(unique(statetype, components)))
    _dim = dim(only(unique(dim, components)))
    _pdim = pdim(only(unique(pdim, components)))

    (statestride, pstride, gbufstride) = register_edges!(im, _statetype, _dim, _pdim, idxs)

    verbose &&
        println(" - EdgeBatch: dim=$(_dim), pdim=$(_pdim), length=$(length(idxs))")

    EdgeBatch{_compT, typeof(_compf)}(idxs, _compf, statestride, pstride, gbufstride)
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

function _find_identical(v::Vector)
    indices = eachindex(v)
    idxs_per_type = Vector{Int}[]
    unique_compf = []
    for i in eachindex(v)
        found = false
        for j in eachindex(unique_compf)
            if compf(v[i]) == unique_compf[j]
                found = true
                push!(idxs_per_type[j], indices[i])
                break
            end
        end
        if !found
            push!(unique_compf, compf(v[i]))
            push!(idxs_per_type, [indices[i]])
        end
    end
    @assert length(unique_compf) == length(idxs_per_type)
    return idxs_per_type
end
