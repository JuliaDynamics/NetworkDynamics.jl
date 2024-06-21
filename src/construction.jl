function Network(g::AbstractGraph,
                 vertexf::Union{VertexFunction,Vector{<:VertexFunction}},
                 edgef::Union{EdgeFunction,Vector{<:EdgeFunction}};
                 execution=SequentialExecution{true}(),
                 edepth=:auto,
                 vdepth=:auto,
                 aggregator=execution isa SequentialExecution ? SequentialAggregator(+) : PolyesterAggregator(+),
                 verbose=false)
    reset_timer!()
    @timeit_debug "Construct Network" begin
        # collect all vertex/edgf to vector
        _vertexf = vertexf isa Vector ? vertexf : [vertexf for _ in vertices(g)]
        _edgef = edgef isa Vector ? edgef : [edgef for _ in edges(g)]
        @argcheck _vertexf isa Vector{<:VertexFunction} "Expected VertexFuncions, got $(eltype(_vertexf))"
        @argcheck _edgef isa Vector{<:EdgeFunction} "Expected EdgeFuncions, got $(eltype(_vertexf))"
        @argcheck length(_vertexf) == nv(g)
        @argcheck length(_edgef) == ne(g)

        verbose &&
            println("Create dynamic network with $(nv(g)) vertices and $(ne(g)) edges:")
        @argcheck execution isa ExecutionStyle "Exectuion type $execution not supportet (choose from $(subtypes(ExecutionStyle)))"

        _maxedepth = isempty(_edgef) ? 0 : mapreduce(depth, min, _edgef)
        if edepth === :auto
            edepth = _maxedepth
            verbose && println(" - auto accumulation depth = $edepth")
        end
        @argcheck edepth<=_maxedepth "For this system accumulation depth is limited to $edepth by the edgefunctions"

        _maxvdepth = mapreduce(depth, min, _vertexf)
        if vdepth === :auto
            vdepth = _maxvdepth
            verbose && println(" - auto edge input depth = $vdepth")
        end
        @argcheck vdepth<=_maxvdepth "For this system edge input depth is limited to $edepth by the vertex dimensions"

        dynstates = mapreduce(+, Iterators.flatten((_vertexf,_edgef))) do t
            isdynamic(t) ? dim(t) : 0
        end

        # create index manager
        im = IndexManager(g, dynstates, edepth, vdepth, _vertexf, _edgef)

        # batch identical edge and vertex functions
        @timeit_debug "batch identical vertexes" begin
            vidxs = _find_identical(vertexf, 1:nv(g))
        end
        @timeit_debug "batch identical edges" begin
            eidxs = _find_identical(edgef, 1:ne(g))
        end

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
            _aggregator = aggregator(im, edgebatches)
        end

        nl = NetworkLayer(im, edgebatches, _aggregator)

        @assert isdense(im)
        mass_matrix = construct_mass_matrix(im)
        nw = Network{typeof(execution),typeof(g),typeof(nl),typeof(vertexbatches),typeof(mass_matrix)}(
            vertexbatches,
            nl, im,
            LazyBufferCache(),
            mass_matrix
        )

    end
    # print_timer()
    return nw
end

function VertexBatch(im::IndexManager, idxs::Vector{Int}; verbose)
    components = @view im.vertexf[idxs]

    try
        _compT = compT(only(unique(compT, components)))
        _compf = compf(only(unique(compf, components)))
        _statetype = statetype(only(unique(statetype, components)))
        _dim = dim(only(unique(dim, components)))
        _pdim = pdim(only(unique(pdim, components)))

        (statestride, pstride, aggbufstride) = register_vertices!(im, _statetype, _dim, _pdim, idxs)

        verbose &&
        println(" - VertexBatch: dim=$(_dim), pdim=$(_pdim), length=$(length(idxs))")

        VertexBatch{_compT, typeof(_compf)}(idxs, _compf, statestride, pstride, aggbufstride)
    catch e
        if e isa ArgumentError && startswith(e.msg, "Collection has multiple elements")
            throw(ArgumentError("Provided vertex functions $idxs use the same function but have different metadata (dim, pdim,type,...)"))
        else
            rethrow(e)
        end
    end
end

function EdgeBatch(im::IndexManager, idxs::Vector{Int}; verbose)
    components = @view im.edgef[idxs]

    try
        _compT = compT(only(unique(compT, components)))
        _compf = compf(only(unique(compf, components)))
        _statetype = statetype(only(unique(statetype, components)))
        _dim = dim(only(unique(dim, components)))
        _pdim = pdim(only(unique(pdim, components)))

        (statestride, pstride, gbufstride) = register_edges!(im, _statetype, _dim, _pdim, idxs)

        verbose &&
        println(" - EdgeBatch: dim=$(_dim), pdim=$(_pdim), length=$(length(idxs))")

        EdgeBatch{_compT, typeof(_compf)}(idxs, _compf, statestride, pstride, gbufstride)
    catch e
        if e isa ArgumentError && startswith(e.msg, "Collection has multiple elements")
            throw(ArgumentError("Provided edge functions $idxs use the same function but have different metadata (dim, pdim,type,...)"))
        else
            rethrow(e)
        end
    end
end

function NetworkLayer(im::IndexManager, batches, agg)
    map = zeros(Int, ne(im.g) * im.vdepth, 2)
    for batch in batches
        for i in 1:length(batch)
            eidx = batch.indices[i]
            e = im.edgevec[eidx]
            dst_range = im.v_data[e.src][1:im.vdepth]
            src_range = im.v_data[e.dst][1:im.vdepth]
            range = gbuf_range(batch, i)
            map[range, 1] .= dst_range
            map[range, 2] .= src_range
        end
    end
    NetworkLayer(im.g, batches, agg, im.edepth, im.vdepth, map)
end

batch_by_idxs(v, idxs::Vector{Vector{Int}}) = [v for batch in idxs]
function batch_by_idxs(v::AbstractVector, batches::Vector{Vector{Int}})
    @assert length(v) == sum(length.(batches))
    [v[batch] for batch in batches]
end

_find_identical(v::ComponentFunction, indices) = [collect(indices)]
function _find_identical(v::Vector, indices)
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

function construct_mass_matrix(im; type=nothing)
    vertexd = filter(pairs(im.vertexf)) do (_, c)
        hasproperty(c,:mass_matrix) && c.mass_matrix != LinearAlgebra.I
    end
    edged = filter(pairs(im.edgef)) do (_, c)
        hasproperty(c,:mass_matrix) && c.mass_matrix != LinearAlgebra.I
    end
    if isempty(vertexd) && isempty(edged)
        return LinearAlgebra.I
    end


    # go through all mass matrices, find type and diagonal structure of massmatrix
    all_diag = true
    _type = Bool
    for comp in Iterators.flatten((values(vertexd), values(edged)))
        _check_massmatrix(comp)
        _type = promote_type(_type, eltype(comp.mass_matrix))
        all_diag = all_diag && LinearAlgebra.isdiag(comp.mass_matrix)
    end
    if type != nothing
        _type = type
    end

    mass_matrix = if all_diag
        UniformScaling{_type}(1)(im.lastidx_dynamic)
    else
        Matrix(UniformScaling{_type}(1)(im.lastidx_dynamic))
    end
    _fill_mass_matrix!(mass_matrix, im, vertexd, edged)
end

function _fill_mass_matrix!(mass_matrix, im, vertexd, edged)
    for (i, c) in vertexd
        range = im.v_data[i]
        # cannot broadcast uniform scaling, collect to diagonal in that case
        mm = c.mass_matrix isa UniformScaling ? c.mass_matrix(length(range)) : c.mass_matrix
        mass_matrix[range, range] .= mm
    end
    for (i, c) in edged
        range = im.e_data[i]
        mm = c.mass_matrix isa UniformScaling ? c.mass_matrix(length(range)) : c.mass_matrix
        mass_matrix[range, range] .= mm
    end
    mass_matrix
end

function _check_massmatrix(c)
    if c.mass_matrix isa UniformScaling
        return
    end
    if length(size(c.mass_matrix)) == 2
        @argcheck size(c.mass_matrix) == (dim(c), dim(c))
        return
    end
    throw(ArgumentError("Mass matrix must be a square matrix,\
                         a uniform scaling, or scalar. Got $(c.mass_matrix) \
                         in component :$(c.name)."))
end

"""
    Network(nw::Network; kwargs...)

Rebuild the Network with same graph and vertex/edge functions but possibly different kwargs.
"""
function Network(nw::Network; kwargs...)
    Network(nw.im.g, nw.im.vertexf, nw.im.edgef; kwargs...)
end
