function Network(g::AbstractGraph,
                 vertexf::Union{VertexFunction,Vector{<:VertexFunction}},
                 edgef::Union{EdgeFunction,Vector{<:EdgeFunction}};
                 execution=SequentialExecution{true}(),
                 edepth=:auto,
                 vdepth=:auto,
                 aggregator=execution isa SequentialExecution ? SequentialAggregator(+) : PolyesterAggregator(+),
                 check_graphelement=true,
                 dealias=false,
                 verbose=false)
    # TimerOutputs.reset_timer!()
    @timeit_debug "Construct Network" begin
        # collect all vertex/edgf to vector
        all_same_v = vertexf isa VertexFunction
        all_same_e = edgef isa EdgeFunction
        maybecopy = dealias ? copy : identity
        _vertexf = all_same_v ? [maybecopy(vertexf) for _ in vertices(g)] : vertexf
        _edgef   = all_same_e ? [maybecopy(edgef) for _ in edges(g)] : edgef

        @argcheck _vertexf isa Vector{<:VertexFunction} "Expected VertexFuncions, got $(eltype(_vertexf))"
        @argcheck _edgef isa Vector{<:EdgeFunction} "Expected EdgeFuncions, got $(eltype(_vertexf))"
        @argcheck length(_vertexf) == nv(g)
        @argcheck length(_edgef) == ne(g)

        # search for vertex functions with feed forward
        if any(hasff, _vertexf)
            throw(ArgumentError("Vertex functions with feed forward are not supported yet! \
                As an intermediate solution, you can call `ff_to_constraint(vf)` on the vertex function\
                to turn feed forward outputs into algebraic states."))
        end

        # check if components alias eachother copy if necessary
        # allready dealiase if provided as single functions
        if dealias && !all_same_v
            dealias!(_vertexf)
        end
        if dealias && !all_same_e
            dealias!(_edgef)
        end

        verbose &&
            println("Create dynamic network with $(nv(g)) vertices and $(ne(g)) edges:")
        @argcheck execution isa ExecutionStyle "Execution type $execution not supported (choose from $(subtypes(ExecutionStyle)))"


        if !allequal(outdim, _vertexf)
            throw(ArgumentError("All vertex functions must have the same output dimension!"))
        end
        if !allequal(outdim_dst, _edgef)
            throw(ArgumentError("All edge functions must have the same output dimension!"))
        end
        vdepth = outdim(first(_vertexf))
        edepth = outdim_dst(first(_edgef))

        dynstates = mapreduce(dim, +, Iterators.flatten((_vertexf,_edgef)))

        # create index manager
        @timeit_debug "Construct Index manager" begin
            valias = dealias ? :none : (all_same_v ? :all : :some)
            ealias = dealias ? :none : (all_same_e ? :all : :some)
            im = IndexManager(g, dynstates, edepth, vdepth, _vertexf, _edgef; valias, ealias)
        end


        # check graph_element metadata and attach if necessary
        if check_graphelement
            @timeit_debug "Check graph element" begin
                if !all_same_v
                    for (i, vf) in pairs(_vertexf)
                        if has_graphelement(vf)
                            if get_graphelement(vf) != i
                                @warn "Vertex function $(vf.name) is placed at node index $i bus has \
                                `graphelement` $(get_graphelement(vf)) stored in metadata. \
                                The wrong data will be ignored! Use `check_graphelement=false` tu supress this warning."
                            end
                        end
                    end
                elseif has_graphelement(vertexf)
                    @warn "Provided vertex function has assigned `graphelement` metadata. \
                    but is used at every vertex. The `graphelement` will be ignored."
                end
                if !all_same_e
                    for (iteredge, ef) in zip(im.edgevec, _edgef)
                        if has_graphelement(ef)
                            ge = get_graphelement(ef)
                            src = get(im.unique_vnames, ge.src, ge.src)
                            dst = get(im.unique_vnames, ge.dst, ge.dst)
                            if iteredge.src != src || iteredge.dst != dst
                                @warn "Edge function $(ef.name) at $(iteredge.src) => $(iteredge.dst) has wrong `:graphelement` $src => $dst). \
                                The wrong data will be ignored! Use `check_graphelement=false` tu supress this warning."
                            end
                        end
                    end
                elseif has_graphelement(edgef)
                    @warn "Provided edge function has assigned `graphelement` metadata. \
                    but is used for all edges. The `graphelement` will be ignored."
                end
            end
        end

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

        nl = NetworkLayer(im.g, edgebatches, _aggregator, im.edepth, im.vdepth)

        @assert isdense(im)
        mass_matrix = construct_mass_matrix(im)
        N = ForwardDiff.pickchunksize(max(im.lastidx_dynamic, im.lastidx_p))
        caches = (; output = DiffCache(zeros(im.lastidx_out), N),
                  aggregation = DiffCache(zeros(im.lastidx_aggr), N))

        gbufprovider = if usebuffer(execution)
            EagerGBufProvider(im, edgebatches)
        else
            LazyGBufProvider(im, edgebatches)
        end

        nw = Network{typeof(execution),typeof(g),typeof(nl), typeof(vertexbatches),
                     typeof(mass_matrix),eltype(caches),typeof(gbufprovider)}(
            vertexbatches,
            nl, im,
            caches,
            mass_matrix,
            gbufprovider
        )

    end
    # TimerOutputs.print_timer()
    return nw
end

function Network(vertexfs, edgefs; kwargs...)
    @argcheck all(has_graphelement, vertexfs) "All vertex functions must have assigned `graphelement` to implicitly construct graph!"
    @argcheck all(has_graphelement, edgefs) "All edge functions must have assigned `graphelement` to implicitly construct graph!"

    vidxs = get_graphelement.(vertexfs)
    allunique(vidxs) || throw(ArgumentError("All vertex functions must have unique `graphelement`!"))
    sort(vidxs) == 1:length(vidxs) || throw(ArgumentError("Vertex functions must have `graphelement` in range 1:length(vertexfs)!"))

    vdict = Dict(vidxs .=> vertexfs)

    # find unique maapings from name => graphelement
    vnamedict = unique_mappings(getproperty.(vertexfs, :name), get_graphelement.(vertexfs))

    simpleedges = map(edgefs) do e
        ge = get_graphelement(e)
        src = get(vnamedict, ge.src, ge.src)
        dst = get(vnamedict, ge.dst, ge.dst)
        if src isa Symbol || dst isa Symbol
            throw(ArgumentError("Edge graphelement $src => $dst continas non-unique or unknown vertex names!"))
        end
        SimpleEdge(src, dst)
    end
    allunique(simpleedges) || throw(ArgumentError("Some edge functions have the same `graphelement`!"))
    edict = Dict(simpleedges .=> edgefs)

    # if all src < dst then we can use SimpleGraph, else digraph
    g = if all(e -> e.src < e.dst, simpleedges)
        SimpleGraph(length(vertexfs))
    else
        SimpleDiGraph(length(vertexfs))
    end
    for edge in simpleedges
        if g isa SimpleDiGraph && has_edge(g, edge.dst, edge.src)
            @warn "Edges $(edge.src) -> $(edge.dst) and $(edge.dst) -> $(edge.src) are both present in the graph!"
        end
        r = add_edge!(g, edge)
        r || error("Could not add edge $(edge) to graph $(g)!")
    end

    vfs_ordered = [vdict[k] for k in vertices(g)]
    efs_ordered = [edict[k] for k in edges(g)]

    Network(g, vfs_ordered, efs_ordered; check_graphelement=false, kwargs...)
end

"""
    dealias!(cfs::Vector{<:ComponentFunction})

Checks if any component functions reference the same metadtata/symmetada fields and
creates copies of them if necessary.
"""
function dealias!(cfs::Vector{<:ComponentFunction}; warn=false)
    ag = aliasgroups(cfs)

    isempty(ag) && return cfs # nothign to do

    copyidxs = reduce(vcat, values(ag))

    cfs[copyidxs] = copy.(cfs[copyidxs])
end


"""
    aliasgroups(cfs::Vector{<:ComponentFunction})

Returns a dict `cf => idxs` which contains all the component functions
which appear multiple times with all their indices.
"""
function aliasgroups(cfs::Vector{T}) where {T<:ComponentFunction}
    d = IdDict{T, Vector{Int}}()

    for (i, cf) in pairs(cfs)
        if haskey(d, cf) # c allready present
            push!(d[cf], i)
        else
            d[cf] = [i]
        end
    end

    filter!(x -> length(x.second) > 1, d)
end

function VertexBatch(im::IndexManager, idxs::Vector{Int}; verbose)
    components = @view im.vertexf[idxs]

    try
        # TODO: those checks seems expensive and redundant
        _compT = dispatchT(only(unique(dispatchT, components)))
        _compf = compf(only(unique(compf, components)))
        _compg = compg(only(unique(compg, components)))
        _ff    = fftype(only(unique(fftype, components)))

        _dim = dim(only(unique(dim, components)))
        _outdim = outdim(only(unique(outdim, components)))
        _pdim = pdim(only(unique(pdim, components)))

        (statestride, outstride, pstride, aggbufstride) =
            register_vertices!(im, _dim, _outdim, _pdim, idxs)

        verbose &&
        println(" - VertexBatch: dim=$(_dim), pdim=$(_pdim), length=$(length(idxs))")

        VertexBatch{_compT, typeof(_compf), typeof(_compg), typeof(_ff), typeof(idxs)}(
            idxs, _compf, _compg, _ff, statestride, outstride, pstride, aggbufstride)
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
        # TODO: those checks seems expensive and redundant
        _compT = dispatchT(only(unique(dispatchT, components)))
        _compf = compf(only(unique(compf, components)))
        _compg = compg(only(unique(compg, components)))
        _ff    = fftype(only(unique(fftype, components)))

        _dim = dim(only(unique(dim, components)))
        _outdim = outdim(only(unique(outdim, components)))
        _pdim = pdim(only(unique(pdim, components)))

        (statestride, outstride, pstride, gbufstride) =
            register_edges!(im, _dim, _outdim, _pdim, idxs)

        verbose &&
        println(" - EdgeBatch: dim=$(_dim), pdim=$(_pdim), length=$(length(idxs))")

        EdgeBatch{_compT, typeof(_compf), typeof(_compg), typeof(_ff), typeof(idxs)}(
            idxs, _compf, _compg, _ff, statestride, outstride, pstride, gbufstride)
    catch e
        if e isa ArgumentError && startswith(e.msg, "Collection has multiple elements")
            throw(ArgumentError("Provided edge functions $idxs use the same function but have different metadata (dim, pdim,type,...)"))
        else
            rethrow(e)
        end
    end
end

batch_by_idxs(v, idxs::Vector{Vector{Int}}) = [v for batch in idxs]
function batch_by_idxs(v::AbstractVector, batches::Vector{Vector{Int}})
    @assert length(v) == sum(length.(batches))
    [v[batch] for batch in batches]
end

_find_identical(v::ComponentFunction, indices) = [collect(indices)]
function _find_identical(v::Vector, indices)
    idxs_per_type = Vector{Int}[]
    unique_comp = []
    for i in eachindex(v)
        found = false
        for j in eachindex(unique_comp)
            if batchequal(v[i], unique_comp[j])
                found = true
                push!(idxs_per_type[j], indices[i])
                break
            end
        end
        if !found
            push!(unique_comp, v[i])
            push!(idxs_per_type, [indices[i]])
        end
    end
    @assert length(unique_comp) == length(idxs_per_type)
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
    Network(nw::Network; g, vertexf, edgef, kwargs...)

Rebuild the Network with same graph and vertex/edge functions but possibly different kwargs.
"""
function Network(nw::Network;
                 g = nw.im.g,
                 vertexf = copy.(nw.im.vertexf),
                 edgef = copy.(nw.im.edgef),
                 kwargs...)

    _kwargs = Dict(:execution => executionstyle(nw),
                   :edepth => :auto,
                   :vdepth => :auto,
                   :aggregator => get_aggr_constructor(nw.layer.aggregator),
                   :check_graphelement => true,
                   :dealias => false,
                   :verbose => false)
    for (k, v) in kwargs
        _kwargs[k] = v
    end

    # check, that we actually provide all of the arguments
    # mainly so we don't forget to add it here if we introduce new kw arg to main constructor
    m = only(methods(Network, [typeof(g), typeof(vertexf), typeof(edgef)]))
    @assert keys(_kwargs) == Set(Base.kwarg_decl(m))

    Network(g, vertexf, edgef; _kwargs...)
end
