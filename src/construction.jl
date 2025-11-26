"""
    Network([g,] vertexf, edgef; kwarg...)

Construct a `Network` object from a graph `g` and edge and component models `vertexf` and `edgef`.

Arguments:
 - `g::AbstractGraph`: The graph on which the network is defined.
    Optional, can be omitted if all component models have a defined `graphelement`.
    See `vidx` and `src`/`dst` keywors for [`VertexModel`](@ref) and [`EdgeModel`](@ref) constructors respectively.
 - `vertexm`:
    A single [`VertexModel`](@ref) or a vector of [`VertexModel`](@ref) objects.
    The order of the vertex models must mirror the order of the `vertices(g)` iterator.

 - `edgem`: A single [`EdgeModel`](@ref) or a vector of [`EdgeModel`](@ref) objects.
    The order of the edge models must mirror the order of the `edges(g)` iterator.

Optional keyword arguments:
 - `execution=SequentialExecution{true}()`:
    Execution model of the network. E.g. [`SequentialExecution`](@ref), [`KAExecution`](@ref), [`PolyesterExecution`](@ref) or [`ThreadedExecution`](@ref).
 - `aggregator=execution is a SequentialExecution ? SequentialAggregator(+) : PolyesterAggregator(+)`: (@Hans : is this correct?)
    Aggregation function applied to the edge models. E.g. [`SequentialAggregator`](@ref), [`PolyesterAggregator`](@ref), [`ThreadedAggregator`](@ref), [`SparseAggregator`](@ref).
 - `check_graphelement=true`:
    Check if the `graphelement` metadata is consistent with the graph.
 - `dealias=false`
    Check if the components alias each other and create copies if necessary.
    This is necessary if the same component model is referenced in multiple places in the Network but you want to
    dynamically assign metadata, such as initialization information to specific instances.
 - `verbose=false`:
    Show additional information during construction.
"""
function Network(g::AbstractGraph,
                 vertexm::Union{VertexModel,Vector{<:VertexModel}},
                 edgem::Union{EdgeModel,Vector{<:EdgeModel}};
                 execution=SequentialExecution{true}(),
                 aggregator=execution isa SequentialExecution ? SequentialAggregator(+) : PolyesterAggregator(+),
                 check_graphelement=true,
                 dealias=false,
                 verbose=false)
    # TimerOutputs.reset_timer!()
    @timeit_debug "Construct Network" begin
        # collect all vertex/edgf to vector
        all_same_v = vertexm isa VertexModel
        all_same_e = edgem isa EdgeModel
        maybecopy = dealias ? copy : identity
        _vertexm = all_same_v ? [maybecopy(vertexm) for _ in vertices(g)] : vertexm
        _edgem   = all_same_e ? [maybecopy(edgem) for _ in edges(g)] : edgem

        @argcheck _vertexm isa Vector{<:VertexModel} "Expected VertexModels, got $(eltype(_vertexm))"
        @argcheck _edgem isa Vector{<:EdgeModel} "Expected EdgeModels, got $(eltype(_vertexm))"
        @argcheck length(_vertexm) == nv(g)
        @argcheck length(_edgem) == ne(g)

        # search for vertex models with feed forward
        if CHECK_COMPONENT[] && any(hasff, _vertexm)
            throw(ArgumentError("Vertex model with feed forward are not supported yet! \
                As an intermediate solution, you can call `ff_to_constraint(vf)` on the vertex model\
                to turn feed forward outputs into algebraic states."))
        end

        # check if components alias each other copy if necessary
        # already dealiase if provided as single functions
        if dealias && !all_same_v
            dealias!(_vertexm)
        end
        if dealias && !all_same_e
            dealias!(_edgem)
        end

        verbose &&
            println("Create dynamic network with $(nv(g)) vertices and $(ne(g)) edges:")
        @argcheck execution isa ExecutionStyle "Execution type $execution not supported (choose from $(subtypes(ExecutionStyle)))"


        if !allequal(outdim, _vertexm)
            throw(ArgumentError("All vertex models must have the same output dimension!"))
        end
        if !allequal(outdim_dst, _edgem)
            throw(ArgumentError("All edge models must have the same output dimension!"))
        end
        vdepth = outdim(first(_vertexm))
        edepth = outdim_dst(first(_edgem))

        dynstates = mapreduce(dim, +, Iterators.flatten((_vertexm,_edgem)))

        # create index manager
        @timeit_debug "Construct Index manager" begin
            valias = dealias ? :none : (all_same_v ? :all : :some)
            ealias = dealias ? :none : (all_same_e ? :all : :some)
            im = IndexManager(g, dynstates, edepth, vdepth, _vertexm, _edgem; valias, ealias)
        end


        # check graph_element metadata and attach if necessary
        if check_graphelement
            @timeit_debug "Check graph element" begin
                if !all_same_v
                    for (i, vf) in pairs(_vertexm)
                        if has_graphelement(vf)
                            if get_graphelement(vf) != i
                                @warn "Vertex model $(vf.name) is placed at node index $i bus has \
                                `graphelement` $(get_graphelement(vf)) stored in metadata. \
                                The wrong data will be ignored! Use `check_graphelement=false` tu supress this warning."
                            end
                        end
                    end
                elseif has_graphelement(vertexm)
                    @warn "Provided vertex model has assigned `graphelement` metadata. \
                    but is used at every vertex. The `graphelement` will be ignored."
                end
                if !all_same_e
                    for (iteredge, ef) in zip(im.edgevec, _edgem)
                        if has_graphelement(ef)
                            ge = get_graphelement(ef)
                            src = get(im.unique_vnames, ge.src, ge.src)
                            dst = get(im.unique_vnames, ge.dst, ge.dst)
                            if iteredge.src != src || iteredge.dst != dst
                                @warn "Edge model $(ef.name) at $(iteredge.src) => $(iteredge.dst) has wrong `:graphelement` $src => $dst). \
                                The wrong data will be ignored! Use `check_graphelement=false` tu supress this warning."
                            end
                        end
                    end
                elseif has_graphelement(edgem)
                    @warn "Provided edge model has assigned `graphelement` metadata. \
                    but is used for all edges. The `graphelement` will be ignored."
                end
            end
        end

        # batch identical edge and vertex model
        @timeit_debug "batch identical vertexes" begin
            vidxs = if all_same_v
                [collect(1:nv(g))]
            else
                _find_identical_components(_vertexm)
            end
        end
        @timeit_debug "batch identical edges" begin
            eidxs = if all_same_e
                [collect(1:ne(g))]
            else
                _find_identical_components(_edgem)
            end
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
                    aggregation = DiffCache(zeros(im.lastidx_aggr), N),
                    external = DiffCache(zeros(im.lastidx_extbuf), N))

        gbufprovider = if usebuffer(execution)
            EagerGBufProvider(im, edgebatches)
        else
            LazyGBufProvider(im, edgebatches)
        end

        # create map for extenral inputs
        extmap = has_external_input(im) ? ExtMap(im) : nothing

        nw = Network{typeof(execution),typeof(g),typeof(nl), typeof(vertexbatches),
                     typeof(mass_matrix),eltype(caches),typeof(gbufprovider),typeof(extmap)}(
            vertexbatches,
            nl, im,
            caches,
            mass_matrix,
            gbufprovider,
            extmap,
            Ref{Union{Nothing,SparseMatrixCSC{Bool,Int}}}(nothing),
        )

    end
    # TimerOutputs.print_timer()
    return nw
end

function _find_identical_components(models)
    # identical components are based on identical _component_hash
    # those can have different metadata but are considered identical when in comes to batching
    hashs = _component_hash.(models)
    find_identical(hashs)
end
# hash condition: components with same hash will end up in the same batch
function _component_hash(c::ComponentModel)
    hash((
        typeof(c), # same type
        compf(c),  # same f-function
        compg(c),  # same g-function
        fftype(c), # same feedforward type
        dim(c),    # same state dimension
        outdim(c), # same output dimension
        pdim(c),   # same parameter dimension
        extdim(c), # same external input dimension
    ))
end

function Network(vertexfs, edgefs; warn_order=true, kwargs...)
    vertexfs = vertexfs isa VertexModel ? [vertexfs] : vertexfs
    edgefs   = edgefs isa EdgeModel     ? [edgefs]   : edgefs
    @argcheck all(has_graphelement, vertexfs) "All vertex models must have assigned `graphelement` to implicitly construct graph!"
    @argcheck all(has_graphelement, edgefs) "All edge models must have assigned `graphelement` to implicitly construct graph!"

    vidxs = get_graphelement.(vertexfs)
    allunique(vidxs) || throw(ArgumentError("All vertex models must have unique `graphelement`!"))
    sort(vidxs) == 1:length(vidxs) || throw(ArgumentError("Vertex models must have `graphelement` in range 1:length(vertexfs)!"))

    vdict = Dict(vidxs .=> vertexfs)

    # find unique mappings from name => graphelement
    vnamedict = unique_mappings(getproperty.(vertexfs, :name), get_graphelement.(vertexfs))

    simpleedges = map(edgefs) do e
        ge = get_graphelement(e)
        src = get(vnamedict, ge.src, ge.src)
        dst = get(vnamedict, ge.dst, ge.dst)
        if src isa Symbol || dst isa Symbol
            throw(ArgumentError("Edge graphelement $src => $dst contains non-unique or unknown vertex names!"))
        end
        SimpleEdge(src, dst)
    end
    if !allunique(simpleedges)
        msg = "Some edge models have the same `graphelement`!"
        alldup = filter(x -> length(x) > 1, find_identical(simpleedges))
        for dup in alldup
            msg *= "\n - Edges $dup refer to same element $(simpleedges[first(dup)])"
        end
        throw(ArgumentError(msg))
    end
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
    if warn_order && any(vfs_ordered .!== vertexfs)
        @warn "Order of vertex models was changed to match the natural ordering of vertices in graph (as in `vertices(g)`)! \
               Concretely, this means that `nw[VIndex(1)]` referencs the vertex with `vidxs=1` \
               **and not necessarily the first vertex model in the provided list**. \
               Disable warning with kw `warn_order=false`."
    end
    if warn_order && any(efs_ordered .!== edgefs)
        @warn "Order of edge models was changed to match the natural ordering of edges in graph (as in `edges(g)`)! \
               Concretely, this means that `nw[EIndex(1)]` references the first edge from `edges(g)` \
               **and not necessarily the first edge model in the provided list**. \
               Disable warning with kw `warn_order=false`."
    end

    Network(g, vfs_ordered, efs_ordered; check_graphelement=false, kwargs...)
end

"""
    dealias!(cfs::Vector{<:ComponentModel})

Checks if any component models reference the same metadtata/symmetada fields and
creates copies of them if necessary.
"""
function dealias!(cfs::Vector{<:ComponentModel}; warn=false)
    ag = aliasgroups(cfs)

    isempty(ag) && return cfs # nothign to do

    copyidxs = reduce(vcat, values(ag))

    cfs[copyidxs] = copy.(cfs[copyidxs])
end


"""
    aliasgroups(cfs::Vector{<:ComponentModel})

Returns a dict `cf => idxs` which contains all the component models
which appear multiple times with all their indices.
"""
function aliasgroups(cfs::Vector{T}) where {T<:ComponentModel}
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
    first_comp = im.vertexm[first(idxs)]
    _compT  = dispatchT(first_comp)
    _compf  = compf(first_comp)
    _compg  = compg(first_comp)
    _ff     = fftype(first_comp)
    _dim    = dim(first_comp)
    _pdim   = pdim(first_comp)
    _outdim = outdim(first_comp)
    _extdim = extdim(first_comp)

    strides = register_vertices!(im, idxs, _dim, _outdim, _pdim, _extdim)

    verbose && println(" - VertexBatch: dim=$(_dim), pdim=$(_pdim), length=$(length(idxs))")

    ComponentBatch(_compT, idxs, _compf, _compg, _ff, strides.state, strides.p, strides.in, strides.out, strides.ext)
end

function EdgeBatch(im::IndexManager, idxs::Vector{Int}; verbose)
    first_comp = im.edgem[first(idxs)]
    _compT  = dispatchT(first_comp)
    _compf  = compf(first_comp)
    _compg  = compg(first_comp)
    _ff     = fftype(first_comp)
    _dim    = dim(first_comp)
    _pdim   = pdim(first_comp)
    _outdim = outdim(first_comp)
    _extdim = extdim(first_comp)

    strides = register_edges!(im, idxs, _dim, _outdim, _pdim, _extdim)

    verbose && println(" - EdgeBatch: dim=$(_dim), pdim=$(_pdim), length=$(length(idxs))")

    ComponentBatch(_compT, idxs, _compf, _compg, _ff, strides.state, strides.p, strides.in, strides.out, strides.ext)
end

batch_by_idxs(v, idxs::Vector{Vector{Int}}) = [v for batch in idxs]
function batch_by_idxs(v::AbstractVector, batches::Vector{Vector{Int}})
    @assert length(v) == sum(length.(batches))
    [v[batch] for batch in batches]
end

function construct_mass_matrix(im; type=nothing)
    vertexd = filter(pairs(im.vertexm)) do (_, c)
        hasproperty(c,:mass_matrix) && c.mass_matrix != LinearAlgebra.I
    end
    edged = filter(pairs(im.edgem)) do (_, c)
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
    Network(nw::Network; g, vertexm, edgem, kwargs...)

Rebuild the Network with same graph and vertex/edge models but possibly different kwargs.
"""
function Network(nw::Network;
                 g = nw.im.g,
                 vertexm = copy.(nw.im.vertexm),
                 edgem = copy.(nw.im.edgem),
                 kwargs...)

    _kwargs = Dict(:execution => executionstyle(nw),
                   :aggregator => get_aggr_constructor(nw.layer.aggregator),
                   :check_graphelement => true,
                   :dealias => false,
                   :verbose => false)
    for (k, v) in kwargs
        _kwargs[k] = v
    end

    # check, that we actually provide all of the arguments
    # mainly so we don't forget to add it here if we introduce new kw arg to main constructor
    m = only(methods(Network, [typeof(g), typeof(vertexm), typeof(edgem)]))
    wrong_kwargs = setdiff(keys(_kwargs), Set(Base.kwarg_decl(m)))
    if !isempty(wrong_kwargs)
        throw(ArgumentError("Got unknown keyword arguments $(collect(wrong_kwargs)) for copy-constructor Network(nw; kwargs...).\n\
                             Possible arguments: $(collect(Set(Base.kwarg_decl(m))))"))
    end
    @assert keys(_kwargs) == Set(Base.kwarg_decl(m))

    Network(g, vertexm, edgem; _kwargs...)
end
