module NetworkDynamics

using DiffEqBase
using Graphs

include("Utilities.jl")
include("ComponentFunctions.jl")
include("NetworkStructures.jl")
include("NetworkDiffEq.jl")


export network_dynamics



"""
This function assembles the arrays that hold the structural information
of the individual vertices and edges: Their dimensions, the symbols of their
variables and their mass matrices. Used internally by the network_dynamics
constructors.
"""
function collect_ve_info(vertices!, edges!, graph)
    if vertices! isa Array
        @assert length(vertices!) == nv(graph)
        v_dims = Int[v.dim for v in vertices!]
        symbols_v = [Symbol(vertices![i].sym[j],"_",i)
                           for i in 1:length(vertices!)
                           for j in 1:v_dims[i]]
        mmv_array = [v.mass_matrix for v in vertices!]
    else
        v_dims = [vertices!.dim for v in vertices(graph)]
        symbols_v = [Symbol(vertices!.sym[j],"_",i)
                     for i in 1:nv(graph)
                     for j in 1:v_dims[i]]
        mmv_array = [vertices!.mass_matrix for v in vertices(graph)]
    end

    if edges! isa Array
        @assert length(edges!) == ne(graph)
        e_dims = Int[e.dim for e in edges!]
        symbols_e = [Symbol(edges![i].sym[j],"_",i)
                           for i in 1:length(edges!)
                           for j in 1:e_dims[i]]
        if eltype(edges!)  <: Union{StaticEdge, StaticDelayEdge}  # improve type hierarchy
            mme_array = nothing
        else
            mme_array = [e.mass_matrix for e in edges!]
        end
    else
        e_dims = [edges!.dim for e in edges(graph)]
        symbols_e = [Symbol(edges!.sym[j],"_",i)
                     for i in 1:ne(graph)
                     for j in 1:e_dims[i]]
        if typeof(edges!) <: Union{StaticEdge, StaticDelayEdge} # improve type hierarchy
            mme_array = nothing
        else
            mme_array = [edges!.mass_matrix for e in edges(graph)]
        end
    end

    v_dims, e_dims, symbols_v, symbols_e, mmv_array, mme_array
end


function collect_unique_components(comp::T, NV) where T <: AbstractArray
    ucomp = Vector(unique(comp))
    ind_ucomp = Vector{Int64}[]
    for v in ucomp
        ind_v = Int64[]
        for i in eachindex(comp)
            if v == comp[i]
                push!(ind_v, i)
            end
        end
        push!(ind_ucomp, ind_v)
    end
    ucomp, ind_ucomp
end

function collect_unique_components(comp, NV)
    ucomp = [comp]
    ind_ucomp = [collect(1:NV)]
    ucomp, ind_ucomp
end


## This is the "startpage" version of ND that checks for incompatibile combinations

"""
    network_dynamics(vertices!, edges!, g; parallel = false)

Assembles the the dynamical equations of the network problem into an `ODEFunction`
compatible with the `DifferentialEquations.jl` solvers. Takes as arguments an array
of VertexFunctions **`vertices!`**, an array of EdgeFunctions **`edges!`** and a
`Graphs.jl` object **`g`**. The optional argument `parallel` is a boolean
value that denotes if the central loop should be executed in parallel with the number of threads set by the environment variable `JULIA_NUM_THREADS`.
"""
function network_dynamics(vertices!::Vector{T}, edges!::Vector{U}, graph; kwargs...) where {T <: VertexFunction, U <: EdgeFunction}
    @assert length(vertices!) == nv(graph)
    @assert length(edges!) == ne(graph)


    hasDelay = any(v -> v isa DDEVertex, vertices!) ||
               any(e -> e isa StaticDelayEdge, edges!)
    hasODEEdge = any(e -> e isa ODEEdge, edges!)


    hasDelay && hasODEEdge ? error(
        ArgumentError("ODEEdges with delay are not supported at the moment.")) : nothing

    # If one edge is an ODEEdge all other edges will be promoted. Eventually we will get rid of promotions.
    if hasODEEdge && any(e -> e isa StaticEdge, edges!)
        edges! = Array{ODEEdge}(edges!)
    end

    return _network_dynamics(vertices!, edges!, graph; kwargs...)
end


function network_dynamics(vertices!,  edges!, graph; parallel=false)
    # If vertices! and/or edges! are individual functions and no other dispatch was
    # triggered, assume all vertices, respectively edges will be of that type
    if typeof(vertices!) <: VertexFunction
        vertices! = [vertices! for i in 1:nv(graph)]
    end
    if typeof(edges!) <: EdgeFunction
        edges! = [edges! for i in 1:ne(graph)]
    end

    try
        Array{VertexFunction}(vertices!)
    catch err
        throw(ArgumentError("Cannot convert the vertices to an Array{VertexFunction}!"))
    end

    try
        Array{EdgeFunction}(edges!)
    catch err
        throw(ArgumentError("Cannot convert the edges to an Array{EdgeFunction}!"))
    end

    vertices! isa Array{Any} ? vertices! = Array{VertexFunction}(vertices!) : nothing
    edges! isa Array{Any} ? edges! = Array{EdgeFunction}(edges!) : nothing
    network_dynamics(vertices!, edges!, graph, parallel = parallel)
end


# catch all version (works only for ODEVertex, DDEVertex, StaticEdge, StaticDelayEdge)
function _network_dynamics(vertices!::Union{Array{T, 1}, T},
                          edges!::Union{Array{U, 1}, U},
                          graph;
                          initial_history=nothing,
                          x_prototype=zeros(1),
                          parallel=false) where {T <: VertexFunction, U <: EdgeFunction}


    warn_parallel(parallel)

    edges! = prepare_edges(edges!, graph)

    v_dims, e_dims, symbols_v, symbols_e, mmv_array, mme_array = collect_ve_info(vertices!, edges!, graph)

    # These arrays are used for initializing the GraphData and will be overwritten
    v_array = similar(x_prototype, sum(v_dims))
    e_array = similar(x_prototype, sum(e_dims))


    unique_vertices!, unique_v_indices = collect_unique_components(vertices!, nv(graph))
    unique_edges!, unique_e_indices = collect_unique_components(edges!, ne(graph))

    hasDelay = any(v -> v isa DDEVertex, unique_vertices!) ||
               any(e -> e isa StaticDelayEdge, unique_edges!)

    if initial_history === nothing && hasDelay
      initial_history = ones(sum(v_dims))
    end

    graph_stucture = GraphStruct(graph, v_dims, e_dims, symbols_v, symbols_e)

    graph_data = GraphData(v_array, e_array, graph_stucture)

    nd! = NetworkDE(unique_vertices!, unique_v_indices, unique_edges!, unique_e_indices, graph, graph_stucture, graph_data, initial_history, parallel)
    mass_matrix = construct_mass_matrix(mmv_array, graph_stucture)

    return hasDelay ? DDEFunction(nd!; mass_matrix = mass_matrix, syms=symbols_v) :
                      ODEFunction(nd!; mass_matrix = mass_matrix, syms=symbols_v)
end

## ODEEdge

function _network_dynamics(vertices!::Union{Vector{T}, T}, edges!::Union{Vector{U}, U}, graph; initial_history=nothing, x_prototype=zeros(1), parallel=false) where {T <: ODEVertex, U <: ODEEdge}

    warn_parallel(parallel)

    # We are abusing prepare_edges here to throw errors for impossible connections of
    # coupling type and graph. We don't save it's return value though, since reconstruction
    # of ODEEdges is forbidden at the moment.

    prepare_edges(edges!, graph)

    v_dims, e_dims, symbols_v, symbols_e, mmv_array, mme_array = collect_ve_info(vertices!, edges!, graph)

    # These arrays are used for initializing the GraphData and will be overwritten
    x_array = similar(x_prototype, sum(v_dims) + sum(e_dims))
    v_array = view(x_array, 1:sum(v_dims))
    e_array = view(x_array, sum(v_dims)+1:sum(v_dims)+sum(e_dims))

    symbols = vcat(symbols_v, symbols_e)

    unique_vertices!, unique_v_indices = collect_unique_components(vertices!, nv(graph))
    unique_edges!, unique_e_indices = collect_unique_components(edges!, ne(graph))

    graph_stucture = GraphStruct(graph, v_dims, e_dims, symbols_v, symbols_e)

    graph_data = GraphData(v_array, e_array, graph_stucture)

    nd! = NetworkDE(unique_vertices!, unique_v_indices, unique_edges!, unique_e_indices, graph, graph_stucture, graph_data, initial_history, parallel)

    mass_matrix = construct_mass_matrix(mmv_array, mme_array, graph_stucture)

    ODEFunction(nd!; mass_matrix = mass_matrix, syms=symbols)
end

# Here is not the best place for this error,
# but requires the least changes before the next refactor.
function prepare_edges(edges, g)
    @assert typeof(edges) <: Union{EdgeFunction, Vector}
    throw(ArgumentError("Graph type not recognized. Currently only SimpleGraph and SimpleDiGraph are supported."))
end
"""
If only a sinlge Function is given, not an Array of EdgeFunctions.
"""
function prepare_edges(edge::EdgeFunction, g::SimpleGraph)
    if edge.coupling == :directed
        throw(ArgumentError("Coupling type of EdgeFunction not available for undirected Graphs"))
    elseif edge.coupling == :undefined
        @info("Reconstructing EdgeFunction with :undefined coupling type...")
        return reconstruct_edge(edge, :undirected)
    end
    return edge
end

function prepare_edges(edge::EdgeFunction, g::SimpleDiGraph)
    if edge.coupling ∈ (:symmetric, :antisymmetric, :undirected, :fiducial)
        throw(ArgumentError("Coupling type of EdgeFunction not available for directed Graphs"))
    elseif edge.coupling == :undefined
        @info("Reconstructing EdgeFunction with :undefined coupling type...")
        return reconstruct_edge(edge, :directed)
    end
    return edge
end

""" prepare_edges(edges, g::SimpleGraph)


"""
function prepare_edges(edges::Vector, g::SimpleGraph)
    # Depending on the coupling we might get different eltypes
    new_edges = Vector{EdgeFunction}(undef, length(edges))
    infobool = true
    for (i, edge) in enumerate(edges)
        if edge.coupling == :directed
            throw(ArgumentError("Coupling type of edge $i not available for undirected Graphs"))
        elseif edge.coupling == :undefined
            if infobool
                @info("Reconstructing EdgeFuntions with :undefined coupling type.")
                infobool = false
            end
            new_edges[i] = reconstruct_edge(edge, :undirected)
        else
            new_edges[i] = edges[i]
        end
    end
    # by recreating the array the eltype gets narrowed down as much as possible
    narrowed_type = [e for e in new_edges]
    return narrowed_type
end

"""
"""
function prepare_edges(edges::Vector, g::SimpleDiGraph)
    # Depending on the coupling we might get different eltypes
    new_edges = Vector{EdgeFunction}(undef, length(edges))
    infobool = true
    for (i, edge) in enumerate(edges)
        if edge.coupling ∈ (:symmetric, :antisymmetric, :undirected, :fiducial)
            throw(ArgumentError("Coupling type of edge $i not available for directed Graphs"))
        elseif edge.coupling == :undefined
            if infobool
                @info("Reconstructing EdgeFuntions with :undefined coupling type.")
                infobool = false
            end
            new_edges[i] = reconstruct_edge(edge, :directed)
        else
            new_edges[i] = edges[i]
        end
    end
    # by recreating the array the eltype gets narrowed down as much as possible
    narrowed_type = [e for e in new_edges]
    return narrowed_type
end



@inline function reconstruct_edge(edge::StaticEdge, coupling::Symbol)
    let f! = edge.f!, dim = edge.dim, sym = edge.sym
        return StaticEdge(f! = f!,
                          dim = dim,
                          coupling = coupling,
                          sym = sym)
    end
end
@inline function reconstruct_edge(edge::StaticDelayEdge, coupling::Symbol)
    let f! = edge.f!, dim = edge.dim, sym = edge.sym
        return StaticDelayEdge(f! = f!,
                               dim = dim,
                               coupling = coupling,
                               sym = sym)
    end
end
@inline function reconstruct_edge(edge::ODEEdge, coupling::Symbol)
    error("Reconstruction of ODEEdges is not implemented at the moment.")
end

end # module
