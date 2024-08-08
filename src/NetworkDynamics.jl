module NetworkDynamics

using SciMLBase
using Graphs
import Graphs: AbstractSimpleGraph

include("Utilities.jl")
include("ComponentFunctions.jl")
include("NetworkStructures.jl")
include("NetworkDiffEq.jl")


export network_dynamics, allocation_report



"""
    collect_ve_info(vertices!, edges!, graph)

This function assembles the arrays that hold the structural information
of the individual vertices and edges: Their dimensions, the symbols of their
variables and their mass matrices. Used internally by the network_dynamics
constructors.
"""
function collect_ve_info(vertices!, edges!, graph)
    if vertices! isa Array
        @assert length(vertices!) == nv(graph)
        v_dims = Int[v.dim for v in vertices!]
        symbols_v = [Symbol(vertices![i].sym[j], "_", i)
                     for i in 1:length(vertices!)
                     for j in 1:v_dims[i]]
        mmv_array = [v.mass_matrix for v in vertices!]
    else
        v_dims = [vertices!.dim for v in vertices(graph)]
        symbols_v = [Symbol(vertices!.sym[j], "_", i)
                     for i in 1:nv(graph)
                     for j in 1:v_dims[i]]
        mmv_array = [vertices!.mass_matrix for v in vertices(graph)]
    end

    if edges! isa Array
        @assert length(edges!) == ne(graph)
        e_dims = Int[e.dim for e in edges!]
        symbols_e = [Symbol(edges![i].sym[j], "_", i)
                     for i in 1:length(edges!)
                     for j in 1:e_dims[i]]
        hasMM = any(e -> hasproperty(e, :mass_matrix), edges!) # ODEEdges are the only ones with MM atm
        if hasMM # 1 ODEEdge should imply all e are ODEEdge
            mme_array = [e.mass_matrix for e in edges!]
        else
            mme_array = nothing
        end
    else
        e_dims = [edges!.dim for e in edges(graph)]
        symbols_e = [Symbol(edges!.sym[j], "_", i)
                     for i in 1:ne(graph)
                     for j in 1:e_dims[i]]
        if hasproperty(edges!, :mass_matrix)
            mme_array = [edges!.mass_matrix for e in edges(graph)]
        else
            mme_array = nothing
        end
    end

    v_dims, e_dims, symbols_v, symbols_e, mmv_array, mme_array
end


function collect_unique_components(comp::T, NV) where {T<:AbstractArray}
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
function network_dynamics(vertices!::Union{T,Vector{T}},
                          edges!::Union{U,Vector{U}},
                          graph::AbstractSimpleGraph;
                          kwargs...) where {T,U}
    if vertices! isa Vector
        @assert length(vertices!) == nv(graph)
        if length(vertices!) == 0
            # empty graph
            vertices! = VertexFunction[]
        end
        if !(eltype(vertices!)<:VertexFunction)
            # narrow type
            vertices! = [v for v in vertices!]
        end
        hasDelayVertex = any(v -> v isa DDEVertex, vertices!)
    else
        hasDelayVertex = vertices! isa DDEVertex
    end

    if edges! isa Vector
        @assert length(edges!) == ne(graph)
        if length(edges!) == 0
            # graph without lines
            edges! = EdgeFunction[]
        end
        if !(eltype(edges!)<:EdgeFunction)
            # narrow type
            edges! = [e for e in edges!]
        end
        hasDelayEdge = any(e -> e isa StaticDelayEdge, edges!)
        hasODEEdge = any(e -> e isa ODEEdge, edges!)
        hasStaticEdge = any(e -> e isa StaticEdge, edges!)
    else
        hasDelayEdge = edges! isa StaticDelayEdge
        hasODEEdge = edges! isa ODEEdge
        hasStaticEdge = edges! isa StaticEdge
    end

    hasDelay = hasDelayVertex || hasDelayEdge

    hasDelay && hasODEEdge && error(ArgumentError("ODEEdges with delay are not supported at the moment."))

    # If one edge is an ODEEdge all other edges will be promoted. Eventually we will get rid of promotions.
    if hasODEEdge && hasStaticEdge
        edges! = Vector{ODEEdge}(edges!)
    end

    return _network_dynamics(vertices!, edges!, graph; kwargs...)
end

# catch all version (works only for ODEVertex, DDEVertex, StaticEdge, StaticDelayEdge)
function _network_dynamics(vertices!::Union{Vector{T},T},
                           edges!::Union{Vector{U},U},
                           graph;
                           initial_history=nothing,
                           x_prototype=zeros(1),
                           parallel=false) where {T<:VertexFunction,U<:EdgeFunction}


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

    graph_stucture = GraphStruct(graph, v_dims, e_dims, symbols_v, symbols_e)

    graph_data = GraphData(v_array, e_array, graph_stucture)

    nd! = NetworkDE(unique_vertices!, unique_v_indices, unique_edges!, unique_e_indices, graph, graph_stucture, graph_data, parallel)
    mass_matrix = construct_mass_matrix(mmv_array, graph_stucture)

    return if hasDelay
        DDEFunction(nd!; mass_matrix=mass_matrix, syms=symbols_v)
    else
        ODEFunction(nd!; mass_matrix=mass_matrix, syms=symbols_v)
    end
end

## ODEEdge

function _network_dynamics(vertices!::Union{Vector{T},T}, edges!::Union{Vector{U},U}, graph; x_prototype=zeros(1), parallel=false) where {T<:ODEVertex, U<:ODEEdge}

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

    nd! = NetworkDE(unique_vertices!, unique_v_indices, unique_edges!, unique_e_indices, graph, graph_stucture, graph_data, parallel)

    mass_matrix = construct_mass_matrix(mmv_array, mme_array, graph_stucture)

    ODEFunction(nd!; mass_matrix=mass_matrix, syms=symbols)
end

# Here is not the best place for this error,
# but requires the least changes before the next refactor.
function prepare_edges(edges, g)
    @assert typeof(edges) <: Union{EdgeFunction,Vector}
    throw(ArgumentError("Graph type not recognized. Currently only SimpleGraph and SimpleDiGraph are supported."))
end

# If only a sinlge Function is given, not an Array of EdgeFunctions.
function prepare_edges(edge::EdgeFunction, g::SimpleGraph)
    if edge.coupling == :directed
        throw(ArgumentError("Coupling type of EdgeFunction not available for undirected Graphs"))
    elseif edge.coupling == :undefined
        @info("Reconstructing EdgeFunction with :undefined coupling type..")
        return reconstruct_edge(edge, :undirected)
    end
    return edge
end

function prepare_edges(edge::EdgeFunction, g::SimpleDiGraph)
    if edge.coupling ∈ (:symmetric, :antisymmetric, :undirected, :fiducial)
        throw(ArgumentError("Coupling type of EdgeFunction not available for directed Graphs"))
    elseif edge.coupling == :undefined
        @info("Reconstructing EdgeFunction with :undefined coupling type..")
        return reconstruct_edge(edge, :directed)
    end
    return edge
end

"""
    prepare_edges(edges, g::SimpleGraph)
"""
function prepare_edges(edges::Vector, g::SimpleGraph)
    ne(g) == length(edges) || error("If edges are given as Vector, size musst equal number of edges in g.")
    # Depending on the coupling we might get different eltypes
    new_edges = Vector{EdgeFunction}(undef, length(edges))
    infobool = true
    for (i, edge) in enumerate(edges)
        if edge.coupling == :directed
            throw(ArgumentError("Coupling type of edge $i not available for undirected Graphs"))
        elseif edge.coupling == :undefined
            if infobool
                @info("Reconstructing EdgeFuntions with :undefined coupling to have :undirected coupling. For optimal performance specify the coupling type during initialization of the edge function.")
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
    ne(g) == length(edges) || error("If edges are given as Vector, size musst equal number of edges in g.")
    # Depending on the coupling we might get different eltypes
    new_edges = Vector{EdgeFunction}(undef, length(edges))
    infobool = true
    for (i, edge) in enumerate(edges)
        if edge.coupling ∈ (:symmetric, :antisymmetric, :undirected, :fiducial)
            throw(ArgumentError("Coupling type of edge $i not available for directed Graphs"))
        elseif edge.coupling == :undefined
            if infobool
                @info("Reconstructing EdgeFuntions with :undefined coupling type..")
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
    let f = edge.f, dim = edge.dim, sym = edge.sym, name=edge.name, psym=edge.psym, obsf=edge.obsf, obssym=edge.obssym
        return StaticEdge(; f, dim, coupling, sym, name, psym, obsf, obssym)
    end
end
@inline function reconstruct_edge(edge::StaticDelayEdge, coupling::Symbol)
    let f = edge.f, dim = edge.dim, sym = edge.sym
        return StaticDelayEdge(; f=f,
                               dim=dim,
                               coupling=coupling,
                               sym=sym)
    end
end
@inline function reconstruct_edge(edge::ODEEdge, coupling::Symbol)
    error("Reconstruction of ODEEdges is not implemented at the moment.")
end

end # module
