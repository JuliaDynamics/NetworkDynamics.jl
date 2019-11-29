module NetworkDynamics

using Reexport

include("Functions.jl")
@reexport using .NDFunctions

include("NetworkStructures.jl")
@reexport using .NetworkStructures

include("nd_ODE_ODE.jl")
@reexport using .nd_ODE_ODE_mod

include("nd_ODE_Static.jl")
@reexport using .nd_ODE_Static_mod

include("Utilities.jl")
@reexport using .Utilities

include("SimpleAPI.jl")

export network_dynamics

using LinearAlgebra
using SparseArrays
using LightGraphs
using DiffEqBase
using DiffEqOperators

#= network_dynamics: The Main Constructor of the Package. It takes Arrays of Vertex- and Edgefunction + a graph and
spits out an ODEFunction or DDEFunction. Others still need to be implemented. =#

"""
This function assembles the arrays that hold the structural information
of the individual vertices and edges: Their dimensions, the symbols of their
variables and their mass matrices. Used internally by the network_dynamics
constructors.
"""
function collect_ve_info(vertices!, edges!, graph)
    if vertices! isa Array
        @assert length(vertices!) == length(vertices(graph))
        v_dims = [v.dim for v in vertices!]
        symbols_v = [Symbol(vertices![i].sym[j],"_",i) for i in 1:length(vertices!) for j in 1:v_dims[i]]
        mmv_array = [v.mass_matrix for v in vertices!]
    else
        v_dims = [vertices!.dim for v in vertices(graph)]
        symbols_v = [Symbol(vertices!.sym[j],"_",i) for i in 1:nv(graph) for j in 1:v_dims[i]]
        mmv_array = [vertices!.mass_matrix for v in vertices(graph)]
    end

    if edges! isa Array
        @assert length(edges!) == length(edges(graph))
        e_dims = [e.dim for e in edges!]
        symbols_e = [Symbol(edges![i].sym[j],"_",i) for i in 1:length(edges!) for j in 1:e_dims[i]]
        mme_array = [e.mass_matrix for e in edges!]
    else
        e_dims = [edges!.dim for e in edges(graph)]
        symbols_e = [Symbol(edges!.sym[j],"_",i) for i in 1:ne(graph) for j in 1:e_dims[i]]
        mme_array = [edges!.mass_matrix for e in edges(graph)]
    end

    v_dims, e_dims, symbols_v, symbols_e, mmv_array, mme_array
end

function network_dynamics(vertices!::Union{Array{T, 1}, T}, edges!::Union{Array{U, 1}, U}, graph; calc_JVOp = (false, nothing), x_prototype=zeros(1)) where {T <: ODEVertex, U <: StaticEdge}
    v_dims, e_dims, symbols_v, symbols_e, mmv_array, mme_array = collect_ve_info(vertices!, edges!, graph)

    v_array = similar(x_prototype, sum(v_dims))
    e_array = similar(x_prototype, sum(e_dims))

    symbols = symbols_v

    graph_stucture = GraphStruct(graph, v_dims, e_dims)

    graph_data = GraphData(v_array, e_array, symbols_v, symbols_e, graph_stucture)

    nd! = nd_ODE_Static(vertices!, edges!, graph, graph_stucture, graph_data)

    if calc_JVOp[1]
        Jv = JacVecOperator(nd!, x_array, calc_JVOp[2], 0.0)
    else
        Jv = nothing
    end

    # Construct mass matrix
    mass_matrix = construct_mass_matrix(mmv_array, graph_stucture)

    ODEFunction(nd!; mass_matrix = mass_matrix, syms=symbols, jac_prototype=Jv)
end


function network_dynamics(vertices!::Union{Array{T, 1}, T}, edges!::Union{Array{U, 1}, U}, graph; calc_JVOp = (false, nothing), x_prototype=zeros(1)) where {T <: ODEVertex, U <: ODEEdge}
    v_dims, e_dims, symbols_v, symbols_e, mmv_array, mme_array = collect_ve_info(vertices!, edges!, graph)

    x_array = similar(x_prototype, sum(v_dims) + sum(e_dims))

    v_array = view(x_array, 1:sum(v_dims))
    e_array = view(x_array, sum(v_dims)+1:sum(v_dims)+sum(e_dims))

    symbols = vcat(symbols_v, symbols_e)

    graph_stucture = GraphStruct(graph, v_dims, e_dims)

    graph_data = GraphData(v_array, e_array, symbols_v, symbols_e, graph_stucture)

    nd! = nd_ODE_ODE(vertices!, edges!, graph, graph_stucture, graph_data)

    if calc_JVOp[1]
        Jv = JacVecOperator(nd!, x_array, calc_JVOp[2], 0.0)
    else
        Jv = nothing
    end

    # Construct mass matrix
    mass_matrix = construct_mass_matrix(mmv_array, mme_array, graph_stucture)

    ODEFunction(nd!; mass_matrix = mass_matrix, syms=symbols, jac_prototype=Jv)
end

function network_dynamics(vertices!,  edges!, graph, p)
    try
        Array{VertexFunction}(vertices!)
    catch err
        println("Cannot convert the vertices to an Array{VertexFunction}!")
        println(err)
        return nothing
    end
    try
        Array{EdgeFunction}(edges!)
    catch err
        println("Cannot convert the edges to an Array{EdgeFunction}!")
        println(err)
        return nothing
    end
    va! = Array{VertexFunction}(vertices!)
    ea! = Array{EdgeFunction}(edges!)
    network_dynamics(va!,  ea!, graph, p)
end

function network_dynamics(vertices!::Array{VertexFunction}, edges!::Array{EdgeFunction}, graph, p)
    @assert length(vertices!) == length(vertices(graph))
    @assert length(edges!) == length(edges(graph))

    contains_dyn_edge = false

    for e in edges!
        if isa(e, ODEEdge)
            contains_dyn_edge = true
        end
    end

    if contains_dyn_edge
        return network_dynamics(Array{ODEVertex}(vertices!),Array{ODEEdge}(edges!), graph, p)
    else
        return network_dynamics(Array{ODEVertex}(vertices!),Array{StaticEdge}(edges!), graph, p)
    end
    nothing
end

end # module
