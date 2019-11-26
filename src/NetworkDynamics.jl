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

function network_dynamics(vertices!::Array{T, 1}, edges!::Array{U, 1}, graph, p; x_prototype=zeros(1)) where {T <: ODEVertex, U <: StaticEdge}
    @assert length(vertices!) == length(vertices(graph))
    @assert length(edges!) == length(edges(graph))

    v_dims = [v.dim for v in vertices!]
    e_dims = [e.dim for e in edges!]

    v_array = similar(x_prototype, sum(v_dims))
    e_array = similar(x_prototype, sum(e_dims))

    symbols_v = [Symbol(vertices![i].sym[j],"_",i) for i in 1:length(vertices!) for j in 1:v_dims[i]]
    symbols_e = [Symbol(edges![i].sym[j],"_",i) for i in 1:length(edges!) for j in 1:e_dims[i]]
    symbols = symbols_v

    graph_stucture = GraphStruct(graph, v_dims, e_dims)

    graph_data = GraphData(v_array, e_array, symbols_v, symbols_e, graph_stucture)

    nd! = nd_ODE_Static(vertices!, edges!, graph, graph_stucture, graph_data)

    # This is disabled for now, because of various bugs in the solvers we need
    # The solvers that work with JacVecOperators don't work with mass matrics
    # right now
    # Jv = JacVecOperator(nd!, v_array, p, 0.0)

    # Construct mass matrix
    mass_matrix = construct_mass_matrix([v.mass_matrix for v in vertices!], graph_stucture)


    ODEFunction(nd!; mass_matrix = mass_matrix, syms=symbols) # jac_prototype=Jv
end


function network_dynamics(vertices!::Array{T, 1}, edges!::Array{U, 1}, graph, p; x_prototype=zeros(1)) where {T <: ODEVertex, U <: ODEEdge}
    @assert length(vertices!) == length(vertices(graph))
    @assert length(edges!) == length(edges(graph))

    v_dims = [v.dim for v in vertices!]
    e_dims = [e.dim for e in edges!]

    x_array = similar(x_prototype, sum(v_dims) + sum(e_dims))

    v_array = view(x_array, 1:sum(v_dims))
    e_array = view(x_array, sum(v_dims)+1:sum(v_dims)+sum(e_dims))

    symbols_v = [Symbol(vertices![i].sym[j],"_",i) for i in 1:length(vertices!) for j in 1:v_dims[i]]
    symbols_e = [Symbol(edges![i].sym[j],"_",i) for i in 1:length(edges!) for j in 1:e_dims[i]]
    symbols = vcat(symbols_v, symbols_e)

    graph_stucture = GraphStruct(graph, v_dims, e_dims)

    graph_data = GraphData(v_array, e_array, symbols_v, symbols_e, graph_stucture)

    nd! = nd_ODE_ODE(vertices!, edges!, graph, graph_stucture, graph_data)

    # This is disabled for now, because of various bugs in the solvers we need
    # The solvers that work with JacVecOperators don't work with mass matrics
    # right now
    # Jv = JacVecOperator(nd!, x_array, p, 0.0)

    # Construct mass matrix
    mass_matrix = construct_mass_matrix([v.mass_matrix for v in vertices!], [e.mass_matrix for e in edges!], graph_stucture)

    ODEFunction(nd!; mass_matrix = mass_matrix, syms=symbols) # jac_prototype=Jv
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
