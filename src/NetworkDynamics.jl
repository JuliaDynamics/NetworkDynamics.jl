module NetworkDynamics

using Reexport

include("Functions.jl")
@reexport using .NDFunctions

include("nd_ODE_ODE.jl")
@reexport using .nd_ODE_ODE_mod

include("nd_ODE_Static.jl")
@reexport using .nd_ODE_Static_mod

include("Utilities.jl")
@reexport using .Utilities

include("NetworkStructures.jl")
@reexport using .NetworkStructures

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

    graph_stucture = nd_ODE_Static_mod.GraphStruct(graph, v_dims, e_dims)

    graph_data = nd_ODE_Static_mod.GraphData{typeof(v_array)}(v_array, e_array, graph_stucture)

    nd! = nd_ODE_Static(vertices!, edges!, graph, graph_stucture, graph_data)

    Jv = JacVecOperator(nd!, v_array, p, 0.0)

    # Construct mass matrix
    mass_matrix = construct_mass_matrix([v.mass_matrix for v in vertices!], sum(v_dims), graph_stucture)

    symbols = [Symbol(vertices![i].sym[j],"_",i) for i in 1:length(vertices!) for j in 1:v_dims[i]]


    ODEFunction(nd!; jac_prototype=Jv, mass_matrix = mass_matrix, syms=symbols)
end


function network_dynamics(vertices!::Array{T, 1}, edges!::Array{U, 1}, graph, p; x_prototype=zeros(1)) where T <: ODEVertex where U <: ODEEdge
    @assert length(vertices!) == length(vertices(graph))
    @assert length(edges!) == length(edges(graph))

    v_dims = [v.dim for v in vertices!]
    e_dims = [e.dim for e in edges!]

    x_array = similar(x_prototype, sum(v_dims) + sum(e_dims))

    v_array = view(x_prototype, 1:sum(v_dims))
    e_array = view(x_prototype, sum(v_dims)+1:sum(v_dims)+sum(e_dims))

    graph_stucture = GraphStruct(graph, v_dims, e_dims)

    graph_data = GraphData{typeof(v_array)}(v_array, e_array, graph_stucture)

    nd! = nd_ODE_ODE(vertices!, edges!, graph, graph_stucture, graph_data)

    Jv = JacVecOperator(nd!, v_array, p, 0.0)

    # Construct mass matrix
    # correct this mass_matrix = construct_mass_matrix([v.mass_matrix for v in vertices!], sum(v_dims), graph_stucture)

    # correct this  symbols = [Symbol(vertices![i].sym[j],"_",i) for i in 1:length(vertices!) for j in 1:v_dims[i]]


    #ODEFunction(nd!; jac_prototype=Jv, mass_matrix = mass_matrix, syms=symbols)
    nothing
end

#
# function network_dynamics(vertices!,  edges!, graph)
#     try
#         va! = Array{NDFunctions.VertexFunction}(vertices!)
#     catch err
#         println("Cannot convert the vertices to an Array{VertexFunction}!")
#         println(err)
#         return nothing
#     end
#     try
#         ea! = Array{NDFunctions.EdgeFunction}(edges!)
#     catch err
#         println("Cannot convert the edges to an Array{EdgeFunction}!")
#         println(err)
#         return nothing
#     end
#     network_dynamics(va!,  ea!, graph)
# end
#
# function network_dynamics(vertices!::Array{VertexFunction}, edges!::Array{EdgeFunction}, graph)
#     @assert length(vertices!) == length(vertices(graph))
#     @assert length(edges!) == length(edges(graph))
#
#     contains_delays = false
#     contains_stochastic = false
#     contains_dyn_edge = false
#
#     for v in vertices!
#         # if typeof(v) == DDEVertex
#         #     contains_delays = true
#         # end
#         # if typeof(v) == SDEVertex
#         #     contains_stochastic = true
#         # end
#     end
#     for e in edges!
#         # if typeof(e) == DDEEdge
#         #     contains_delays = true
#         # end
#         if typeof(e) == ODEEdge
#             contains_dyn_edge = true
#         end
#         # if typeof(e) == SDEEdge
#         #     contains_stochastic = true
#         # end
#     end
#
#     # ToDo... more logic about what to construct when... A lot of this should
#     # be covered by the casts.
#     if contains_delay && contains_stochastic
#         println("Stochasticity and delay are not supported together.")
#         return nothing
#     elseif contains_delay
#         # return network_dynamics(Array{DDEVertex}(vertices!),Array{DDEEdge}(edges!),graph)
#     end
#     nothing
# end

end # module
