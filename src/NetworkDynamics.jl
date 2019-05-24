module NetworkDynamics

using Reexport

include("Functions.jl")
@reexport using .NDFunctions

include("nd_ODE_ODE.jl")
@reexport using .nd_ODE_ODE_mod

include("nd_ODE_Static.jl")
@reexport using .nd_ODE_Static_mod

include("nd_DDE_DDE.jl")
@reexport using .nd_DDE_DDE_mod

include("Utilities.jl")
@reexport using .Utilities

include("NetworkStructures.jl")
@reexport using .NetworkStructures

export network_dynamics

using LinearAlgebra
using SparseArrays
using LightGraphs
using DifferentialEquations

#= network_dynamics: The Main Constructor of the Package. It takes Arrays of Vertex- and Edgefunction + a graph and
spits out an ODEFunction or DDEFunction. Others still need to be implemented. =#

function network_dynamics(vertices!::Array{ODEVertex,1}, edges!::Array{StaticEdge,1}, graph)
    @assert length(vertices!) == length(vertices(graph))
    @assert length(edges!) == length(edges(graph))

    dim_v = [v.dim for v in vertices!]
    dim_e = [e.dim for e in edges!]
    dim_nd = sum(dim_v)

    e_int = zeros(sum(dim_e))

    graph_stucture = create_graph_structure(graph, dim_v, dim_e, e_int)

    nd! = nd_ODE_Static(vertices!, edges!, graph, graph_stucture)

    # Construct mass matrix
    mass_matrix = construct_mass_matrix([v.mass_matrix for v in vertices!], dim_nd, graph_stucture)

    symbols = [Symbol(vertices![i].sym[j],"_",i) for i in 1:length(vertices!) for j in 1:dim_v[i]]

    ODEFunction(nd!,mass_matrix = mass_matrix,syms=symbols)
end


function network_dynamics(vertices!::Array{ODEVertex}, edges!::Array{ODEEdge}, graph)
    @assert length(vertices!) == length(vertices(graph))
    @assert length(edges!) == length(edges(graph))

    dim_v = [v.dim for v in vertices!]
    dim_e = [e.dim for e in edges!]
    dim_nd = sum(dim_e) + sum(dim_v)

    e_int = zeros(sum(dim_e))

    graph_stucture = create_graph_structure(graph, dim_v, dim_e, e_int)

    nd! = nd_ODE_ODE(vertices!, edges!, graph, graph_stucture)

    # Construct mass matrix
    mass_matrix = construct_mass_matrix([v.mass_matrix for v in vertices!], [e.mass_matrix for e in edges!], dim_nd, graph_stucture)

    vsymbols = [Symbol(vertices![i].sym[j],"_",i) for i in 1:length(vertices!) for j in 1:dim_v[i]]
    esymbols = [Symbol(edges![i].sym[j],"_",i) for i in 1:length(edges!) for j in 1:dim_e[i]]
    symbols = vcat(vsymbols,esymbols)

    ODEFunction(nd!,mass_matrix = mass_matrix,syms=symbols)
end

function network_dynamics(vertices!::Array{DDEVertex}, edges!::Array{DDEEdge}, graph)
    @assert 1==0
    # This isn't really ready for testing/production yet.
    vertex_functions = [v.f! for v in vertices!]
    dim_v = [v.dim for v in vertices!]
    edge_functions = [e.f! for e in edges!]
    tau_s = [v.tau_s for v in vertices!]
    tau_d = [v.tau_d for v in vertices!]
    dim_e = [e.dim for e in edges!]
    dim_nd = sum(dim_e) + sum(dim_v)

    vsymbols = [Symbol(vertices![i].sym[j],"_",i) for i in 1:length(vertices!) for j in 1:dim_v[i]]
    esymbols = [Symbol(edges![i].sym[j],"_",i) for i in 1:length(edges!) for j in 1:dim_e[i]]
    symbols = append!(vsymbols,esymbols)

    nd! = nd_DDE_DDE(vertex_functions, edge_functions, graph, dim_v, dim_e, tau_s, tau_d)

    # Construct mass matrix
    mmv_array = [v.mass_matrix for v in vertices!]
    mme_array = [e.mass_matrix for e in edges!]
    if all([mm == I for mm in mmv_array]) && all([mm == I for mm in mme_array])
        mass_matrix = I
    else
        mass_matrix = sparse(1.0I,dim_nd,dim_nd)
        for i in 1:length(vertex_functions)
            if vertices![i].mass_matrix != I
                mass_matrix[nd!.v_idx[i],nd!.v_idx[i]] .= vertices![i].mass_matrix
            end
        end
        for i in 1:length(edge_functions)
            if edges![i].mass_matrix != I
                mass_matrix[nd!.e_x_idx[i],nd!.e_x_idx[i]] = edges![i].mass_matrix
            end
        end
    end

    DDEFunction(nd!,mass_matrix = mass_matrix, syms = symbols)
end


function network_dynamic(vertices!,  edges!, graph)
    try
        va! = Array{NDFunctions.VertexFunction}(vertices!)
    catch err
        println("Cannot convert the vertices to an Array{VertexFunction}!")
        println(err)
        return nothing
    end
    try
        ea! = Array{NDFunctions.EdgeFunction}(edges!)
    catch err
        println("Cannot convert the edges to an Array{EdgeFunction}!")
        println(err)
        return nothing
    end
    network_dynamic(va!,  ea!, graph)
end

function network_dynamics(vertices!::Array{VertexFunction}, edges!::Array{EdgeFunction}, graph)
    @assert length(vertices!) == length(vertices(graph))
    @assert length(edges!) == length(edges(graph))

    contains_delays = false
    contains_stochastic = false
    contains_dyn_edge = false

    for v in vertices!
        if typeof(v) == DDEVertex
            contains_delays = true
        end
        # if typeof(v) == SDEVertex
        #     contains_stochastic = true
        # end
    end
    for e in edges!
        if typeof(e) == DDEEdge
            contains_delays = true
        end
        if typeof(e) == ODEEdge
            contains_dyn_edge = true
        end
        # if typeof(e) == SDEEdge
        #     contains_stochastic = true
        # end
    end

    # ToDo... more logic about what to construct when... A lot of this should
    # be covered by the casts.
    if contains_delay && contains_stochastic
        println("Stochasticity and delay are not supported together.")
        return nothing
    elseif contains_delay
        return network_dynamics(Array{DDEVertex}(vertices!),Array{DDEEdge}(edges!),graph)
    end
    nothing
end

end # module
