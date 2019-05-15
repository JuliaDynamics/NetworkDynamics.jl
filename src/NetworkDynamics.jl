module NetworkDynamics

include("nd_ODE_ODE.jl")
using .nd_ODE_ODE_mod
export nd_ODE_ODE

include("nd_ODE_Static.jl")
using .nd_ODE_Static_mod
export nd_ODE_Static
export StaticEdgeFunction

include("nd_DDE_DDE.jl")
using .nd_DDE_DDE_mod
export nd_DDE_DDE

include("Functions.jl")
using .NDFunctions
export StaticVertex
export StaticEdge
export ODEVertex
export ODEEdge
export VertexFunction
export EdgeFunction
export DDEVertex
export DDEEdge

include("Utilities.jl")
using .Utilities
export RootRhs
export find_valid_ic

export network_dynamics

using LinearAlgebra
using SparseArrays
using LightGraphs
using DifferentialEquations

#= network_dynamics: The Main Constructor of the Package. It takes Arrays of Vertex- and Edgefunction + a graph and
spits out an ODEFunction or DDEFunction. Others still need to be implemented. =#

function network_dynamics(vertices!::Array{VertexFunction}, edges!::Array{EdgeFunction}, graph)
    @assert length(vertices!) = length(vertices(graph))
    for i in 1:length(vertices!)
        if typeof(vertices![i]) == DDEVertex
            for i in 1:length(vertices!)
                vertices![i] = DDEVertex(vertices![i])
            end
            break
        end
    end
    for i in 1:length(edges!)
        if typeof(edges![i]) == DDEEdge
            for i in 1:length(edges!)
                edges![i] = DDEEdge(edges![i])
            end
            break
        end
    end
    network_dynamics(vertices!,edges!,graph)
end

function network_dynamics(vertices!::Array{ODEVertex,1}, edges!::Array{StaticEdge,1}, graph)

    vertex_functions = [v.f! for v in vertices!]
    dim_v = [v.dim for v in vertices!]
    edge_functions = [e.f! for e in edges!]
    dim_e = [e.dim for e in edges!]
    dim_nd = sum(dim_v)

    symbols = [Symbol(vertices![i].sym[j],"_",i) for i in 1:length(vertices!) for j in 1:dim_v[i]]

    nd! = nd_ODE_Static(vertex_functions, edge_functions, graph, dim_v, dim_e)

    # Construct mass matrix
    mm_array = [v.mass_matrix for v in vertices!]
    if all(mm_array .== I)
        mass_matrix = I
    else
        mass_matrix = sparse(1.0I,dim_nd,dim_nd)
        for i in 1:length(vertex_functions)
            if vertices![i].mass_matrix != I
                mass_matrix[nd!.v_idx[i],nd!.v_idx[i]] .= vertices![i].mass_matrix
            end
        end
    end

    ODEFunction(nd!,mass_matrix = mass_matrix,syms=symbols)
end

function network_dynamics(vertices!::Array{ODEVertex}, edges!::Array{ODEEdge}, graph)

    vertex_functions = [v.f! for v in vertices!]
    dim_v = [v.dim for v in vertices!]
    edge_functions = [e.f! for e in edges!]
    dim_e = [e.dim for e in edges!]
    dim_nd = sum(dim_e) + sum(dim_v)

    vsymbols = [Symbol(vertices![i].sym[j],"_",i) for i in 1:length(vertices!) for j in 1:dim_v[i]]
    esymbols = [Symbol(edges![i].sym[j],"_",i) for i in 1:length(edges!) for j in 1:dim_e[i]]
    symbols = append!(vsymbols,esymbols)

    nd! = nd_ODE_ODE(vertex_functions, edge_functions, graph, dim_v, dim_e)

    # Construct mass matrix
    mmv_array = [v.mass_matrix for v in vertices!]
    mme_array = [e.mass_matrix for e in edges!]
    if all(mme_array .== I) && all(mmv_array .== I)
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

    ODEFunction(nd!,mass_matrix = mass_matrix,syms = symbols)
end

function network_dynamics(vertices!::Array{DDEVertex}, edges!::Array{DDEEdge}, graph)

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
    if all(mme_array .== I) && all(mmv_array .== I)
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
