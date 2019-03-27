module NetworkDynamics

include("nd_ODE_ODE_scalar.jl")
using .nd_ODE_ODE_scalar_mod
export nd_ODE_ODE_scalar

include("nd_ODE_Static_scalar.jl")
using .nd_ODE_Static_scalar_mod
export nd_ODE_Static_scalar

include("nd_ODE_ODE.jl")
using .nd_ODE_ODE_mod
export nd_ODE_ODE

include("nd_ODE_Static.jl")
using .nd_ODE_Static_mod
export nd_ODE_Static



include("Functions.jl")
using .NDFunctions
export StaticVertex
export StaticEdge
export ODEVertex
export ODEEdge
export VertexFunction
export EdgeFunction

export network_dynamics

export network_dynamics

using LinearAlgebra
using SparseArrays
using LightGraphs
using DifferentialEquations

function network_dynamics(vertices!::Array{ODEVertex,1}, edges!::Array{StaticEdge,1}, graph)
    @assert length(vertices!) == length(vertices(graph))
    @assert length(edges!) == length(edges(graph))
    massmatrix = nothing # Construct Mass Matrix from vertices! and edges!
    vertex_functions = [v.f! for v in vertices!]
    dim_v = [v.dim for v in vertices!]
    edge_functions = [e.f! for e in edges!]
    dim_e = [e.dim for e in edges!]
    dim_nd = sum(dim_v)
    massmatrix = sparse(1.0I,dim_nd,dim_nd) # Construct Mass Matrix from vertices! and edges!

    if all(dim_v .== 1) && all(dim_e .== 1)
        nd! = nd_ODE_Static_scalar(vertex_functions, edge_functions, graph)
        ODEFunction(nd!)
    else
        nd! = nd_ODE_Static(vertex_functions, edge_functions, graph, dim_v, dim_e)
        for i in 1:length(vertex_functions)
            for idx in nd!.v_idx
                massmatrix[idx,idx] = vertices![i].massmatrix
            end
        end
        ODEFunction(nd!,mass_matrix = massmatrix)
    end
end



function network_dynamics(vertices!::Array{ODEVertex}, edges!::Array{ODEEdge}, graph)
    @assert length(vertices!) == length(vertices(graph))
    @assert length(edges!) == length(edges(graph))
    vertex_functions = [v.f! for v in vertices!]
    dim_v = [v.dim for v in vertices!]
    edge_functions = [e.f! for e in edges!]
    dim_e = [e.dim for e in edges!]
    dim_nd = sum(dim_e) + sum(dim_v)
    massmatrix = sparse(1.0I,dim_nd,dim_nd) # Construct Mass Matrix from vertices! and edges!

    if all(dim_v .== 1) && all(dim_e .== 1)
        nd! = nd_ODE_ODE_scalar(vertex_functions, edge_functions, graph)
        ODEFunction(nd!)
    else
        nd! = nd_ODE_ODE(vertex_functions, edge_functions, graph, dim_v, dim_e)
        for i in 1:length(vertex_functions)
            for idx in nd!.v_idx
                massmatrix[idx,idx] = vertices![i].massmatrix
            end
        end
        for i in 1:length(edge_functions)
            for idx in nd!.e_x_idx
                massmatrix[idx,idx] = edges![i].massmatrix
            end
        end
        ODEFunction(nd!,mass_matrix = massmatrix)
    end
end

end
