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

include("nd_DDE_DDE_scalar.jl")
using .nd_DDE_DDE_scalar_mod
export nd_DDE_DDE_scalar

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

    massmatrix = nothing # Construct Mass Matrix from vertices! and edges!
    vertex_functions = [v.f! for v in vertices!]
    dim_v = [v.dim for v in vertices!]
    edge_functions = [e.f! for e in edges!]
    dim_e = [e.dim for e in edges!]
    dim_nd = sum(dim_v)
    massmatrix = sparse(1.0I,dim_nd,dim_nd)
    symbols = [Symbol(vertices![i].sym[j],"_",i) for i in 1:length(vertices!) for j in 1:dim_v[i]]

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
        ODEFunction(nd!,mass_matrix = massmatrix,syms=symbols)
    end
end

function network_dynamics(vertices!::Array{ODEVertex}, edges!::Array{ODEEdge}, graph)

    vertex_functions = [v.f! for v in vertices!]
    dim_v = [v.dim for v in vertices!]
    edge_functions = [e.f! for e in edges!]
    dim_e = [e.dim for e in edges!]
    dim_nd = sum(dim_e) + sum(dim_v)
    massmatrix = sparse(1.0I,dim_nd,dim_nd) # Construct Mass Matrix from vertices! and edges!
    vsymbols = [Symbol(vertices![i].sym[j],"_",i) for i in 1:length(vertices!) for j in 1:dim_v[i]]
    esymbols = [Symbol(edges![i].sym[j],"_",i) for i in 1:length(edges!) for j in 1:dim_e[i]]
    symbols = append!(vsymbols,esymbols)
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
        ODEFunction(nd!,mass_matrix = massmatrix,syms = symbols)
    end
end

function network_dynamics(vertices!::Array{DDEVertex}, edges!::Array{DDEEdge}, graph)

    vertex_functions = [v.f! for v in vertices!]
    dim_v = [v.dim for v in vertices!]
    edge_functions = [e.f! for e in edges!]
    tau_s = [v.tau_s for v in vertices!]
    tau_d = [v.tau_d for v in vertices!]
    dim_e = [e.dim for e in edges!]
    dim_nd = sum(dim_e) + sum(dim_v)
    massmatrix = sparse(1.0I,dim_nd,dim_nd) # Construct Mass Matrix from vertices! and edges!
    vsymbols = [Symbol(vertices![i].sym[j],"_",i) for i in 1:length(vertices!) for j in 1:dim_v[i]]
    esymbols = [Symbol(edges![i].sym[j],"_",i) for i in 1:length(edges!) for j in 1:dim_e[i]]
    symbols = append!(vsymbols,esymbols)

    if all(dim_v .== 1) && all(dim_e .== 1)
        nd! = nd_DDE_DDE_scalar(vertex_functions, edge_functions, graph, tau_s, tau_d)
        DDEFunction(nd!)
    else
        nd! = nd_DDE_DDE(vertex_functions, edge_functions, graph, dim_v, dim_e, tau_s, tau_d)
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
        DDEFunction(nd!,mass_matrix = massmatrix, syms = symbols)
    end
end
end
