module nd_ODE_Static_mod

include("NetworkStructures.jl")
using .NetworkStructures

using Parameters
using LightGraphs
using LinearAlgebra

export nd_ODE_Static

#= nd_ODE_Static constructs a (dx,x,p,t)-function from an Array of functions for the vertices,
 edges as well as a graph.
The arguments of the vertex functions must be of the form (dv,v,e_s,e_d,p,t),
where dv is the vertex variable derivative, v the vertex variable and e_s and e_d Arrays of edge variables that
have the vertex as source and destination respectively. p and t are as usual.
The arguments of the edge functions must be of the form (e,v_s,v_d,p,t),
where e is the edge variable and v_s and v_d the vertex variables of the vertices
the edge has as source and destination respectively.
This works for multi-dimensional variables as well. =#


@with_kw struct nd_ODE_Static
    g
    edges!
    vertices!
    num_v # Number of vertices
    num_e # Number of edges
    e_int # Variables living on edges
    e_idx # Array of Array of indices of variables in e_int belonging to edges
    s_idx # Array of Array of indices of variables in x belonging to source vertex of edge
    d_idx # Array of Array of indices of variables in x belonging to destination vertex of edge
    v_idx # Array of Array of indices of variables in x belonging to vertex
    e_s # Array of Array of views on the variables in e_int of the edges that are source of a vertex
    e_d # Array of Array of views on the variables in e_int of the edges that are destination of a vertex
end

function (d::nd_ODE_Static)(dx, x, p::Nothing, t)
    @views begin
    for i in 1:d.num_e
        d.edges![i](d.e_int[d.e_idx[i]], x[d.s_idx[i]], x[d.d_idx[i]], p, t)
    end
    for i in 1:d.num_v
        d.vertices![i](dx[d.v_idx[i]], x[d.v_idx[i]], d.e_s[i], d.e_d[i], p, t)
    end
    end
    nothing
end

function (d::nd_ODE_Static)(dx, x, p, t)
    @views begin
    for i in 1:d.num_e
        d.edges![i](d.e_int[d.e_idx[i]], x[d.s_idx[i]], x[d.d_idx[i]], p[d.num_v + i], t)
    end
    for i in 1:d.num_v
        d.vertices![i](dx[d.v_idx[i]], x[d.v_idx[i]], d.e_s[i], d.e_d[i], p[i], t)
    end
    end
    nothing
end

function nd_ODE_Static(vertices!, edges!, g, dim_v, dim_e; T=Float64)

    e_int=zeros(T, sum(dim_e))

    (num_v, num_e, e_int, e_idx, e_x_idx, s_idx, d_idx, v_idx, e_s, e_d) =
        create_network_data(g, dim_v, dim_e, e_int)

    nd_ODE_Static{T}(
    g,
    edges!,
    vertices!,
    num_v, # Number of vertices
    num_e, # Number of edges
    e_int, # Variables living on edges
    e_idx, # Array of Array of indices of variables in e_int belonging to edges
    s_idx, # Array of Array of indices of variables in x belonging to source vertex of edge
    d_idx, # Array of Array of indices of variables in x belonging to destination vertex of edge
    v_idx, # Array of Array of indices of variables in x belonging to vertex
    e_s, # Array of Array of views on the variables in e_int of the edges that are source of a vertex
    e_d # Array of Array of views on the variables in e_int of the edges that are destination of a vertex
    )
end




#= The struct comes in a version where edges and vertices are not arrays but
homogenous across the network.
=#

@with_kw struct nd_ODE_Static_hom{T}
    g
    edges!
    vertices!
    num_v # Number of vertices
    num_e # Number of edges
    e_int::AbstractArray{T} # Variables living on edges
    e_idx # Array of Array of indices of variables in e_int belonging to edges
    s_idx # Array of Array of indices of variables in x belonging to source vertex of edge
    d_idx # Array of Array of indices of variables in x belonging to destination vertex of edge
    v_idx # Array of Array of indices of variables in x belonging to vertex
    e_s # Array of Array of views on the variables in e_int of the edges that are source of a vertex
    e_d # Array of Array of views on the variables in e_int of the edges that are destination of a vertex
end

function (d::nd_ODE_Static_hom{T})(dx, x::AbstractArray{T}, p::Nothing, t) where T
    @views begin
    for i in 1:d.num_e
        d.edges!(d.e_int[d.e_idx[i]], x[d.s_idx[i]], x[d.d_idx[i]], p, t)
    end
    for i in 1:d.num_v
        d.vertices!(dx[d.v_idx[i]], x[d.v_idx[i]], d.e_s[i], d.e_d[i], p, t)
    end
    end
    nothing
end

function (d::nd_ODE_Static_hom{T})(dx, x::AbstractArray{T}, p, t) where T
    @views begin
    for i in 1:d.num_e
        d.edges!(d.e_int[d.e_idx[i]], x[d.s_idx[i]], x[d.d_idx[i]], p[d.num_v + i], t)
    end
    for i in 1:d.num_v
        d.vertices!(dx[d.v_idx[i]], x[d.v_idx[i]], d.e_s[i], d.e_d[i], p[i], t)
    end
    end
    nothing
end



struct StaticEdgeFunction
    nd_ODE_Static::nd_ODE_Static
end

function StaticEdgeFunction(vertices!, edges!, g::AbstractGraph)
    dim_v = [v.dim for v in vertices!]
    dim_e = [e.dim for e in edges!]
    vertex_functions = [v.f! for v in vertices!]
    edge_functions = [e.f! for e in edges!]
    StaticEdgeFunction(nd_ODE_Static(vertex_functions, edge_functions, g, dim_v, dim_e))
end

function (sef::StaticEdgeFunction)(x, p::Nothing, t)
    d = sef.nd_ODE_Static
    @views begin
        for i in 1:d.num_e
            d.edges![i](d.e_int[d.e_idx[i]], x[d.s_idx[i]], x[d.d_idx[i]], p, t)
        end
    end
    (d.e_s, d.e_d)
end

function (sef::StaticEdgeFunction)(x, p, t)
    d = sef.nd_ODE_Static
    @views begin
        for i in 1:d.num_e
            d.edges![i](d.e_int[d.e_idx[i]], x[d.s_idx[i]], x[d.d_idx[i]], p[d.num_v + i], t)
        end
    end
    (d.e_s, d.e_d)
end

end #module
