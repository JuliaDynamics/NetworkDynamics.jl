module nd_ODE_Static_mod

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


@with_kw struct nd_ODE_Static{T}
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

function (d::nd_ODE_Static{T})(dx, x::AbstractArray{T}, p::Nothing, t) where T
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

function (d::nd_ODE_Static{T})(dx, x::AbstractArray{T}, p, t) where T
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

#= If the type of x doesn't match that of the internal variables we construct a
new internal variable system.
=#

function (d::nd_ODE_Static){T}(dx, x::AbstractArray{U}, p::Nothing, t) where T where U
    e_int = zeros(U, d.dim_nd)
    e_s = [[view(e_int, d.e_idx[i_e]) for i_e in 1:d.num_e if i_v == d.s_e[i_e]] for i_v in 1:d.num_v]
    e_d = [[view(e_int, d.e_idx[i_e]) for i_e in 1:d.num_e if i_v == d.d_e[i_e]] for i_v in 1:d.num_v]

    @views begin
    for i in 1:d.num_e
        d.edges![i](e_int[d.e_idx[i]], x[d.s_idx[i]], x[d.d_idx[i]], p, t)
    end
    for i in 1:d.num_v
        d.vertices![i](dx[d.v_idx[i]], x[d.v_idx[i]], e_s[i], e_d[i], p, t)
    end
    end
    nothing
end

function (d::nd_ODE_Static){T}(dx, x::AbstractArray{U}, p, t) where T where U
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

"""
    dim_v is an array of the number of variables per vertex
    dim_e is an array of the number of variables per edge"""

function nd_ODE_Static(vertices!, edges!, s_e, d_e, dim_v, dim_e)
    num_v = length(dim_v)
    num_e = length(dim_e)

    e_int = zeros(sum(dim_e))

    # x has length sum(dim_v)
    # v_idx is an Array of Array of indices of variables in x belonging to vertex

    counter = 1
    v_idx = [zeros(Int32, dim) for dim in dim_v]
    for i in 1:num_v
        v_idx[i] .= collect(counter:counter + dim_v[i] - 1)
        counter += dim_v[i]
    end

    counter = 1
    e_idx = [zeros(Int32, dim) for dim in dim_e]
    for i in 1:num_e
        e_idx[i] .= collect(counter:counter + dim_e[i] - 1)
        counter += dim_e[i]
    end

    # For every vertex, and for every edge, if the source of the edge is that vertex, take the view on the variables of that edge.
    # Thus e_s[i] is an array of views onto the variables of the edges for which i is the source.
    e_s = [[view(e_int, e_idx[i_e]) for i_e in 1:num_e if i_v == s_e[i_e]] for i_v in 1:num_v]
    e_d = [[view(e_int, e_idx[i_e]) for i_e in 1:num_e if i_v == d_e[i_e]] for i_v in 1:num_v]

    s_idx = [v_idx[s_e[i_e]] for i_e in 1:num_e]
    d_idx = [v_idx[d_e[i_e]] for i_e in 1:num_e]


    nd_ODE_Static(
    g
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

function nd_ODE_Static(vertices!, edges!, g::AbstractGraph, dim_v, dim_e)
    s_e = [src(e) for e in edges(g)]
    d_e = [dst(e) for e in edges(g)]
    nd_ODE_Static(vertices!, edges!, s_e, d_e, dim_v, dim_e)
end

struct StaticEdgeFunction
    nd_ODE_Static:: nd_ODE_Static
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
