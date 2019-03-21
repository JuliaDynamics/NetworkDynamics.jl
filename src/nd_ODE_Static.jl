module StaticLines

using Parameters
using LightGraphs
using LinearAlgebra

export static_lines

@with_kw struct static_lines
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
    no_parameters
end

function (d::static_lines)(dx, x, p, t)
    @views begin
        if d.no_parameters == true
            for i in 1:d.num_e
                d.edges![i](d.e_int[d.e_idx[i]], x[d.s_idx[i]], x[d.d_idx[i]], p, t)
            end
            for i in 1:d.num_v
                d.vertices![i](dx[d.v_idx[i]], x[d.v_idx[i]], d.e_s[i], d.e_d[i], p, t)
            end
        else
            for i in 1:d.num_e
                d.edges![i](d.e_int[d.e_idx[i]], x[d.s_idx[i]], x[d.d_idx[i]], p[d.num_v + i], t)
            end
            for i in 1:d.num_v
                d.vertices![i](dx[d.v_idx[i]], x[d.v_idx[i]], d.e_s[i], d.e_d[i], p[i], t)
            end
        end
    end
    nothing
end

"""
    dim_v is an array of the number of variables per vertex
    dim_e is an array of the number of variables per edge"""

function static_lines(vertices!, edges!, s_e, d_e, dim_v, dim_e, no_parameters)
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


    static_lines(
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
    e_d, # Array of Array of views on the variables in e_int of the edges that are destination of a vertex
    no_parameters)
end

function static_lines(vertices!, edges!, g::AbstractGraph, dim_v, dim_e; no_parameters = true)
    s_e = [src(e) for e in edges(g)]
    d_e = [dst(e) for e in edges(g)]
    static_lines(vertices!, edges!, s_e, d_e, dim_v, dim_e, no_parameters)
end

end #module
