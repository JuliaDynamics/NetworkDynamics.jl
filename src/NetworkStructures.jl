module NetworkStructures

using LightGraphs
#= This module contains helper functions that calculate the index structures
and view arrays that NEtwork Dynamics makes use of.
=#

export create_idxs, create_network_data


function create_idxs(dims, counter=1)
    idxs = [zeros(Int32, dim) for dim in dims]
    for (i, dim) in enumerate(dims)
        idxs[i] .= collect(counter:(counter + dim - 1))
        counter += dim
    end
    idxs
end

function create_network_data(g, dim_v, dim_e, e_int)
    s_e = [src(e) for e in edges(g)]
    d_e = [dst(e) for e in edges(g)]

    num_v = nv(g)
    num_e = ne(g)

    # x has length sum(dim_v)

    # v_idx is an Array of Array of indices of variables in x belonging to vertices
    # e_idx is an Array of Array of indices of variables in e_int belonging to edges
    # e_x_idx is an Array of Array of indices of variables in x belonging to edges if edges are dynamic.

    v_idx = create_idxs(dim_v)
    e_idx = create_idxs(dim_e)
    e_x_idx = [e_idx[i] .+ sum(dim_v) for i in 1:num_e]

    # For every vertex, and for every edge, if the source of the edge is that vertex, take the view on the variables of that edge.
    # Thus e_s[i] is an array of views onto the variables of the edges for which i is the source.
    e_s = [[view(e_int, e_idx[i_e]) for i_e in 1:num_e if i_v == s_e[i_e]] for i_v in 1:num_v]
    e_d = [[view(e_int, e_idx[i_e]) for i_e in 1:num_e if i_v == d_e[i_e]] for i_v in 1:num_v]

    s_idx = [v_idx[s_e[i_e]] for i_e in 1:num_e]
    d_idx = [v_idx[d_e[i_e]] for i_e in 1:num_e]

    (num_v, # Number of vertices
    num_e, # Number of edges
    e_int, # Variables living on edges
    e_idx, # Array of Array of indices of variables in e_int belonging to edges
    e_x_idx, # Array of Array of indices of variables in x belonging to edges if edges are dynamic
    s_idx, # Array of Array of indices of variables in x belonging to source vertex of edge
    d_idx, # Array of Array of indices of variables in x belonging to destination vertex of edge
    v_idx, # Array of Array of indices of variables in x belonging to vertex
    e_s, # Array of Array of views on the variables in e_int of the edges that are source of a vertex
    e_d) # Array of Array of views on the variables in e_int of the edges that are destination of a vertex
end

end # module
