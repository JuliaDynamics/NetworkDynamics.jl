module NetworkStructures

using LightGraphs
using LinearAlgebra
using SparseArrays
#= This module contains helper functions that calculate the index structures
and view arrays that Network Dynamics makes use of.

The key structure is the GraphStructure that contains the internal variables of
the edges, as well as all the views into them, and the complete set of indices.
=#

export create_idxs, create_graph_structure, GraphStructure, construct_mass_matrix, GraphStructure_like

import Base.getindex


struct GS_e_int_View{G}
    gs::G
    idx_offset::Int
end

function getindex(gev::GS_e_int_View, idx::Int)
    gev.gs.e_int[idx + gs_x_view.idx_offset]
end


struct GS_x_View{G}
    gs::G
    idx_offset::Int
end

function getindex(gxv::GS_x_View, idx::Int)
    gxv.gs.x[idx + gs_x_view.idx_offset]
end

#
# mutable struct GS{T}
#     x::T
#     gs_x_view::GS_x_View{GS{T}}
#     function GS(arr::T) where T
#         gs = new{T}(arr)
#         gs.gs_x_view = GS_x_View{GS{T}}(gs, 1)
#         gs
#     end
# end

const Idx = UnitRange{Int}

struct GraphStructure{T_e_int}
    num_v::Int # Number of vertices
    num_e::Int # Number of edges
    e_int::Array{T_e_int, 1} # Variables living on edges
    e_idx::Array{Idx, 1} # Array of indices of variables in e_int belonging to edges
    e_x_idx::Array{Idx, 1} # Array of indices of variables in x belonging to edges if edges are dynamic
    s_idx::Array{Idx, 1} # Array of indices of variables in x belonging to source vertex of edge
    d_idx::Array{Idx, 1} # Array of indices of variables in x belonging to destination vertex of edge
    v_idx::Array{Idx, 1} # Array of indices of variables in x belonging to vertex
    # Array of Array of views on the variables in e_int of the edges that are source of a vertex
    e_s::Array{Array{SubArray{T_e_int,1,Array{T_e_int,1},Tuple{Idx},true},1},1}
    e_d::Array{Array{SubArray{T_e_int,1,Array{T_e_int,1},Tuple{Idx},true},1},1}
end

function GraphStructure_like(GS, x, g)
    s_e = [src(e) for e in edges(g)]
    d_e = [dst(e) for e in edges(g)]

    e_int = zeros(eltype(x), length(GS.e_int))

    e_s = [[view(e_int, GS.e_idx[i_e]) for i_e in 1:GS.num_e if i_v == s_e[i_e]] for i_v in 1:GS.num_v]
    e_d = [[view(e_int, GS.e_idx[i_e]) for i_e in 1:GS.num_e if i_v == d_e[i_e]] for i_v in 1:GS.num_v]

    GraphStructure(GS.num_v, # Number of vertices
    GS.num_e, # Number of edges
    e_int, # Variables living on edges
    GS.e_idx, # Array of Array of indices of variables in e_int belonging to edges
    GS.e_x_idx, # Array of Array of indices of variables in x belonging to edges if edges are dynamic
    GS.s_idx, # Array of Array of indices of variables in x belonging to source vertex of edge
    GS.d_idx, # Array of Array of indices of variables in x belonging to destination vertex of edge
    GS.v_idx, # Array of Array of indices of variables in x belonging to vertex
    e_s, # Array of Array of views on the variables in e_int of the edges that are source of a vertex
    e_d) # Array of Array of views on the variables in e_int of the edges that are destination of a vertex
end

function create_idxs(dims, counter=1)::Array{Idx, 1}
    idxs = [1:1 for dim in dims]
    for (i, dim) in enumerate(dims)
        idxs[i] = counter:(counter + dim - 1)
        counter += dim
    end
    idxs
end

function create_graph_structure(g, dim_v, dim_e, e_int)
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

    GraphStructure(num_v, # Number of vertices
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

function construct_mass_matrix(mmv_array, dim_nd, gs::GraphStructure)
    if all([mm == I for mm in mmv_array])
        mass_matrix = I
    else
        mass_matrix = sparse(1.0I,dim_nd,dim_nd)
        for (i, mm) in enumerate(mmv_array)
            if mm != I
                mass_matrix[gs.v_idx[i],gs.v_idx[i]] .= mm
            end
        end
    end
    mass_matrix
end

function construct_mass_matrix(mmv_array, mme_array, dim_nd, gs::GraphStructure)
    if all([mm == I for mm in mmv_array]) && all([mm == I for mm in mme_array])
        mass_matrix = I
    else
        mass_matrix = sparse(1.0I,dim_nd,dim_nd)
        for (i, mm) in enumerate(mmv_array)
            if mm != I
                mass_matrix[gs.v_idx[i],gs.v_idx[i]] .= mm
            end
        end
        for (i, mm) in enumerate(mme_array)
            if mm != I
                mass_matrix[gs.e_x_idx[i],gs.e_x_idx[i]] .= mm
            end
        end
    end
    mass_matrix
end

end # module
