"""
    This module contains the logic that calculate the index structures
    and data access structs that Network Dynamics makes use of.

    The key structure is the GraphData structure that allows accessing data on
    vertices and edges of the graph in an efficient manner. The neccessary indices
    are precomputed in GraphStructure.
"""
module NetworkStructures

using LightGraphs
using LinearAlgebra
using SparseArrays


# We need rather complicated sets of indices into the arrays that hold the
# vertex and the edge variables. We precompute everything we can and store it
# in GraphStruct.

export create_idxs, create_offsets, GraphStruct, GraphData, EdgeData, VertexData, construct_mass_matrix

const Idx = UnitRange{Int}

"""
Create indices for stacked array of dimensions dims
"""
function create_idxs(dims; counter=1)::Array{Idx, 1}
    idxs = [1:1 for dim in dims]
    for (i, dim) in enumerate(dims)
        idxs[i] = counter:(counter + dim - 1)
        counter += dim
    end
    idxs
end

"""
Create offsets for stacked array of dimensions dims
"""
function create_offsets(dims; counter=0)::Array{Int, 1}
    offs = [1 for dim in dims]
    for (i, dim) in enumerate(dims)
        offs[i] = counter
        counter += dim
    end
    offs
end

"""
Create indexes for stacked array of dimensions dims using the offsets offs
"""
function create_idxs(offs, dims)::Array{Idx, 1}
    idxs = [1+off:off+dim for (off, dim) in zip(offs, dims)]
end

"""
This struct holds the offsets and indices for all relevant aspects of the graph
The assumption is that there will be two arrays, one for the vertex variables
and one for the edge variables.

The graph structure is encoded in the source and destination relationships s_e
and d_e. These are arrays that hold the node that is the source/destination of
the indexed edge. Thus ``e_i = (s_e[i], d_e[i])``
"""
struct GraphStruct
    num_v::Int
    num_e::Int
    v_dims::Array{Int, 1}
    e_dims::Array{Int, 1}
    v_syms::Array{Symbol, 1}
    e_syms::Array{Symbol, 1}
    dim_v::Int
    dim_e::Int
    s_e::Array{Int, 1}
    d_e::Array{Int, 1}
    v_offs::Array{Int, 1}
    e_offs::Array{Int, 1}
    v_idx::Array{Idx, 1}
    e_idx::Array{Idx, 1}
    s_e_offs::Array{Int, 1}
    d_e_offs::Array{Int, 1}
    s_e_idx::Array{Idx, 1}
    d_e_idx::Array{Idx, 1}
    e_s_v_dat::Array{Array{Tuple{Int,Int}, 1}}
    e_d_v_dat::Array{Array{Tuple{Int,Int}, 1}}
end
function GraphStruct(g, v_dims, e_dims, v_syms, e_syms)
    num_v = nv(g)
    num_e = ne(g)

    s_e = [src(e) for e in edges(g)]
    d_e = [dst(e) for e in edges(g)]

    v_offs = create_offsets(v_dims)
    e_offs = create_offsets(e_dims)

    v_idx = create_idxs(v_offs, v_dims)
    e_idx = create_idxs(e_offs, e_dims)

    s_e_offs = [v_offs[s_e[i_e]] for i_e in 1:num_e]
    d_e_offs = [v_offs[d_e[i_e]] for i_e in 1:num_e]

    s_e_idx = [v_idx[s_e[i_e]] for i_e in 1:num_e]
    d_e_idx = [v_idx[d_e[i_e]] for i_e in 1:num_e]

    e_s_v_dat = [[(offset, dim) for (i_e, (offset, dim)) in enumerate(zip(e_offs, e_dims)) if i_v == s_e[i_e]] for i_v in 1:num_v]
    e_d_v_dat = [[(offset, dim) for (i_e, (offset, dim)) in enumerate(zip(e_offs, e_dims)) if i_v == d_e[i_e]] for i_v in 1:num_v]

    GraphStruct(
    num_v,
    num_e,
    v_dims,
    e_dims,
    v_syms,
    e_syms,
    sum(v_dims),
    sum(e_dims),
    s_e,
    d_e,
    v_offs,
    e_offs,
    v_idx,
    e_idx,
    s_e_offs,
    d_e_offs,
    s_e_idx,
    d_e_idx,
    e_s_v_dat,
    e_d_v_dat)
end

# In order to access the data in the arrays efficiently we create views that
# allow us to efficiently index into the underlying arrays.

import Base.getindex, Base.setindex!, Base.length, Base.IndexStyle, Base.size, Base.eltype, Base.dataids

#= Once forward declaration of types is implemented properly we can replace this
 by

incomplete struct GraphData{T}
end

struct EdgeData{T} <: AbstractArray{T, 1}
 gd::GraphData{T}
 idx_offset::Int
 len::Int
end

Then the recursive type signatures in GraphData that look like:

EdgeData{GraphData{T}, T}

will simplify to

EdgeData{T}

We should create a branch with this variant that can run on this branch of julia:

https://github.com/JuliaLang/julia/pull/32658
=#
struct EdgeData{G,T} <: AbstractArray{T, 1} # This is a lie. The elements of EdgeData are not of type T but of type eltype(T).
    gd::G
    idx_offset::Int
    len::Int
end

@inline Base.@propagate_inbounds function getindex(e_dat::EdgeData, idx)
    e_dat.gd.e_array[idx + e_dat.idx_offset]
end

@inline Base.@propagate_inbounds function setindex!(e_dat::EdgeData, x, idx)
    e_dat.gd.e_array[idx + e_dat.idx_offset] = x
    nothing
end

@inline function Base.length(e_dat::EdgeData)
    e_dat.len
end

@inline function Base.size(e_dat::EdgeData)
    (e_dat.len, )
end

@inline function Base.eltype(e_dat::EdgeData{G, T}) where {G, T}
    eltype(T)
end

Base.IndexStyle(::Type{<:EdgeData}) = IndexLinear()

@inline Base.dataids(e_dat::EdgeData) = dataids(e_dat.gd.e_array)

struct VertexData{G, T} <: AbstractArray{T, 1}
    gd::G
    idx_offset::Int
    len::Int
end

@inline Base.@propagate_inbounds function getindex(v_dat::VertexData, idx)
    v_dat.gd.v_array[idx + v_dat.idx_offset]
end

@inline Base.@propagate_inbounds function setindex!(v_dat::VertexData, x, idx)
    v_dat.gd.v_array[idx + v_dat.idx_offset] = x
    nothing
end

@inline function Base.length(v_dat::VertexData)
    v_dat.len
end

@inline function Base.size(e_dat::VertexData)
    (e_dat.len, )
end

@inline function Base.eltype(e_dat::VertexData{G, T}) where {G, T}
    eltype(T)
end

Base.IndexStyle(::Type{<:VertexData}) = IndexLinear()

@inline Base.dataids(v_dat::VertexData) = dataids(v_dat.gd.v_array)

# Putting the above together we create a GraphData object:

# An alternative design that needs to be evaluated for performance is to create
# only one array of VertexData and EdgeData and index into that, possibly with a
# new set of access types...

# We require potentially different data types for vertices and edges because
# there are situations with autodifferentiation that require one of them to be
# dual and the other not.

mutable struct GraphData{Tv, Te}
    v_array::Tv
    e_array::Te
    v::Array{VertexData{GraphData{Tv, Te}, Tv}, 1}
    e::Array{EdgeData{GraphData{Tv, Te}, Te}, 1}
    v_s_e::Array{VertexData{GraphData{Tv, Te}, Tv}, 1} # the vertex that is the source of e
    v_d_e::Array{VertexData{GraphData{Tv, Te}, Tv}, 1} # the vertex that is the destination of e
    e_s_v::Array{Array{EdgeData{GraphData{Tv, Te}, Te}, 1}, 1} # the edges that have v as source
    e_d_v::Array{Array{EdgeData{GraphData{Tv, Te}, Te}, 1}, 1} # the edges that have v as destination
    function GraphData{Tv, Te}(v_array::Tv, e_array::Te, gs::GraphStruct) where {Tv, Te}
        gd = new{Tv, Te}(v_array, e_array, )
        gd.v = [VertexData{GraphData{Tv, Te}, Tv}(gd, offset, dim) for (offset,dim) in zip(gs.v_offs, gs.v_dims)]
        gd.e = [EdgeData{GraphData{Tv, Te}, Te}(gd, offset, dim) for (offset,dim) in zip(gs.e_offs, gs.e_dims)]
        gd.v_s_e = [VertexData{GraphData{Tv, Te}, Tv}(gd, offset, dim) for (offset,dim) in zip(gs.s_e_offs, gs.v_dims[gs.s_e])]
        gd.v_d_e = [VertexData{GraphData{Tv, Te}, Tv}(gd, offset, dim) for (offset,dim) in zip(gs.d_e_offs, gs.v_dims[gs.d_e])]
        gd.e_s_v = [[EdgeData{GraphData{Tv, Te}, Te}(gd, offset, dim) for (offset,dim) in e_s_v] for e_s_v in gs.e_s_v_dat]
        gd.e_d_v = [[EdgeData{GraphData{Tv, Te}, Te}(gd, offset, dim) for (offset,dim) in e_d_v] for e_d_v in gs.e_d_v_dat]
        gd
    end
end

function GraphData(v_array, e_array, gs)
    GraphData{typeof(v_array), typeof(e_array)}(v_array, e_array, gs)
end

#= In order to manipulate initial conditions using this view of the underlying
array we provide view functions that give access to the arrays. =#

export view_v
export view_e

function view_v(gd::GraphData, gs::GraphStruct, sym="")
    v_idx = [i for (i, s) in enumerate(gs.v_syms) if occursin(string(sym), string(s))]
    view(gd.v_array, v_idx)
end

function view_e(gd::GraphData, gs::GraphStruct, sym="")
    e_idx = [i for (i, s) in enumerate(gs.e_syms) if occursin(string(sym), string(s))]
    view(gd.e_array, e_idx)
end


function view_v(nd, x, p, t, sym="")
    gd = nd(x, p, t, GetGD)
    gs = nd(GetGS)
    v_idx = [i for (i, s) in enumerate(gs.v_syms) if occursin(string(sym), string(s))]
    view(gd.v_array, v_idx)
end

function view_e(nd, x, p, t, sym="")
    gd = nd(x, p, t, GetGD)
    gs = nd(GetGS)
    e_idx = [i for (i, s) in enumerate(gs.e_syms) if occursin(string(sym), string(s))]
    view(gd.e_array, e_idx)
end

export swap_v_array!, swap_e_array!

"""
    swap_v_array!(gd:GraphData, array)

Swaps the underlying vertex data array of an GraphData type with a new one.
"""
@inline function swap_v_array!(gd::GraphData{Tv, Te}, array::Tv) where {Tv, Te}
    gd.v_array = array
end

"""
    swap_e_array!(gd:GraphData, array)

Swaps the underlying edge data array of an GraphData type with a new one.
"""
@inline function swap_e_array!(gd::GraphData{Tv, Te}, array::Tv) where {Tv, Te}
    gd.e_array = array
end

export get_vertex, get_edge, get_src_vertex, get_dst_vertex, get_out_edges, get_in_edges

"""
    get_vertex(gd::GraphData, idx::Int) -> View

Returns a view-like access to the underlying data of the i-th vertex.
"""
@inline get_vertex(gd::GraphData, i::Int) = gd.v[i]

"""
    get_edge(gd::GraphData, idx::Int) -> View

Returns a view-like access to the underlying data of the i-th edge.
"""
@inline get_edge(gd::GraphData, i::Int) = gd.e[i]

"""
    get_src_vertex(gd::GraphData, idx::Int) -> View

Returns a view-like access to the underlying data of source vertex of the i-th edge.
"""
@inline get_src_vertex(gd::GraphData, i::Int) = gd.v_s_e[i]

"""
    get_dst_vertex(gd::GraphData, idx::Int) -> View

Returns a view-like access to the underlying data of destination vertex of the i-th edge.
"""
@inline get_dst_vertex(gd::GraphData, i::Int) = gd.v_d_e[i]

"""
    get_out_edges(gd::GraphData, i::Int)

Returns an Vector of view-like accesses to all the outgoing edges of the i-th vertex.
"""
@inline get_out_edges(gd::GraphData, i::Int) = gd.e_s_v[i]

"""
    get_in_edges(gd::GraphData, i::Int)

Returns an Vector of view-like accesses to all the incoming edges of the i-th vertex.
"""
@inline get_in_edges(gd::GraphData, i) = gd.e_d_v[i]

function construct_mass_matrix(mmv_array, gs)
    if all([mm == I for mm in mmv_array])
        mass_matrix = I
    else
        mass_matrix = sparse(1.0I,gs.dim_v, gs.dim_v)
        for (i, mm) in enumerate(mmv_array)
            ind = gs.v_idx[i]
            if ndims(mm) == 0
                copyto!(mass_matrix[ind, ind], mm*I)
            elseif ndims(mm) == 1
                copyto!(mass_matrix[ind, ind], Diagonal(mm))
            elseif ndims(mm) == 2 # ndims(I) = 2
                # `I` does not support broadcasting but copyto!
                copyto!(mass_matrix[ind, ind], mm)
            else
                error("The mass matrix needs to be interpretable as a 2D matrix.")
            end
        end
    end
    mass_matrix
end

function construct_mass_matrix(mmv_array, mme_array, gs)
    if all([mm == I for mm in mmv_array]) && all([mm == I for mm in mme_array])
        mass_matrix = I
    else
        dim_nd = gs.dim_v + gs.dim_e
        mass_matrix = sparse(1.0I,dim_nd,dim_nd)
        for (i, mm) in enumerate(mmv_array)
            ind = gs.v_idx[i]
            if ndims(mm) == 0
                copyto!(mass_matrix[ind, ind], mm*I)
            elseif ndims(mm) == 1
                copyto!(mass_matrix[ind, ind], Diagonal(mm))
            elseif ndims(mm) == 2 # ndims(I) = 2
                # `I` does not support broadcasting but copyto!
                copyto!(mass_matrix[ind, ind], mm)
            else
                error("The mass matrix needs to be interpretable as a 2D matrix.")
            end
        end
        for (i, mm) in enumerate(mme_array)
            ind = gs.e_idx[i]
            if ndims(mm) == 0
                copyto!(mass_matrix[ind, ind], mm*I)
            elseif ndims(mm) == 1
                copyto!(mass_matrix[ind, ind], Diagonal(mm))
            elseif ndims(mm) == 2 # ndims(I) = 2
                # `I` does not support broadcasting but copyto!
                copyto!(mass_matrix[ind, ind], mm)
            else
                error("The mass matrix needs to be interpretable as a 2D matrix.")
            end
        end
    end
    mass_matrix
end


#= These types are used to dispatch the network dynamics functions to provide
access to the underlying GraphData and GraphStruct objects. =#
export GetGD

struct GetGD
end

export GetGS

struct GetGS
end


#= Experimental and untested: Wrap a solution object so we get back a GraphData
object at every time. =#
export ND_Solution

struct ND_Solution
    nd
    p
    sol
end
function (nds::ND_Solution)(t)
    nds.nd(nds.sol(t), nds.p, t, GetGD)
end

end # module
