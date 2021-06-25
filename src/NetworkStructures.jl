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
export GraphStruct, GraphData, EdgeData, VertexData

const Idx = UnitRange{Int}

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

    num_v::Int                                 # number of vertices
    num_e::Int                                 # number of edges

    v_dims::Array{Int, 1}                      # dimensions per vertex
    e_dims::Array{Int, 1}                      # dimensions per edge

    v_syms::Array{Symbol, 1}                   # symbol per vertex
    e_syms::Array{Symbol, 1}                   # symbol per edge

    dim_v::Int                                 # total vertex dimensions
    dim_e::Int                                 # total edge dimensions

    s_e::Array{Int, 1}                         # src-vertex idx per edge
    d_e::Array{Int, 1}                         # dst-vertex idx per edge

    s_v::Array{Array{Int,1}}                   # indices of source edges per vertex
    d_v::Array{Array{Int,1}}                   # indices of destination edges per vertex

    v_offs::Array{Int, 1}                      # linear offset per vertex
    e_offs::Array{Int, 1}                      # linear offset per edge

    v_idx::Array{Idx, 1}                       # lin. idx-range per vertex
    e_idx::Array{Idx, 1}                       # lin. idx-range per edge

    s_e_offs::Array{Int, 1}                    # offset of src-vertex per edge
    d_e_offs::Array{Int, 1}                    # offset of dst-vertex per edge

    s_e_idx::Array{Idx, 1}                     # idx-range of src-vertex per edge
    d_e_idx::Array{Idx, 1}                     # idx-range of dst-vertex per edge

    dst_edges_dat::Vector{Vector{Tuple{Int,Int}}}
end
function GraphStruct(g, v_dims, e_dims, v_syms, e_syms)
    num_v = nv(g)
    num_e = ne(g)

    s_e = src.(edges(g))
    d_e = dst.(edges(g))

    # Initialize Vector of empty Vector{Int} for the edge indices per vertex
    s_v = [Int[] for i=1:num_v]
    d_v = [Int[] for i=1:num_v]
    for i = 1:num_e
        push!(s_v[s_e[i]], i)
        push!(d_v[d_e[i]], i)
    end


    v_offs = create_offsets(v_dims)
    e_offs = create_offsets(e_dims)

    v_idx = create_idxs(v_offs, v_dims)
    e_idx = create_idxs(e_offs, e_dims)

    s_e_offs = [v_offs[s_e[i_e]] for i_e in 1:num_e]
    d_e_offs = [v_offs[d_e[i_e]] for i_e in 1:num_e]

    s_e_idx = [v_idx[s_e[i_e]] for i_e in 1:num_e]
    d_e_idx = [v_idx[d_e[i_e]] for i_e in 1:num_e]

    dst_edges_dat = Vector{Vector{Tuple{Int,Int}}}(undef, nv(g))

    for i_v in 1:nv(g)
        offsdim_arr = Tuple{Int,Int}[]
        for i_e in d_v[i_v]
            # dims is a multiple of 2 for SimpleGraph by design of VertexFunction
            if !is_directed(g)
                push!(offsdim_arr, (e_offs[i_e], e_dims[i_e] / 2))
            else
                push!(offsdim_arr, (e_offs[i_e], e_dims[i_e]))
            end
        end
        for i_e in s_v[i_v]
            if !is_directed(g)
                push!(offsdim_arr, (e_offs[i_e] +  e_dims[i_e] / 2, e_dims[i_e] / 2))
            end
        end
        dst_edges_dat[i_v] = offsdim_arr
    end

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
    s_v,
    d_v,
    v_offs,
    e_offs,
    v_idx,
    e_idx,
    s_e_offs,
    d_e_offs,
    s_e_idx,
    d_e_idx,
    dst_edges_dat)
end


import Base.getindex, Base.setindex!, Base.length, Base.IndexStyle, Base.size, Base.eltype, Base.dataids

"""
    struct EdgeData{GDB, elE} <: AbsstractArray{elE, 1}

The EdgeData object behaves like an array and allows access to the underlying
data of a specific edge (like a View). Unlike a View, the parent array is stored
in a mutable GraphDataBuffer object and can be swapped.
"""
struct EdgeData{GDB, elE} <: AbstractArray{elE, 1}
    gdb::GDB
    idx_offset::Int
    len::Int
end

@inline Base.@propagate_inbounds function getindex(e_dat::EdgeData, idx)
    e_dat.gdb.e_array[idx + e_dat.idx_offset]
end

@inline Base.@propagate_inbounds function setindex!(e_dat::EdgeData, x, idx)
    e_dat.gdb.e_array[idx + e_dat.idx_offset] = x
    nothing
end

@inline function Base.length(e_dat::EdgeData)
    e_dat.len
end

@inline function Base.size(e_dat::EdgeData)
    (e_dat.len, )
end

@inline function Base.eltype(e_dat::EdgeData{GDB, elE}) where {GDB, elE}
    elE
end

Base.IndexStyle(::Type{<:EdgeData}) = IndexLinear()

@inline Base.dataids(e_dat::EdgeData) = dataids(e_dat.gdb.e_array)

"""
    struct VertexData{GDB, elV} <: AbsstractArray{elV, 1}

The VertexData object behaves like an array and allows access to the underlying
data of a specific vertex (like a View). Unlike a View, the parent array is stored
in a mutable GraphDataBuffer object and can be swapped.
"""
struct VertexData{GDB, elV} <: AbstractArray{elV, 1}
    gdb::GDB
    idx_offset::Int
    len::Int
end

@inline Base.@propagate_inbounds function getindex(v_dat::VertexData, idx)
    v_dat.gdb.v_array[idx + v_dat.idx_offset]
end

@inline Base.@propagate_inbounds function setindex!(v_dat::VertexData, x, idx)
    v_dat.gdb.v_array[idx + v_dat.idx_offset] = x
    nothing
end

@inline function Base.length(v_dat::VertexData)
    v_dat.len
end

@inline function Base.size(e_dat::VertexData)
    (e_dat.len, )
end

@inline function Base.eltype(e_dat::VertexData{G, elV}) where {G, elV}
    elV
end

Base.IndexStyle(::Type{<:VertexData}) = IndexLinear()

@inline Base.dataids(v_dat::VertexData) = dataids(v_dat.gdb.v_array)


"""
    mutable struct GraphDataBuffer{Tv, Te}

Is a composite type which holds two Arrays for the underlying data of a graph.
The type is mutable, therfore the v_array and e_array can be changed.
"""
mutable struct GraphDataBuffer{Tv, Te}
    v_array::Tv
    e_array::Te
end

"""
    GraphData{GDB, elV, elE}

The GraphData object contains a reference to the GraphDataBuffer object and to all the
view-like EdgeData/VertexData objects. It is used to access the underlying linear data
of a graph in terms of edges and vertices. The underlying data kann be swapped using the
    swap_v_array
    swap_e_array
methods.
The data for specific edges/vertices can be accessed using the
    get_vertex, get_edge
    get_src_vertex, get_dst_vertex
    get_src_edges, get_dst_edges
methods.
"""
struct GraphData{GDB, elV, elE}
    gdb::GDB
    v::Array{VertexData{GDB, elV}, 1}
    e::Array{EdgeData{GDB, elE}, 1}
    v_s_e::Array{VertexData{GDB, elV}, 1} # the vertex that is the source of e
    v_d_e::Array{VertexData{GDB, elV}, 1} # the vertex that is the destination of e
    dst_edges::Array{Array{EdgeData{GDB, elE}, 1}, 1} # the half-edges that have v as destination
end

function GraphData(v_array::Tv, e_array::Te, gs::GraphStruct; global_offset = 0) where {Tv, Te}
    gdb = GraphDataBuffer{Tv, Te}(v_array, e_array)
    GDB = typeof(gdb)
    elV = eltype(v_array)
    elE = eltype(e_array)
    v = [VertexData{GDB, elV}(gdb, offset + global_offset, dim) for (offset,dim) in zip(gs.v_offs, gs.v_dims)]
    e = [EdgeData{GDB, elE}(gdb, offset + global_offset, dim) for (offset,dim) in zip(gs.e_offs, gs.e_dims)]
    v_s_e = [VertexData{GDB, elV}(gdb, offset + global_offset, dim) for (offset,dim) in zip(gs.s_e_offs, gs.v_dims[gs.s_e])]
    v_d_e = [VertexData{GDB, elV}(gdb, offset + global_offset, dim) for (offset,dim) in zip(gs.d_e_offs, gs.v_dims[gs.d_e])]
    dst_edges = [[EdgeData{GDB, elE}(gdb, offset + global_offset, dim) for (offset,dim) in in_edge] for in_edge in gs.dst_edges_dat]
    GraphData{GDB, elV, elE}(gdb, v, e, v_s_e, v_d_e, dst_edges)
end


#= In order to manipulate initial conditions using this view of the underlying
array we provide view functions that give access to the arrays. =#

export view_v
export view_e

function view_v(gd::GraphData, gs::GraphStruct, sym="")
    v_idx = [i for (i, s) in enumerate(gs.v_syms) if occursin(string(sym), string(s))]
    view(gd.gdb.v_array, v_idx)
end

function view_e(gd::GraphData, gs::GraphStruct, sym="")
    e_idx = [i for (i, s) in enumerate(gs.e_syms) if occursin(string(sym), string(s))]
    view(gd.gdb.e_array, e_idx)
end


function view_v(nd, x, p, t, sym="")
    gd = nd(x, p, t, GetGD)
    gs = nd(GetGS)
    v_idx = [i for (i, s) in enumerate(gs.v_syms) if occursin(string(sym), string(s))]
    view(gd.gdb.v_array, v_idx)
end

function view_e(nd, x, p, t, sym="")
    gd = nd(x, p, t, GetGD)
    gs = nd(GetGS)
    e_idx = [i for (i, s) in enumerate(gs.e_syms) if occursin(string(sym), string(s))]
    view(gd.gdb.e_array, e_idx)
end

export swap_v_array!, swap_e_array!

"""
    swap_v_array!(gd:GraphData, array)

Swaps the underlying vertex data array of an GraphData type with a new one.
"""
@inline function swap_v_array!(gd::GraphData{GDB, elV, elE}, array::AbstractArray{elV}) where {GDB, elV, elE}
    gd.gdb.v_array = array
end

"""
    swap_e_array!(gd:GraphData, array)

Swaps the underlying edge data array of an GraphData type with a new one.
"""
@inline function swap_e_array!(gd::GraphData{GDB, elV, elE}, array::AbstractArray{elE}) where {GDB, elV, elE}
    gd.gdb.e_array = array
end

export get_vertex, get_edge, get_src_vertex, get_dst_vertex, get_src_edges, get_dst_edges

"""
    get_vertex(gd::GraphData, idx::Int) -> View

Returns a view-like access to the underlying data of the i-th vertex.
"""
@inline @Base.propagate_inbounds get_vertex(gd::GraphData, i::Int) = gd.v[i]

"""
    get_edge(gd::GraphData, idx::Int) -> View

Returns a view-like access to the underlying data of the i-th edge.
"""
@inline @Base.propagate_inbounds get_edge(gd::GraphData, i::Int) = gd.e[i]

"""
    get_src_vertex(gd::GraphData, idx::Int) -> View

Returns a view-like access to the underlying data of source vertex of the i-th edge.
"""
@inline @Base.propagate_inbounds get_src_vertex(gd::GraphData, i::Int) = gd.v_s_e[i]

"""
    get_dst_vertex(gd::GraphData, idx::Int) -> View

Returns a view-like access to the underlying data of destination vertex of the i-th edge.
"""
@inline @Base.propagate_inbounds get_dst_vertex(gd::GraphData, i::Int) = gd.v_d_e[i]

"""
    get_src_edges(gd::GraphData, i::Int)

Returns a Vector of view-like accesses to all the (half-)edges that have the i-th vertex as source (for directed graphs these are the out-edges).
"""
@inline @Base.propagate_inbounds get_src_edges(gd::GraphData, i::Int) = error("Not implemented.")

"""
    get_dst_edges(gd::GraphData, i::Int)

Returns a Vector of view-like accesses to all the (half-)edges that have the i-th vertex as destination (for directed graphs these are the in-edges).
"""
@inline @Base.propagate_inbounds get_dst_edges(gd::GraphData, i) = gd.dst_edges[i]




export GetGD
export GetGS

"""This type is used to dispatch the network dynamics functions to provide
access to the underlying GraphData object."""
struct GetGD
end

"""This type is used to dispatch the network dynamics functions to provide
access to the underlying GraphStruct object."""
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
