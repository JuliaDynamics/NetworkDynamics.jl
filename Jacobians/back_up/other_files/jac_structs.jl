"""
Together with Constructors this module forms the backbone of the core API.
It provide the basic types to construct Arrays of VertexFunction and
EdgeFunction which can be handled by network_dynamics.
"""
#module ComponentFunctions
module jac_structs

using LinearAlgebra

import Base.convert
import Base.promote_rule


export VertexFunction
export EdgeFunction
export StaticVertex
export StaticEdge
export ODEVertex
export ODEEdge
export DDEVertex
export StaticDelayEdge

"""
Abstract supertype for all vertex functions.
"""
abstract type VertexFunction end
"""
Abstract supertype for all edge functions.
"""
abstract type EdgeFunction end

"""
    StaticVertex(f!, dim, sym)

Wrapper that ensures compatibility of a **mutating** function **`f!`** with
the key constructor `network_dynamics`.

**`f!`**  describes the local behaviour at a static node and has to respect
the following calling syntax

```julia
f!(v, edges, p, t) -> nothing
```

Here  `v`, `p` and `t` are the usual arguments, while
`edges` is an arrays containing the edges for which the
described vertex the destination (in-edges for directed graphs).

`dim` is the number of independent variables in the vertex equations and
`sym` is an array of symbols for these variables.

For more details see the documentation.
"""
@Base.kwdef struct StaticVertex{T} <: VertexFunction
    f!::T
    dim::Int
    sym=[:v for i in 1:dim]
    vertex_jacobian! = :F # signature (J::AbstractMatrix, v, p, t) -> nothing
    # vertex_jacobian! ist optionales fieldname, was als Funktion vom user angegeben wird und bei
    # Objekterstellung von Vertex-Funktion angegeben werden kann
end

"""
    StaticEdge(f!, dim, sym)

Wrapper that ensures compatibility of a **mutating** function `f!` with
the key constructor `network_dynamics`.

`f!`  describes the local behaviour at a static edge and has to respect
the following calling syntax

```julia
f!(e, v_s, v_t, p, t) -> nothing
```

Here  `e`, `p` and `t` are the usual arguments, while
`v_s` and `v_d` are arrays containing the vertices which are
the source and destination of the described edge.

- `dim` is the number of independent variables in the edge equations and
- `sym` is an array of symbols for these variables.
- `coupling` is a Symbol describing if the EdgeFunction is intended for a directed graph (`:directed`) or for an undirected graph (`{:undirected, :symmetric, :antisymmetric, :fiducial}`). `:directed` is intended for directed graphs. `:undirected` is the default option and is only compatible with SimpleGraph. in this case f! should specify the coupling from a source vertex to a destination vertex. `:symmetric` and `:antisymmetric` trigger performance optimizations, if `f!` has that symmetry property. `:fiducial` lets the user specify both the coupling from src to dst, as well as the coupling from dst to src and is intended for advanced users.

For more details see the documentation.
"""
@Base.kwdef struct StaticEdge{T} <: EdgeFunction
    f!::T # (e, v_s, v_t, p, t) -> nothing
    dim::Int # number of dimensions of e
    coupling = :undefined # :directed, :symmetric, :antisymmetric, :fiducial, :undirected
    sym=[:e for i in 1:dim] # Symbols for the dimensions
    edge_jacobian! = :F # signature (J_s::AbstractMatrix, J_d::AbstractMatric, v, p, t) -> nothing


    function StaticEdge(user_f!::T,
                           dim::Int,
                           coupling::Symbol,
                           sym::Vector{Symbol},
                           edge_jacobian!) where T

        coupling_types = (:undefined, :directed, :fiducial, :undirected, :symmetric,
                          :antisymmetric)

        coupling ∈ coupling_types ? nothing :
            error("Coupling type not recognized. Choose from $coupling_types.")

        dim > 0 ? nothing : error("dim has to be a positive number.")

        dim == length(sym) ? nothing : error("Please specify a symbol for every dimension.")

        if coupling ∈ [:undefined, :directed]
            return new{T}(user_f!, dim, coupling, sym)

        elseif coupling == :fiducial
            dim % 2 == 0 ? nothing : error("Fiducial edges are required to have even dim.
                                            The first dim args are used for src -> dst,
                                            the second for dst -> src coupling.")
            return new{T}(user_f!, dim, coupling, sym)

        elseif coupling == :undirected
            # This might cause unexpected behaviour if source and destination vertex don't
            # have the same internal arguments.
            # Make sure to explicitly define the edge is :fiducial in that case.
            f! = @inline (e, v_s, v_d, p, t) -> begin
                @inbounds user_f!(view(e,1:dim), v_s, v_d, p, t)
                @inbounds user_f!(view(e,dim+1:2dim), v_d, v_s, p, t)
                nothing
            end
        elseif coupling == :antisymmetric
            f! = @inline (e, v_s, v_d, p, t) -> begin
                @inbounds user_f!(view(e,1:dim), v_s, v_d, p, t)
                @inbounds for i in 1:dim
                    e[dim + i] = -1.0 * e[i]
                end
                nothing
            end
        elseif coupling == :symmetric
            f! = @inline (e, v_s, v_d, p, t) -> begin
                @inbounds user_f!(view(e,1:dim), v_s, v_d, p, t)
                @inbounds for i in 1:dim
                    e[dim + i] = e[i]
                end
                nothing
            end
        end
        # For edges with mass matrix this will be a little more complicated
        return new{typeof(f!)}(f!, 2dim, coupling, repeat(sym, 2))
    end
end




"""
    ODEVertex(f!, dim, mass_matrix, sym)

Wrapper that ensures compatibility of a **mutating** function **`f!`** with
the key constructor `network_dynamics`.

**`f!`**  describes the local behaviour at a dynamic node and has to respect
the following calling syntax

```julia
f!(dv, v, edges, p, t) -> nothing
```

Here `dv`, `v`, `p` and `t` are the usual ODE arguments, while
`edges` is an Array containing the edges for which the vertex is the destination (in-edges for directed graphs).

**`dim`** is the number of independent variables in the vertex equations and
**`sym`** is an array of symbols for these variables.
**`mass_matrix`** is an optional argument that defaults to the identity
matrix `I`. If a mass matrix M is given the system `M * dv = f!` will be
solved.

For more details see the documentation.
"""
@Base.kwdef struct ODEVertex{T} <: VertexFunction
    f!::T # signature (dx, x, edges, p, t) -> nothing
    dim::Int
    mass_matrix=I
    sym=[:v for i in 1:dim]
    vertex_jacobian! = :F # signature (J::AbstractMatrix, v, p, t)
end

"""
    ODEEdge(f!, dim, mass_matrix, sym)

Wrapper that ensures compatibility of a **mutating** function **`f!`** with
the key constructor `network_dynamics`.

**`f!`**  describes the local behaviour at a dynamic edge and has to respect
the following calling syntax

```julia
f!(de, e, v_s, v_t, p, t) -> nothing
```

Here  `de`, `e`, `p` and `t` are the usual arguments, while
`v_s` and `v_d` are arrays containing the vertices which are
the source and destination of the described edge.

**`dim`** is the number of independent variables in the edge equations and
**`sym`** is an array of symbols for these variables. For more details see
the documentation.
**`mass_matrix`** is an optional argument that defaults to the identity
matrix `I`. If a mass matrix M is given the system `M * de = f!` will be
solved.

For more details see the documentation.
"""
@Base.kwdef struct ODEEdge{T} <: EdgeFunction
    f!::T # (de, e, v_s, v_t, p, t) -> nothing
    dim::Int # number of dimensions of e
    coupling = :undefined # :directed, :symmetric, :antisymmetric, :fiducial, :undirected
    mass_matrix=I # Mass matrix for the equation
    sym=[:e for i in 1:dim] # Symbols for the dimensions
    edge_jacobian! = :F # signature (J_s::AbstractMatrix, J_d::AbstractMatric, v, p, t) -> nothing


    function ODEEdge(user_f!::T,
                     dim::Int,
                     coupling::Symbol,
                     mass_matrix,
                     sym::Vector{Symbol},
                     edge_jacobian!) where T

        coupling_types = (:directed, :fiducial, :undirected)

        coupling == :undefined ?
            error("ODEEdges with undefined coupling type are not implemented at the "*
            "moment. Choose `coupling` from $coupling_types.") : nothing

        coupling ∈ (:symmetric, :antisymmetric) ?
             error("Coupling type $coupling is not available for ODEEdges.") : nothing



        coupling ∈ coupling_types ? nothing :
            error("Coupling type not recognized. Choose from $coupling_types.")


        dim > 0 ? nothing : error("dim has to be a positive number.")

        dim == length(sym) ? nothing : error("Please specify a symbol for every dimension.")

        if coupling == :directed
            return new{T}(user_f!, dim, coupling, mass_matrix, sym)
        elseif coupling == :fiducial
            dim % 2 == 0 ? nothing : error("Fiducial edges are required to have even dim.
                                            The first dim args are used for src -> dst,
                                            the second for dst -> src coupling.")
            return new{T}(user_f!, dim, coupling, mass_matrix, sym)

        elseif coupling == :undirected
            # This might cause unexpected behaviour if source and destination vertex don't
            # have the same internal arguments.
            # Make sure to explicitly define the edge is :fiducial in that case.
            f! = @inline (de, e, v_s, v_d, p, t) -> begin
                @inbounds user_f!(view(de,1:dim), view(e,1:dim), v_s, v_d, p, t)
                @inbounds user_f!(view(de,dim+1:2dim), view(e,dim+1:2dim), v_d, v_s, p, t)
                nothing
            end
            let M = mass_matrix
                if M === I
                    newM = M
                elseif M isa Number
                    newM = M
                elseif M isa Vector
                    newM = repeat(M,2)
                elseif M isa Matrix
                    newM = [M zeros(size(M)); zeros(size(M)) M]
                end
                return new{typeof(f!)}(f!, 2dim, coupling, newM, repeat(sym, 2))
            end
        end
    end
end

######################## Graph Data Block ######################################
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

    s_e = [src(e) for e in edges(g)]
    d_e = [dst(e) for e in edges(g)]

    s_v = [findall(isequal(i), src.(edges(g))) for i = 1:nv(g)]
    d_v = [findall(isequal(i), dst.(edges(g))) for i = 1:nv(g)]

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

################################################################################

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

################################################################################

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

################################################################################

mutable struct GraphDataBuffer{Tv, Te}
    v_array::Tv
    e_array::Te
end

mutable struct JacGraphDataBuffer{Tvj, Tej, Tep}
    v_Jac_array::Tvj
    e_Jac_array::Tej
    e_Jac_product_array::Tep
end

struct GraphData{GDB, elV, elE}
    gdb::GDB
    v::Array{VertexData{GDB, elV}, 1}
    e::Array{EdgeData{GDB, elE}, 1}
    v_s_e::Array{VertexData{GDB, elV}, 1} # the vertex that is the source of e
    v_d_e::Array{VertexData{GDB, elV}, 1} # the vertex that is the destination of e
    dst_edges::Array{Array{EdgeData{GDB, elE}, 1}, 1} # the half-edges that have v as destination
end

struct JacGraphData{JGDB}
    jgdb::JGDB
    v_jac_array::Array{Array{Float64, 2}, 1}
    e_jac_array::Array{Array{Array{Float64, 2}, 1}, 1}
    e_jac_product::Array{Float64, 2}
end

#=struct JacGraphData{JGDB, elV, elE, elEP}
#    jgdb::JGDB
#    v_jac_array::Array{VertexData{JGDB, elV}, 1}
#    e_jac_array::Array{Array{EdgeData{JGDB, elE}, 1}, 1}
#    e_jac_product::EdgeData{JGDB, elEP}
end=#

function GraphData(v_array::Tv, e_array::Te, gs::GraphStruct; global_offset = 0) where {Tv, Te}
    gdb = GraphDataBuffer{Tv, Te}(v_array, e_array)
    GDB = typeof(gdb)
    elV = eltype(v_array)
    elE = eltype(e_array)
    v = [VertexData{GDB, elV}(gdb, offset + global_offset, dim) for (offset,dim) in zip(gs.v_offs, gs.v_dims)]
    e = [EdgeData{GDB, elE}(gdb, offset + global_offset, dim) for (offset,dim) in zip(gs.e_offs, gs.e_dims)]
    v_s_e = [VertexData{GDB, elV}(gdb, offset + global_offset, dim) for (offset,dim) in zip(gs.s_e_offs, gs.v_dims[gs.s_e])] # source vertex of ith edge
    v_d_e = [VertexData{GDB, elV}(gdb, offset + global_offset, dim) for (offset,dim) in zip(gs.d_e_offs, gs.v_dims[gs.d_e])]
    dst_edges = [[EdgeData{GDB, elE}(gdb, offset + global_offset, dim) for (offset,dim) in in_edge] for in_edge in gs.dst_edges_dat]
    GraphData{GDB, elV, elE}(gdb, v, e, v_s_e, v_d_e, dst_edges)
end

function JacGraphData(v_Jac_array::Tvj, e_Jac_array::Tej, e_Jac_product_array::Tep, gs::GraphStruct) where {Tvj, Tej, Tep}
    jgdb = JacGraphDataBuffer{Tvj, Tej, Tep}(v_Jac_array, e_Jac_array, e_Jac_product_array)

    v_jac = [Array{Float64,2}(undef, dim, dim) for dim in gs.v_dims]
    e_jac = [[zeros(dim, srcdim), zeros(dim, dstdim)] for (dim, srcdim, dstdim) in zip(gs.e_dims, gs.v_dims, gs.v_dims)] # homogene Netzwerke: v_src_dim = v_dst_dim = v_dim
    e_jac_product =  zeros(gs.e_dims[1], gs.num_e) # Annahme: homogene edges

    JacGraphData{JGDB}(jgdb, v_jac, e_jac, e_jac_product)
end

#=function JacGraphData(v_Jac_array::Tvj, e_Jac_array::Tej, e_Jac_product_array::Tep, gs::GraphStruct) where {Tvj, Tej}
    jgdb = JacGraphDataBuffer{Tvj, Tej, Tep}(v_Jac_array, e_Jac_array, e_Jac_product_array)
    JGDB = typeof(jgdb)
    elV = eltype(v_jac_array)
    elE = eltype(e_Jac_array)
    elEP = eltype(e_Jac_product_array)

    v_jac_array = [Array{Float64,2}(undef, dim, dim) for dim in gs.v_dims]
    e_jac_array = [[zeros(dim, srcdim), zeros(dim, dstdim)] for (dim, srcdim, dstdim) in zip(gs.e_dims, gs.v_dims, gs.v_dims)] # homogene Netzwerke: v_src_dim = v_dst_dim = v_dim
    e_jac_product =  zeros(gs.e_dims[1], gs.num_e) # Annahme: homogene edges

    JacGraphData{JGDB}(jgdb, v_jac_array, e_jac_array, e_jac_product)
end =#

### TO DO: in NetworkStructures.jl

@inline get_src_edge_jacobian(jgd::JacGraphData, i::Int) = jgd.e_jac_array[i][1]

@inline get_dst_edge_jacobian(jgd::JacGraphData, i::Int) = jgd.e_jac_array[i][2]

@inline get_vertex_jacobian(jgd::JacGraphData, i::Int) = jgd.v_jac_array[i]


############################ NDJacVecOperator ##################################
using DifferentialEquations

mutable struct NDJacVecOperator{uType, tType, G, GD, JGD, T} <: DiffEqBase.AbstractDiffEqLinearOperator{T} # mutable da x, p, t geupdated werden
    x::uType
    p
    t::tType
    graph::G
    graph_structure::GraphStruct
    graph_data::GD
    jac_graph_data::JGD # - dann muss oben noch JGD hin
    parallel::Bool
end


############################ NetworkDynamics Structs ###########################

"""
    network_dynamics(vertices!, edges!, g; parallel = false, jac = false)

Assembles the the dynamical equations of the network problem into an `ODEFunction`
compatible with the `DifferentialEquations.jl` solvers. Takes as arguments an array
of VertexFunctions **`vertices!`**, an array of EdgeFunctions **`edges!`** and a
`LightGraph.jl` object **`g`**. The optional argument `parallel` is a boolean
value that denotes if the central loop should be executed in parallel with the number of threads set by the environment variable `JULIA_NUM_THREADS`.
"""
function network_dynamics(vertices!::Union{Array{T, 1}, T},
                          edges!::Union{Array{U, 1}, U},
                          graph;
                          x_prototype=zeros(1),
                          parallel=false,
                          jac = false) where {T <: ODEVertex, U <: StaticEdge}

    warn_parallel(parallel)

    # user_edges! = copy(edges!)
    edges! = prepare_edges(edges!, graph)


    v_dims, e_dims, symbols_v, symbols_e, mmv_array, mme_array = collect_ve_info(vertices!, edges!, graph)

    # These arrays are used for initializing the GraphData and will be overwritten
    v_array = similar(x_prototype, sum(v_dims))
    e_array = similar(x_prototype, sum(e_dims))

    symbols = symbols_v

    graph_stucture = GraphStruct(graph, v_dims, e_dims, symbols_v, symbols_e) # Funktion

    graph_data = GraphData(v_array, e_array, graph_stucture) # Funktion

    nd! = nd_ODE_Static(vertices!, edges!, graph, graph_stucture, graph_data, parallel) # Objekterstellung, nd_ODE_Static hier struct

    mass_matrix = construct_mass_matrix(mmv_array, graph_stucture)

    if jac == true
        # These additional arrays are used for initializing the GraphData and will be overwritten
        v_jac_array = [Array{Float64,2}(undef, dim, dim) for dim in v_dims]
        e_jac_array = [[zeros(dim, srcdim), zeros(dim, dstdim)] for (dim, srcdim, dstdim) in zip(e_dims, v_dims, v_dims)] # homogene Netzwerke: v_src_dim = v_dst_dim = v_dim
        e_jac_product =  zeros(e_dims[1], graph_stucture.num_e) # Annahme: homogene edges

        #v_jac_array = zeros(for i in vertices Array{Array{Float64, 2}, 1}) # benutze list comprehension
        # initialisiere v_jac_array erstmal mit zeros, aber Dimension muss stimmen
        # benutze v_dims, falls du num_v brauchst, kannst du dir das einfach mit graph_structure.num_v etc. holen

        jac_graph_data = JacGraphData(v_jac_array, e_jac_array, e_jac_product, graph_structure) # Funktion

        t = 0.0 # any Float64
        # p später erstmal nothing
        nd_jac_vec_operator = NDJacVecOperator(similar(v_array), nothing, t, graph, graph_structure, graph_data, jac_graph_data, parallel) # x, p, t werden in update_coefficients geändert

        return ODEFunction(nd!; mass_matrix = mass_matrix, jac = nd_jac_vec_operator, syms = symbols)
    end

    return ODEFunction(nd!; mass_matrix = mass_matrix, syms = symbols)
end

end # module
