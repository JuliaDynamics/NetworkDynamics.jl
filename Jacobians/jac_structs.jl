
################################################################################

"""
Together with Constructors this module forms the backbone of the core API.
It provide the basic types to construct Arrays of VertexFunction and
EdgeFunction which can be handled by network_dynamics.
"""
module ComponentFunctions

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
                           edge_jacobian!) where T # beliebiger Typ

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

######################## Delay Block ###########################################
"""
    DDEVertex(f!, dim, mass_matrix, sym)

Wrapper that ensures compatibility of a **mutating** function **`f!`** with
the key constructor `network_dynamics`.

**`f!`**  describes the local behaviour at a dynamic node and has to respect
the following calling syntax

```julia
f!(dv, v, edges, h, p, t) -> nothing
```

Here `dv`, `v`, `p` and `t` are the usual ODE arguments, while
`edges` is an Arry of incoming edges. `h` is the history array for `v`.

- `dim` is the number of independent variables in the edge equations and
- `sym` is an array of symbols for these variables.
- `coupling` is a Symbol describing if the EdgeFunction is intended for a directed graph (`:directed`) or for an undirected graph (`{:undirected, :fiducial}`). `:directed` is intended for directed graphs. `:undirected` is the default option and is only compatible with SimpleGraph. in this case f! should specify the coupling from a source vertex to a destination vertex.  `:fiducial` lets the user specify both the coupling from src to dst, as well as the coupling from dst to src and is intended for advanced users.
- `mass_matrix` is an optional argument that defaults to the identity
matrix `I`. If a mass matrix M is given the system `M * de = f!` will be
solved.

For more details see the documentation.
"""
@Base.kwdef struct DDEVertex{T} <: VertexFunction
    f!::T # The function with signature (dx, x, edges, h, p, t) -> nothing
    dim::Int # number of dimensions of x
    mass_matrix=I # Mass matrix for the equation
    sym=[:v for i in 1:dim] # Symbols for the dimensions
end

function DDEVertex(ov::ODEVertex)
    let _f! = ov.f!, dim = ov.dim, sym = ov.sym, mass_matrix = ov.mass_matrix
        f! = @inline (dv, v, dst_edges, h_v, p, t) -> _f!(dv, v, dst_edges, p, t)
        return DDEVertex(f!, dim, mass_matrix, sym)
    end
end


function DDEVertex(sv::StaticVertex)
    let _f! = ODE_from_Static(sv.f!), dim = sv.dim, sym = sv.sym
        f! = @inline (dv, v, dst_edges, h_v, p, t) -> _f!(dv, v, dst_edges, p, t)
        return DDEVertex(f!, dim, 0., sym)
    end
end

# Promotion rules [eventually there might be too many types to hardcode everyhting]

convert(::Type{DDEVertex}, x::StaticVertex) = DDEVertex(x)
promote_rule(::Type{DDEVertex}, ::Type{StaticVertex}) = DDEVertex


convert(::Type{DDEVertex}, x::ODEVertex) = DDEVertex(x)
promote_rule(::Type{DDEVertex}, ::Type{ODEVertex}) = DDEVertex


"""
Like a static edge but with extra arguments for the history of the source and destination vertices. This is NOT a DDEEdge.
"""
@Base.kwdef struct StaticDelayEdge{T} <: EdgeFunction
    f!::T # (e, v_s, v_t, h_v_s, h_v_d, p, t) -> nothing
    dim::Int # number of dimensions of e
    coupling = :undefined # :directed, :symmetric, :antisymmetric, :fiducial, :undirected
    sym=[:e for i in 1:dim] # Symbols for the dimensions


    function StaticDelayEdge(user_f!::T,
                           dim::Int,
                           coupling::Symbol,
                           sym::Vector{Symbol}) where T

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
            f! = @inline (e, v_s, v_d, h_v_s, h_v_d, p, t) -> begin
                @inbounds user_f!(view(e,1:dim), v_s, v_d, h_v_s, h_v_d, p, t)
                @inbounds user_f!(view(e,dim+1:2dim), v_d, v_s, h_v_d, h_v_s, p, t)
                nothing
            end
        elseif coupling == :antisymmetric
            f! = @inline (e, v_s, v_d, h_v_s, h_v_d, p, t) -> begin
                @inbounds user_f!(view(e,1:dim), v_s, v_d, h_v_s, h_v_d, p, t)
                @inbounds for i in 1:dim
                    e[dim + i] = -1.0 * e[i]
                end
                nothing
            end
        elseif coupling == :symmetric
            f! = @inline (e, v_s, v_d, h_v_s, h_v_d, p, t) -> begin
                @inbounds user_f!(view(e,1:dim), v_s, v_d, h_v_s, h_v_d, p, t)
                @inbounds for i in 1:dim
                    e[dim + i] = e[i]
                end
                nothing
            end
        end
        return new{typeof(f!)}(f!, 2dim, coupling, repeat(sym, 2))
    end
end

function StaticDelayEdge(se::StaticEdge)
    let _f! = se.f!, dim = se.dim, coupling = se.coupling, sym = se.sym
        f! = @inline (e, v_s, v_d, h_v_s, h_v_d, p, t) -> _f!(e, v_s, v_d, p, t)
        if coupling ∈ (:undefined, :fiducial, :directed)
            return StaticDelayEdge(f!, dim, coupling, sym)
        else
            return StaticDelayEdge(f!, dim, :fiducial, sym)
        end
    end
end

# Promotion rules

convert(::Type{StaticDelayEdge}, x::StaticEdge) = StaticDelayEdge(x)
promote_rule(::Type{StaticDelayEdge}, ::Type{StaticEdge}) = StaticDelayEdge



convert(::Type{ODEVertex}, x::StaticVertex) = ODEVertex(x)
promote_rule(::Type{ODEVertex}, ::Type{StaticVertex}) = ODEVertex

# Not sure if the next line does something?
promote_rule(::Type{ODEVertex{T}}, ::Type{ODEVertex{U}}) where {T, U} = ODEVertex

convert(::Type{ODEEdge}, x::StaticEdge) = ODEEdge(x)
promote_rule(::Type{ODEEdge}, ::Type{StaticEdge}) = ODEEdge
"""
Promotes a StaticVertex to an ODEVertex.
"""
struct ODE_from_Static{T}
    f!::T
end
function (ofs::ODE_from_Static)(dx,x,args...)
    # If mass matrix = 0 the differential equation sets dx = 0.
    # To set x to the value calculated by f! we first write the value calculated
    # by f! into dx, then subtract x. This leads to the  constraint
    # 0 = - x + f(...)
    # where f(...) denotes the value that f!(a, ...) writes into a.
    ofs.f!(dx, args...)
    @inbounds for i in eachindex(dx)
        dx[i] -= x[i]
    end
    nothing
end

function ODEVertex(sv::StaticVertex)
    let _f! = sv.f!, dim = sv.dim, sym = sv.sym
        return ODEVertex(ODE_from_Static(_f!), dim, 0., sym)
    end
end

function ODEEdge(se::StaticEdge)
    let _f! = se.f!, dim = se.dim, coupling = se.coupling, sym = se.sym
        if coupling == :undirected
            # undirected means the reconstruction has already happend and the edge may now
            # be considered fiducial
            return ODEEdge(ODE_from_Static(_f!), dim, :fiducial, 0., sym)
        else
            return ODEEdge(ODE_from_Static(_f!), dim, coupling, 0., sym)
        end
    end
end

end

######################## Delay Block Ende ######################################

######################## Graph Data Block ######################################
### TO DO: in NetworkStructures.jl
# field names von graph data


# Array{T,N}(undef, dims), N ist Dimension des arrays, T ist Typ des arrays
struct EdgeData{GDB, elE} <: AbstractArray{elE, 1}
    gdb::GDB
    idx_offset::Int
    len::Int
end

struct VertexData{GDB, elV} <: AbstractArray{elV, 1}
    gdb::GDB
    idx_offset::Int
    len::Int
end

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


### auslagern
    v_jac_array::Array{Array{Float, 2}, 1} # Dimension richtig? da Jacobians (maximal) 2 dim. # oder elVJac
    e_jac_array::Array{Array{Array{Float, 2}, 1}, 1}
    e_jac_product::Array{Float, 2}
end

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

### auslagern
    v_jac_array = [Array{Float64,2}(undef, dim, dim) for dim in gs.v_dims]
    # For edges there is another array layer in the buffer since each edge has two Jacobians
    #e_jac_array = [[zeros(dim, srcdim), zeros(dim, dstdim)] for (dim, srcdim, dstdim) in zip(gs.e_dims, v_src_dims, v_dst_dims)]
    e_jac_array = [[zeros(dim, srcdim), zeros(dim, dstdim)] for (dim, srcdim, dstdim) in zip(gs.e_dims, gs.v_dims, gs.v_dims)] # homogene Netzwerke: v_src_dim = v_dst_dim = v_dim
    e_jac_product =  zeros(gs.e_dims[1], gs.num_e) # Annahme: homogene edges
end

### TO DO: in NetworkStructures.jl

#@inline function get_src_edge_jacobian(gd::GraphData, i::Int) = gd.e_jac_array[i, 1] ## Micha: gd.e_jac_array[i][1]
@inline function get_src_edge_jacobian(gd::JacGraphData, i::Int) = gd.e_jac_array[i][1]

@inline function get_dst_edge_jacobian(gd::JacGraphData, i::Int) = gd.e_jac_array[i][2]

@inline function get_vertex_jacobian(gd::JacGraphData, i::Int) = gd.v_jac_array[i]


############################ NDJacVecOperator ##################################
mutable struct NDJacVecOperator{uType, P, tType, G, GD} <: DiffEqBase.AbstractDiffEqLinearOperator{T} # mutable da x, p, t geupdated werden
    x::uType # Punkt u, für welchen das jvp berechnet wird
    p::P
    t::tType
    graph::G
    graph_structure::GraphStruct
    graph_data::GD
    parallel::Bool
end

## solle wir wie chris eine mutable struct machen?
mutable struct NDJacVecOperator{uType, P, tType,G,GraphStruct,GD,Bool}# <: DiffEqBase.AbstractDiffEqLinearOperator{T}
    x::uType # Punkt u, für welchen das jvp berechnet wird
    p::P
    t::tType
    graph::G
    graph_structure::GraphStruct
    graph_data::GD
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

    if jac == true
        # These additional arrays are used for initializing the GraphData and will be overwritten
        # v_jac_array = similar(x_prototype, sum(v_dims))
        # e_jac_array = similar(x_prototype, sum(e_dims))
        # e_jac_product = similar(x_prototype, sum(num_e))

        graph_stucture = GraphStruct(graph, v_dims, e_dims, symbols_v, symbols_e)

        graph_data = GraphData(v_array, e_array, graph_stucture, v_jac_array, e_jac_array, e_jac_product)

        nd! = nd_ODE_Static(vertices!, edges!, graph, graph_stucture, graph_data, parallel)
        mass_matrix = construct_mass_matrix(mmv_array, graph_stucture)
        nd_jac_operator = NDJacVecOperator

        ODEFunction(nd!; mass_matrix = mass_matrix, jac = nd_jac_operator, syms = symbols)
    end

    graph_stucture = GraphStruct(graph, v_dims, e_dims, symbols_v, symbols_e)

    graph_data = GraphData(v_array, e_array, graph_stucture)


    nd! = nd_ODE_Static(vertices!, edges!, graph, graph_stucture, graph_data, parallel)
    mass_matrix = construct_mass_matrix(mmv_array, graph_stucture)

    ODEFunction(nd!; mass_matrix = mass_matrix, syms = symbols)
end
