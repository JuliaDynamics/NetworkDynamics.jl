using DifferentialEquations
using Reexport
using LinearAlgebra
using NetworkDynamics
using LightGraphs
using DiffEqBase
using BenchmarkTools
# import DiffEqBase.update_coefficients!

N = 4
#k = 2
#g = barabasi_albert(N, k)
# building simple ring network consisting of 4 nodes
g = SimpleGraph(N)
add_edge!(g, 1, 2)
add_edge!(g, 2, 3)
add_edge!(g, 3, 4)
add_edge!(g, 4, 1)

function diffusionedge!(e, v_s, v_d, p, t)
    e[1] = v_s[1] - v_d[1]
    e[2] = v_s[2] - v_d[2]
    nothing
end

function diffusionvertex!(dv, v, edges, p, t)
    dv[1] = 0.
    dv[2] = 0.
    for e in edges
        dv[1] += e[1]
        dv[2] += e[2]
    end
    nothing
end

function jac_vertex!(J::AbstractMatrix, v, p, t)
    #J = internal_jacobian(J, v, p, t)
    J[1, 1] = 1.0
    J[2, 2] = 1.0
end

function jac_edge!(J_s::AbstractMatrix, J_d::AbstractMatrix, v_s, v_d, p, t)
   J_s[1, 1] = 1.0
   J_s[2, 1] = 1.0

   J_d[1, 1] = -1.0
   J_d[2, 1] = -1.0
end

nd_jac_vertex = ODEVertex(f! = diffusionvertex!, dim = 2, vertex_jacobian! = jac_vertex!)

nd_jac_edge = StaticEdge(f! = diffusionedge!, dim = 2, edge_jacobian! = jac_edge!)

nd_jac = network_dynamics(nd_jac_vertex, nd_jac_edge, g, jac = true)
nd = network_dynamics(nd_jac_vertex, nd_jac_edge, g, jac = false)

@btime network_dynamics(nd_jac_vertex, nd_jac_edge, g, jac = true) # 705 allocations
@btime network_dynamics(nd_jac_vertex, nd_jac_edge, g, jac = false) # 598 allocations

@btime graph_structure = nd_jac.f.graph_structure # 3 allocations
@btime graph_data = nd_jac.f.graph_data #3

@btime v_dims = nd_jac.f.graph_structure.v_dims #4
@btime e_dims = nd_jac.f.graph_structure.e_dims #4
@btime num_e = nd_jac.f.graph_structure.num_e #4
@btime vertices! = nd_jac.f.vertices! #3
@btime edges! = nd_jac.f.edges! #3

# Vgl NetworkDynamics
@btime graph_structure = nd.f.graph_structure #3
@btime graph_data = nd.f.graph_data #3
@btime v_dims = nd.f.graph_structure.v_dims #4
@btime e_dims = nd.f.graph_structure.e_dims #4
@btime num_e = nd.f.graph_structure.num_e #4
@btime vertices! = nd.f.vertices! #3
@btime edges! = nd.f.edges!#3

#3
graph_structure = nd_jac.f.graph_structure # 3 allocations
graph_data = nd_jac.f.graph_data #3

v_dims = nd_jac.f.graph_structure.v_dims #4
e_dims = nd_jac.f.graph_structure.e_dims #4
num_e = nd_jac.f.graph_structure.num_e #4
vertices! = nd_jac.f.vertices! #3
edges! = nd_jac.f.edges!

# for vertices! and edges! we need to get the index
@inline Base.@propagate_inbounds function maybe_idx(p::T, i) where T <: AbstractArray
    p[i]
end

@inline function maybe_idx(p, i)
    p
end

struct JacGraphData1
    v_jac_array::Array{Array{Float64, 2}, 1}
    e_jac_array::Array{Array{Array{Float64, 2}, 1}, 1}
    e_jac_product::Array{Array{Float64, 1}, 1}
end

function JacGraphData1(v_jac_array, e_jac_array, e_jac_product_array, gs::GraphStruct)
    v_jac_array = [Array{Float64,2}(undef, dim, dim) for dim in gs.v_dims]
    e_jac_array = [[zeros(dim, srcdim), zeros(dim, dstdim)] for (dim, srcdim, dstdim) in zip(gs.e_dims, gs.v_dims, gs.v_dims)] # homogene Netzwerke: v_src_dim = v_dst_dim = v_dim
    e_jac_product = [zeros(gs.e_dims[1]) for i in 1:gs.num_e]
    JacGraphData1(v_jac_array, e_jac_array, e_jac_product)
end

## Idee: umschreiben in SArrays oder SVector
v_jac_array = [Array{Float64,2}(undef, dim, dim) for dim in v_dims]
e_jac_array = [[zeros(dim, srcdim), zeros(dim, dstdim)] for (dim, srcdim, dstdim) in zip(e_dims, v_dims, v_dims)] # homogene Netzwerke: v_src_dim = v_dst_dim = v_dim
e_jac_product =  [zeros(e_dims[1]) for i in 1:num_e]

@btime jac_graph_data_object = JacGraphData1(v_jac_array, e_jac_array, e_jac_product, graph_structure) #24
jac_graph_data_object = JacGraphData1(v_jac_array, e_jac_array, e_jac_product, graph_structure) #24
#### NDJacVecOperator1

mutable struct NDJacVecOperator1{T, uType, tType, T1, T2, G, GD, JGD} <: DiffEqBase.AbstractDiffEqLinearOperator{T} # mutable da x, p, t geupdated werden
    x::uType
    p
    t::tType
    vertices!::T1
    edges!::T2
    graph::G
    graph_structure::GraphStruct
    graph_data::GD
    jac_graph_data::JGD
    parallel::Bool

    function NDJacVecOperator1{T}(x, p, t, vertices!, edges!, graph, graph_structure, graph_data, jac_graph_data, parallel) where T
        new{T, typeof(x), typeof(t), typeof(vertices!), typeof(edges!),typeof(graph), typeof(graph_data), typeof(jac_graph_data)}(x, p, t, vertices!, edges!, graph, graph_structure, graph_data, jac_graph_data, parallel)
    end

    function NDJacVecOperator1(x, p, t, vertices!, edges!, graph, graph_structure, graph_data, jac_graph_data, parallel)
        NDJacVecOperator1{eltype(x)}(x, p, t, vertices!, edges!, graph, graph_structure, graph_data, jac_graph_data, parallel)
    end
end

x = similar(zeros(1), sum(v_dims))
p = nothing
t = 0.0
parallel = false

NDJacVecOp = NDJacVecOperator1(x, p, t, vertices!, edges!, g, graph_structure, graph_data, jac_graph_data_object, parallel) # 51 allocations
@btime  NDJacVecOperator1(x, p, t, vertices!, edges!, g, graph_structure, graph_data, jac_graph_data_object, parallel) #51
### get functions for update_coefficients

@inline get_src_edge_jacobian(jgd::JacGraphData1, i::Int) = jgd.e_jac_array[i][1]
@inline get_dst_edge_jacobian(jgd::JacGraphData1, i::Int) = jgd.e_jac_array[i][2]
@inline get_vertex_jacobian(jgd::JacGraphData1, i::Int) = jgd.v_jac_array[i]

@btime @inline get_vertex_jacobian(jgd::JacGraphData1, i::Int) = jgd.v_jac_array[i] #0
# functions needed for update_coefficients1!

@inline function prep_gd(dx::AbstractArray{T}, x::AbstractArray{T}, gd::GraphData{GDB, T, T}, gs) where {GDB, T}
    if size(x) == (gs.dim_v,)
         swap_v_array!(gd, x)
         return gd
    else
         error("Size of x does not match the dimension of the system.")
    end
end

@inline function prep_gd(dx, x, gd, gs)
    # Type mismatch
    if size(x) == (gs.dim_v,)
        e_array = similar(dx, gs.dim_e)
        return GraphData(x, e_array, gs)
    else
        error("Size of x does not match the dimension of the system.")
   end
end


function update_coefficients1!(Jac::NDJacVecOperator1, x, p, t)

    gs = Jac.graph_structure
    checkbounds_p(p, gs.num_v, gs.num_e)
    gd = prep_gd(x, x, Jac.graph_data, Jac.graph_structure)
    jgd = Jac.jac_graph_data

    for i in 1:gs.num_v
        maybe_idx(Jac.vertices!, i).vertex_jacobian!(
          get_vertex_jacobian(jgd, i),
          get_vertex(gd, i),
          p_v_idx(p, i),
          t)
    end

    for i in 1:gs.num_e
          maybe_idx(Jac.edges!, i).edge_jacobian!(
              get_src_edge_jacobian(jgd, i),
              get_dst_edge_jacobian(jgd, i),
              get_src_vertex(gd, i),
              get_dst_vertex(gd, i),
              p_e_idx(p, i),
              t)
      end

    Jac.x = x
    Jac.p = p
    Jac.t = t
end



function jac_vec_prod(Jac::NDJacVecOperator1, z)

    gs = Jac.graph_structure
    p = Jac.p
    #x = Jac.x
    checkbounds_p(p, gs.num_v, gs.num_e)
    jgd = Jac.jac_graph_data

    # first for loop that considers the mutliplication of each edge jacobians with the corresponding component of z
    for i in 1:gs.num_e
        jgd.e_jac_product[i] .= get_src_edge_jacobian(jgd, i) * view(z, gs.s_e_idx[i]) + get_dst_edge_jacobian(jgd, i) * view(z, gs.d_e_idx[i])
        #view(jgd.e_jac_product,i, i) .= get_src_edge_jacobian(jgd, i) * view(z, gs.s_e_idx[i]) + get_dst_edge_jacobian(jgd, i) * view(z, gs.d_e_idx[i])
    end

    # in this function there is no dx in which the Jacobian can be stored, so an extra array must be created and returned
    dx = zeros(gs.v_dims[1], gs.num_v)

    # second for loop in which the multiplication of vertex jacobian and the corresponding component of z is done with addition of the e_jac_product to dx
    for i in 1:gs.num_v
        #view(dx, gs.v_idx[i]) .= get_vertex_jacobian(jgd, i) * view(z, gs.v_idx[i])
        #view(dx, gs.v_idx[i]) .+= sum!([0.0], view(jgd.e_jac_product, gs.d_e_idx[i])[1])
        view(dx, gs.v_idx[i]) .= get_vertex_jacobian(jgd, i) * view(z, gs.v_idx[i]) + sum!([0.0], view(jgd.e_jac_product, gs.d_e_idx[i])[1])
        # view(dx, gs.v_idx[i]) .+= sum(sum(view(jgd.e_jac_product, gs.d_e_idx[i])))
    end
    return dx
end

gs = NDJacVecOp.graph_structure
gs.s_e_idx[2]

@btime jac_vec_prod(NDJacVecOp,x_test) #37
@btime jac_vec_prod(NDJacVecOp,x_test) # 33
function e_jac_prod!(z, gs, jgd)
    for i in 1:gs.num_e
    #for i in enumerate(gs.num_v)
        #println(enumerate(gs.num_v))
        jgd.e_jac_product[i] .= get_src_edge_jacobian(jgd, i) * view(z, gs.s_e_idx[i]) + get_dst_edge_jacobian(jgd, i) * view(z, gs.d_e_idx[i])
    end
end;
 gs.num_e
jgd.e_jac_product
jgd.e_jac_product[1]
view(jgd.e_jac_product, 1 ,1)
gs.e_idx[1]
gs.e_idx[2]
gs.d_e_idx[1]
gs.v_idx[1]
function v_jac_prod!(dx, z, gs, jgd)
    for i in 1:gs.num_v
        #view(dx, gs.v_idx[i]) .= get_vertex_jacobian(jgd, i) * view(z, gs.v_idx[i])
        #view(dx, gs.v_idx[i]) .+= sum!([0.0], view(jgd.e_jac_product, gs.d_e_idx[i])[1])
        view(dx, gs.v_idx[i]) .= get_vertex_jacobian(jgd, i) * view(z, gs.v_idx[i]) + sum!([0.0], view(jgd.e_jac_product, gs.d_e_idx[i])[1])
        # view(dx, gs.v_idx[i]) .+= sum(sum(view(jgd.e_jac_product, gs.d_e_idx[i])))
    end
end;

gs = NDJacVecOp.graph_structure
p = NDJacVecOp.p
x = NDJacVecOp.x
@btime checkbounds_p(p, gs.num_v, gs.num_e)
jgd = NDJacVecOp.jac_graph_data # 2 allocations
dx = zeros(gs.v_dims[1], gs.num_v)
gs.num_v
z = x_test

@btime e_jac_prod!(x_test,gs,jgd) # 12 allocations
i = 1
@btime jgd.e_jac_product[i] #1
@btime get_src_edge_jacobian(jgd, i) # 0 allocations
@btime view(z, gs.s_e_idx[i]) # 3 allocations
@btime get_dst_edge_jacobian(jgd, i) # 0 allocations
@btime view(z, gs.d_e_idx[i]) # 3
@btime (get_dst_edge_jacobian(jgd, i) * view(z, gs.d_e_idx[i])) # 4 allocations
@btime (jgd.e_jac_product[i] .= get_src_edge_jacobian(jgd, i) * view(z, gs.s_e_idx[i]) + get_dst_edge_jacobian(jgd, i) * view(z, gs.d_e_idx[i])) # 12 allocations
@btime (get_src_edge_jacobian(jgd, i) * view(z, gs.s_e_idx[i]) + get_dst_edge_jacobian(jgd, i) * view(z, gs.d_e_idx[i]))
@btime (get_src_edge_jacobian(jgd, i) * z[i]) # 2
z[gs.s_e_idx[i]]
@btime (get_src_edge_jacobian(jgd, i) * z[gs.s_e_idx[i]]) # 4
@btime gs.s_e_idx[i]
@btime e_index = gs.s_e_idx
@btime e_index[i]


@btime v_jac_prod!(dx, x_test, gs, jgd) # 24 allocations
v_jac_prod!(dx, x_test, gs, jgd)
function jac_vec_prod(Jac::NDJacVecOperator1, z)

    gs = Jac.graph_structure
    p = Jac.p
    x = Jac.x
    checkbounds_p(p, gs.num_v, gs.num_e)
    jgd = Jac.jac_graph_data

    #ejp = zeros(gs.num_e, gs.e_dims[1])
    e_jac_prod!(z, gs, jgd)

    dx = zeros(gs.v_dims[1], gs.num_v)
    v_jac_prod!(dx, z, gs, jgd)

    return dx

end
@btime jac_vec_prod(NDJacVecOp,x) #37

function jac_vec_prod!(dx, Jac::NDJacVecOperator1, z)

    gs = Jac.graph_structure
    p = Jac.p
    x = Jac.x
    checkbounds_p(p, gs.num_v, gs.num_e)
    jgd = Jac.jac_graph_data

    for i in 1:gs.num_e
        display(@allocated jgd.e_jac_product[i] .= get_src_edge_jacobian(jgd, i) * view(z, gs.s_e_idx[i]) + get_dst_edge_jacobian(jgd, i) * view(z, gs.d_e_idx[i]))
    end

    for i in 1:gs.num_v
        view(dx, gs.v_idx[i]) .= get_vertex_jacobian(jgd, i) * view(z, gs.v_idx[i])
        view(dx, gs.v_idx[i]) .+= sum!([0.0], view(jgd.e_jac_product, gs.d_e_idx[i])[1])
        # view(dx, gs.v_idx[i]) .+= sum(sum(view(jgd.e_jac_product, gs.d_e_idx[i])))
    end
end
gs.s_e_idx[i]

# functions for callable structs

Base.:*(Jac::NDJacVecOperator1, z::AbstractVector) = jac_vec_prod(Jac, z)

function LinearAlgebra.mul!(dx::AbstractVector, Jac::NDJacVecOperator1, z::AbstractVector)
    jac_vec_prod!(dx, Jac, z)
end

# callable structs

function (Jac::NDJacVecOperator1)(x, p, t) # auch Number bei t?
    update_coefficients!(Jac, x, p, t)
    Jac*x
end

function (Jac::NDJacVecOperator1)(dx, x, p, t)
    update_coefficients!(Jac, x, p, t)
    mul!(dx, Jac, x)
end


## check computation time and allocations

x_test = [1.0, 2.0, 3.0, 4.0]
dx_test = similar(x_test)
#x = similar(zeros(1), sum(v_dims))
p_test = nothing
t_test = 0.0
parallel = false

@btime update_coefficients!(NDJacVecOp, x_test, p_test, t_test)
@allocated call_update_coefficients! = update_coefficients!(NDJacVecOp, x_test, p_test, t_test)

@btime NDJacVecOp(x_test, p_test, t_test)
@allocated call_callable_struct_1 = NDJacVecOp(x_test, p_test, t_test)

@btime NDJacVecOp(dx_test, x_test, p_test, t_test)
@allocated call_callable_struct_2 = NDJacVecOp(dx_test, x_test, p_test, t_test)




@time call_update_coefficients! = update_coefficients!(NDJacVecOp, x_test, p_test, t_test)
@btime update_coefficients!(NDJacVecOp, x_test, p_test, t_test)
@time call_callable_struct_1 = NDJacVecOp(x_test, p_test, t_test)

@time call_callable_struct_2 = NDJacVecOp(dx_test, x_test, p_test, t_test)


@allocated NDJacVecOp(x_test, p_test, t_test)
@allocated NDJacVecOp(dx_test, x_test, p_test, t_test)


# outsourced for loops + corresponding functions

function e_jac_prod!(z, ejp, gs, jgd)
    for i in 1:gs.num_e
        ejp[i, :] .= get_src_edge_jacobian(jgd, i) * view(z, gs.s_e_idx[i]) + get_dst_edge_jacobian(jgd, i) * view(z, gs.d_e_idx[i])
    end
end;

function v_jac_prod!(dx, z, ejp, gs, jgd)
    for i in 1:gs.num_v
        view(dx, gs.v_idx[i]) .= get_vertex_jacobian(jgd, i) * view(z, gs.v_idx[i])
        view(dx, gs.v_idx[i]) .+= vec(sum(view(ejp, gs.d_v[i], :), dims = 1))
    end
end;


function jac_vec_prod(Jac::NDJacVecOperator1, z)

    gs = Jac.graph_structure
    checkbounds_p(p, gs.num_v, gs.num_e)
    gd = prep_gd(x, x, Jac.graph_data, Jac.graph_structure)
    jgd = Jac.jac_graph_data

    ejp = zeros(gs.num_e, gs.e_dims[1])
    e_jac_prod!(z, ejp, gs, jgd)

    dx = zeros(gs.v_dims[1]*gs.num_v)
    v_jac_prod!(dx, z, ejp, gs, jgd)

    return dx

end

function jac_vec_prod!(dx, Jac::NDJacVecOperator1, z)

    gs = Jac.graph_structure
    checkbounds_p(p, gs.num_v, gs.num_e)
    gd = prep_gd(x, x, Jac.graph_data, Jac.graph_structure)
    jgd = Jac.jac_graph_data

    e_jac_prod!(z, e_jac_product, gs, jgd)
    v_jac_prod!(dx, z, e_jac_product, gs, jgd)

end
