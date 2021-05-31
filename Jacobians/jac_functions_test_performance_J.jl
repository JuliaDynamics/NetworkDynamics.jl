using DifferentialEquations
using Reexport
using LinearAlgebra
using NetworkDynamics
using LightGraphs
using DiffEqBase
using BenchmarkTools
# using Profile
# using Traceur
# import DiffEqBase.update_coefficients!

N = 4
#k = 2
#g = barabasi_albert(N, k)

# building simple ring network consisting of 4 nodes
g = SimpleGraph(N)
add_edge!(g, 1, 2)
add_edge!(g, 2, 3)
add_edge!(g, 3, 4)

g
#add_edge!(g, 1, 3)


function diffusionedge!(e, v_s, v_d, p, t)
    e[1] = v_s[1] - v_d[1]
#    e[2] = v_s[2] - v_d[2]
    nothing
end

function diffusionvertex!(dv, v, edges, p, t)
    dv[1] = 0.
    dv[2] = 0.
    for e in edges
        dv[1] += e[1]
        dv[2] += e[1]
    end
    nothing
end

function jac_vertex!(J::AbstractMatrix, v, p, t)
    #J = internal_jacobian(J, v, p, t)
    J[1, 1] = 0.0
    J[2, 1] = 0.0
end

function jac_edge!(J_s::AbstractMatrix, J_d::AbstractMatrix, v_s, v_d, p, t)
   J_s[1, 1] = 1.0
#   J_s[2, 1] = 0.0

   J_d[1, 1] = -1.0
#   J_d[2, 1] = 0.0
end

nd_jac_vertex = ODEVertex(f! = diffusionvertex!, dim = 2, vertex_jacobian! = jac_vertex!)

nd_jac_edge = StaticEdge(f! = diffusionedge!, dim = 1, coupling = :undirected, edge_jacobian! = jac_edge!)

nd_jac = network_dynamics(nd_jac_vertex, nd_jac_edge, g, jac = true)



graph_structure = nd_jac.f.graph_structure
graph_data = nd_jac.f.graph_data

v_dims = nd_jac.f.graph_structure.v_dims
e_dims = nd_jac.f.graph_structure.e_dims
num_e = nd_jac.f.graph_structure.num_e

vertices! = nd_jac.f.vertices!
edges! = nd_jac.f.edges!

graph_structure.d_v

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
    e_jac_array = [[zeros(dim, srcdim), zeros(dim, dstdim)] for (dim, srcdim, dstdim) in zip(gs.v_dims, gs.v_dims, gs.v_dims)] # homogene Netzwerke: v_src_dim = v_dst_dim = v_dim
    e_jac_product = [zeros(gs.v_dims[1]) for i in 1:gs.num_e]
    JacGraphData1(v_jac_array, e_jac_array, e_jac_product)
end


## Idee: umschreiben in SArrays oder SVector
v_jac_array = [Array{Float64,2}(undef, dim, dim) for dim in v_dims]
e_jac_array = [[zeros(dim, srcdim), zeros(dim, dstdim)] for (dim, srcdim, dstdim) in zip(v_dims, v_dims, v_dims)] # homogene Netzwerke: v_src_dim = v_dst_dim = v_dim
e_jac_product =  [zeros(v_dims[1]) for i in 1:num_e] # e_dims

jac_graph_data_object = JacGraphData1(v_jac_array, e_jac_array, e_jac_product, graph_structure)
@btime jac_graph_data_object = JacGraphData1(v_jac_array, e_jac_array, e_jac_product, graph_structure)

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
        new{T, typeof(x), typeof(t), typeof(vertices!), typeof(edges!), typeof(graph), typeof(graph_data), typeof(jac_graph_data)}(x, p, t, vertices!, edges!, graph, graph_structure, graph_data, jac_graph_data, parallel)
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
@btime  NDJacVecOperator1(x, p, t, vertices!, edges!, g, graph_structure, graph_data, jac_graph_data_object, parallel)


### get functions for update_coefficients

@inline get_src_edge_jacobian(jgd::JacGraphData1, i::Int) = jgd.e_jac_array[i][1]
@inline get_dst_edge_jacobian(jgd::JacGraphData1, i::Int) = jgd.e_jac_array[i][2]
@inline get_vertex_jacobian(jgd::JacGraphData1, i::Int) = jgd.v_jac_array[i]


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

    @nd_threads Jac.parallel for i in 1:gs.num_v
        maybe_idx(Jac.vertices!, i).vertex_jacobian!(
          get_vertex_jacobian(jgd, i),
          get_vertex(gd, i),
          p_v_idx(p, i),
          t)
    end

    @nd_threads Jac.parallel for i in 1:gs.num_e
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

x_test2 = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]

dx_test = similar(x_test2)
#x = similar(zeros(1), sum(v_dims))
p_test = nothing
t_test = 0.0
parallel = false


function jac_vec_prod(Jac::NDJacVecOperator1, z)

    gs = Jac.graph_structure
    p = Jac.p
    # x = Jac.x
    checkbounds_p(p, gs.num_v, gs.num_e)
    jgd = Jac.jac_graph_data
#    zwischen_speicher = Jac.jac_graph_data.e_jac_product

    # first for loop that considers the mutliplication of each edge jacobians with the corresponding component of z
    for i in 1:gs.num_e
        jgd.e_jac_product[i] .= get_src_edge_jacobian(jgd, i) * view(z, gs.s_e_idx[i]) + get_dst_edge_jacobian(jgd, i) * view(z, gs.d_e_idx[i])

        #view(zwischen_speicher, i) = get_src_edge_jacobian(jgd, i) * view(z, gs.s_e_idx[i]) + get_dst_edge_jacobian(jgd, i) * view(z, gs.d_e_idx[i])

    end

    # in this function there is no dx in which the Jacobian can be stored, so an extra array must be created and returned
    dx = zeros(gs.v_dims[1], gs.num_v)

    # second for loop in which the multiplication of vertex jacobian and the corresponding component of z is done with addition of the e_jac_product to dx
    for i in 1:gs.num_v
        #view(dx, gs.v_idx[i]) .= get_vertex_jacobian(jgd, i) * view(z, gs.v_idx[i]) + sum!([0.0], view(jgd.e_jac_product, gs.d_e_idx[i])[1])
        view(dx, gs.v_idx[i]) .= get_vertex_jacobian(jgd, i) * view(z, gs.v_idx[i])
        view_e_jac_product = view(jgd.e_jac_product, gs.d_v[i])
        if view_e_jac_product == []
            view(dx, gs.v_idx[i]) .+= zeros(gs.v_dims[1])
        else
            view(dx, gs.v_idx[i]) .+= sum(view_e_jac_product)
        # view(dx, gs.v_idx[i]) .+= sum(view(jgd.e_jac_product, gs.d_v[i]))
        end
    end
    return dx
end

update_coefficients1!(NDJacVecOp, x_test2, p_test, t_test)

jac_vec_prod(NDJacVecOp, x_test2)

test_z_src = view(x_test2, NDJacVecOp.graph_structure.s_e_idx[1])
test_z_dst = view(x_test2, NDJacVecOp.graph_structure.d_e_idx[1])


for i in 1:NDJacVecOp.graph_structure.num_e
    NDJacVecOp.jac_graph_data.e_jac_product[i] .= get_src_edge_jacobian(NDJacVecOp.jac_graph_data, i) * view(x_test2, NDJacVecOp.graph_structure.s_e_idx[i]) + get_dst_edge_jacobian(NDJacVecOp.jac_graph_data, i) * view(x_test2, NDJacVecOp.graph_structure.d_e_idx[i])
end

print(NDJacVecOp.jac_graph_data.e_jac_product)

function jac_vec_prod!(dx, Jac::NDJacVecOperator1, z)

    gs = Jac.graph_structure
    p = Jac.p
    x = Jac.x
    checkbounds_p(p, gs.num_v, gs.num_e)
    jgd = Jac.jac_graph_data
    #zwischen_speicher = Jac.jac_graph_data.e_jac_product

    for i in 1:gs.num_e
        jgd.e_jac_product[i] .= get_src_edge_jacobian(jgd, i) * view(z, gs.s_e_idx[i]) + get_dst_edge_jacobian(jgd, i) * view(z, gs.d_e_idx[i])
    end

    for i in 1:gs.num_v
        # view(dx, gs.v_idx[i]) .= get_vertex_jacobian(jgd, i) * view(z, gs.v_idx[i])
        # view(dx, gs.v_idx[i]) .+= sum!([0.0], view(jgd.e_jac_product, gs.d_e_idx[i])[1])
        #view(dx, gs.v_idx[i]) .= get_vertex_jacobian(jgd, i) * view(z, gs.v_idx[i]) + sum!([0.0], view(jgd.e_jac_product, gs.d_e_idx[i])[1])
        view(dx, gs.v_idx[i]) .= get_vertex_jacobian(jgd, i) * view(z, gs.v_idx[i]) + sum(view(jgd.e_jac_product, gs.d_e_idx[i]))
    end
end

@btime jac_vec_prod!(dx_test, NDJacVecOp, x_test)

# functions for callable structs

Base.:*(Jac::NDJacVecOperator1, z::AbstractVector) = jac_vec_prod(Jac, z)

function LinearAlgebra.mul!(dx::AbstractVector, Jac::NDJacVecOperator1, z::AbstractVector)
    jac_vec_prod!(dx, Jac, z)
end

# callable structs

function (Jac::NDJacVecOperator1)(x, p, t) # auch Number bei t?
    update_coefficients1!(Jac, x, p, t)
    Jac*x
end

function (Jac::NDJacVecOperator1)(dx, x, p, t)
    update_coefficients1!(Jac, x, p, t)
    mul!(dx, Jac, x)
end



## check computation time and allocations

x_test = [1.0, 2.0, 3.0, 4.0]
dx_test = similar(x_test)
#x = similar(zeros(1), sum(v_dims))
p_test = nothing
t_test = 0.0
parallel = false


update_coefficients1!(NDJacVecOp, x_test, p_test, t_test)

NDJacVecOp.jac_graph_data.e_jac_array

@btime update_coefficients1!(NDJacVecOp, x_test, p_test, t_test) # first run: 0 allocations

@btime NDJacVecOp(x_test, p_test, t_test) # first run: 37 allocations - new: 33 alocations



# @trace(NDJacVecOp(x_test, p_test, t_test), modules=[Main])


# @profile NDJacVecOp(x_test, p_test, t_test) # first run: 37 allocations


@btime NDJacVecOp(dx_test, x_test, p_test, t_test) # first run: 20 allocations
#@allocated call_callable_struct_2 = NDJacVecOp(dx_test, x_test, p_test, t_test)

##

e_jac_product1 = NDJacVecOp.jac_graph_data.e_jac_product
e_jac_product2 = NDJacVecOp.jac_graph_data.e_jac_product

for i in 1:NDJacVecOp.graph_structure.num_e
    e_jac_product1[i] .= get_src_edge_jacobian(NDJacVecOp.jac_graph_data, i) * view(x_test, NDJacVecOp.graph_structure.s_e_idx[i]) + get_dst_edge_jacobian(NDJacVecOp.jac_graph_data, i) * view(x_test, NDJacVecOp.graph_structure.d_e_idx[i])
end

e_jac_product1

for i in 1:NDJacVecOp.graph_structure.num_e
    zwischen_speicher2 = e_jac_product2
    #display(@allocated jgd.e_jac_product[i] .= get_src_edge_jacobian(jgd, i) * view(z, gs.s_e_idx[i]) + get_dst_edge_jacobian(jgd, i) * view(z, gs.d_e_idx[i]))
    view(zwischen_speicher2, i) = get_src_edge_jacobian(NDJacVecOp.jac_graph_data, i) * view(x_test, NDJacVecOp.graph_structure.s_e_idx[i]) + get_dst_edge_jacobian(NDJacVecOp.jac_graph_data, i) * view(x_test, NDJacVecOp.graph_structure.d_e_idx[i])
    println(zwischen_speicher2[i])
    #view(NDJacVecOp.jac_graph_data.e_jac_product, i) = get_src_edge_jacobian(NDJacVecOp.jac_graph_data, i) * view(x_test, NDJacVecOp.graph_structure.s_e_idx[i]) + get_dst_edge_jacobian(NDJacVecOp.jac_graph_data, i) * view(x_test, NDJacVecOp.graph_structure.d_e_idx[i])

end




##

@time call_update_coefficients! = update_coefficients1!(NDJacVecOp, x_test, p_test, t_test)

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
