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
    e[1] = v_s[1] - v_d[2]
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

#function jac_vertex!(J::AbstractMatrix, v, p, t)
#    J[1, 1] = 0.0
#    J[2, 1] = 0.0
#end

@Base.propagate_inbounds function jac_vertex!(J_v::AbstractMatrix, J_e::AbstractMatrix, v, p, t)
    # J_v = internal_jacobian(J, v, p, t)
    J_v[1, 1] = 0.0
    J_v[2, 1] = 0.0

    J_e[1, 1] = 1.
    J_e[2, 1] = 1.
    nothing
end

function jac_edge!(J_s::AbstractMatrix, J_d::AbstractMatrix, v_s, v_d, p, t)
   J_s[1, 1] = 1.0
#   J_s[1, 2] = 0.0

   J_d[1, 1] = -1.0
#   J_d[1, 2] = 0.0
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
    v_jac_array::Array{Array{Array{Float64, 2}, 1}, 1} # contains the jacobians for each vertex
    e_jac_array::Array{Array{Array{Float64, 2}, 1}, 1} # contains the jacobians for each edge
    e_jac_product::Array{Array{Float64, 1}, 1} # is needed later in jac_vec_prod(!) as a storage for the products of edge jacobians and vectors z
end

function JacGraphData1(gs::GraphStruct)
    # For prototyping assume homogeneous edges and remove the interface dim doubling
    v_jac_array = [[zeros(vdim, vdim), zeros(vdim, edim)] for (vdim, edim) in zip(gs.v_dims, Iterators.repeated(gs.e_dims[1] ÷ 2))]
    e_jac_array = [[zeros(edim, srcdim), zeros(edim, dstdim)] for (edim, srcdim, dstdim) in zip( Iterators.repeated(gs.e_dims[1] ÷ 2), gs.v_dims[gs.s_e], gs.v_dims[gs.d_e])] # homogene Netzwerke: v_src_dim = v_dst_dim = v_dim
    e_jac_product = [zeros(dim) for dim in Iterators.repeated(gs.e_dims[1] ÷ 2, gs.num_e)]
    JacGraphData1(v_jac_array, e_jac_array, e_jac_product)
end

# anscheinend braucht man diese arrays nicht mehr für JacGraphData1
#v_jac_array = [Array{Float64,2}(undef, dim, dim) for dim in v_dims]
#e_jac_array = [[zeros(dim, srcdim), zeros(dim, dstdim)] for (dim, srcdim, dstdim) in zip(v_dims, v_dims, v_dims)] # homogene Netzwerke: v_src_dim = v_dst_dim = v_dim
#e_jac_product =  [zeros(v_dims[1]) for i in 1:num_e] # e_dims

jac_graph_data_object = JacGraphData1(graph_structure)

#### NDJacVecOperator1

mutable struct NDJacVecOperator1{T, uType, tType, T1, T2, GD, JGD} <: DiffEqBase.AbstractDiffEqLinearOperator{T} # mutable da x, p, t geupdated werden
    x::uType
    p
    t::tType
    vertices!::T1
    edges!::T2
    graph_structure::GraphStruct
    graph_data::GD
    jac_graph_data::JGD
    parallel::Bool

    function NDJacVecOperator1{T}(x, p, t, vertices!, edges!, graph_structure, graph_data, jac_graph_data, parallel) where T
        new{T, typeof(x), typeof(t), typeof(vertices!), typeof(edges!), typeof(graph_data), typeof(jac_graph_data)}(x, p, t, vertices!, edges!, graph_structure, graph_data, jac_graph_data, parallel)
    end

    function NDJacVecOperator1(x, p, t, vertices!, edges!, graph_structure, graph_data, jac_graph_data, parallel)
        NDJacVecOperator1{eltype(x)}(x, p, t, vertices!, edges!, graph_structure, graph_data, jac_graph_data, parallel)
    end
end

x = similar(zeros(1), sum(v_dims))
p = nothing
t = 0.0
parallel = false

NDJacVecOp = NDJacVecOperator1(x, p, t, vertices!, edges!, graph_structure, graph_data, jac_graph_data_object, parallel) # 51 allocations


### get functions for update_coefficients

@inline get_src_edge_jacobian(jgd::JacGraphData1, i::Int) = jgd.e_jac_array[i][1]
@inline get_dst_edge_jacobian(jgd::JacGraphData1, i::Int) = jgd.e_jac_array[i][2]
@inline get_internal_vertex_jacobian(jgd::JacGraphData1, i::Int) = jgd.v_jac_array[i][1]
@inline get_aggregation_vertex_jacobian(jgd::JacGraphData1, i::Int) = jgd.v_jac_array[i][2]


function update_coefficients1!(Jac::NDJacVecOperator1, x, p, t)

    gs = Jac.graph_structure
    checkbounds_p(p, gs.num_v, gs.num_e)
    jgd = Jac.jac_graph_data

    for i in 1:gs.num_v
        @inbounds maybe_idx(Jac.vertices!, i).vertex_jacobian!(
          get_internal_vertex_jacobian(jgd, i),
          get_aggregation_vertex_jacobian(jgd, i),
          view(x, gs.v_idx[i]),
          p_v_idx(p, i),
          t)
    end

    for i in 1:gs.num_e
    @inbounds maybe_idx(Jac.edges!, i).edge_jacobian!(
              get_src_edge_jacobian(jgd, i),
              get_dst_edge_jacobian(jgd, i),
              view(x, gs.s_e_idx[i]),
              view(x, gs.d_e_idx[i]),
              p_e_idx(p, i),
              t)
      end

    Jac.x = x
    Jac.p = p
    Jac.t = t
end


x_test2 = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
dx_test = similar(x_test2)
#dx_test = similar(x_test2)

fill!(dx_test, zero(dx_test[1]))

p_test = nothing
t_test = 0.0
parallel = false

update_coefficients1!(NDJacVecOp, x_test2, p_test, t_test)

# testing single components of the new jac_vec_prod!

gs = NDJacVecOp.graph_structure
jgd = NDJacVecOp.jac_graph_data

get_dst_edges(graph_data, 4)


z = x_test2

for i in 1:gs.num_e
    # Store Edge_jac_src * v_src
    @inbounds mul!(jgd.e_jac_product[i], get_src_edge_jacobian(jgd, i), view(z, gs.s_e_idx[i]))
end

jgd.e_jac_product

for i in 1:gs.num_e
    # in-place Add Edge_jac_dst * v_dst
    # mul!(C,A,B,α,β): C = A B α + C β
    @inbounds mul!(jgd.e_jac_product[i],
         get_dst_edge_jacobian(jgd, i),
         view(z, gs.d_e_idx[i]), 1, 1) # α = 1, β = 1
end

jgd.e_jac_product

for i in 1:gs.num_e
    @inbounds mul!(view(dx_test, gs.d_e_idx[i]),
                   get_aggregation_vertex_jacobian(jgd, gs.d_e[i]),
                   jgd.e_jac_product[i], 1, 1)
end

dx_test

for i in 1:gs.num_e
    @inbounds mul!(view(dx_test, gs.s_e_idx[i]),
                   get_aggregation_vertex_jacobian(jgd, gs.s_e[i]),
                   jgd.e_jac_product[i], -1, 1) # note the -1
end

dx_test

# 2. for loop
for i in 1:gs.num_v
    @inbounds mul!(view(dx_test, gs.v_idx[i]),
                   get_internal_vertex_jacobian(jgd, i), view(z, gs.v_idx[i]), 1, 1)
end

dx_test



get_src_edge_jacobian(jgd, 1) * view(z, gs.s_e_idx[1])
get_src_edge_jacobian(jgd, 1)
get_dst_edge_jacobian(jgd, 1)
get_internal_vertex_jacobian(jgd, 1)
get_aggregation_vertex_jacobian(jgd, 1)

NDJacVecOp.jac_graph_data.v_jac_array

# starting forwarddiff


using ForwardDiff
using Test


# vertex_jacobian

function pseudo_diffusionvertex!(dv, v, e, p, t)
    dv[1] = e[1]
    dv[2] = e[1]
    nothing
end

dv = Array{Float64,2}(undef, 2, 1)
jac_vertex = zeros(2, 3)
v_jac_array_entry = [[zeros(vdim, vdim), zeros(vdim, edim)] for (vdim, edim) in zip(gs.v_dims, Iterators.repeated(gs.e_dims[1] ÷ 2))][1]


d = [1.0, 1.0, 1.0]
p = 1.0
t = 0.0


# richtige Matrix kommt raus
cfg_vertex = ForwardDiff.JacobianConfig((dv, x) -> pseudo_diffusionvertex!(dv, x[1:2], x[3], p, t), dv, d)
ForwardDiff.jacobian!(jac_vertex, (dv, x) -> pseudo_diffusionvertex!(dv, x[1:2], x[3], p, t), dv, d, cfg_vertex, Val{false}())

ForwardDiff.jacobian!(jac_vertex, (dv, x) -> pseudo_diffusionvertex!(dv, x[1:2], x[3], p, t), dv, d)


# extract vertex jacobians and put it into the right format
v_jac_array_entry[1] = jac_vertex[:, 1:2][:, :] # ohne das geht es hier auch
v_jac_array_entry[2] = jac_vertex[:, 3][:, :]

# bzw.
J_v = jac_vertex[:, 1:2][:, :]
J_e = jac_vertex[:, 3][:, :]

@test NDJacVecOp.jac_graph_data.v_jac_array[1] == v_jac_array_entry


# edge jacobian

function pseudo_diffusionedge!(e, v_s, v_d, p, t)
    e[1] = v_s[1] - v_d[1]
    nothing
end

dv6 = Array{Float64,2}(undef, 1, 1)
jac_edge = zeros(1, 4)
e_jac_array_entry = [[zeros(edim, srcdim), zeros(edim, dstdim)] for (edim, srcdim, dstdim) in zip( Iterators.repeated(gs.e_dims[1] ÷ 2), gs.v_dims[gs.s_e], gs.v_dims[gs.d_e])][1]
h6 = [1.0, 1.0, 1.0, 1.0]

cfg_edge = ForwardDiff.JacobianConfig((dv6, x) -> pseudo_diffusionedge!(dv6, x[1:2], x[3:4], p, t), dv6, h6)
ForwardDiff.jacobian!(jac_edge, (dv6, x) -> pseudo_diffusionedge!(dv6, x[1:2], x[3:4], p, t), dv6, h6, cfg_edge, Val{false}())


# extract edge jacobians and put it into the right format
e_jac_array_entry[1] = transpose(jac_edge[1:2][:, :])
e_jac_array_entry[2] = transpose(jac_edge[3:4][:, :])


# bzw.
J_src = transpose(jac_edge[1:2][:, :])
J_dst = transpose(jac_edge[3:4][:, :])

@test NDJacVecOp.jac_graph_data.e_jac_array[1] == e_jac_array_entry















# vertex jacobian - subarrays berechnen

pseudo_jac_v = zeros(2, 2)
pseudo_jac_e = zeros(2, 1)

e = [1.0 0.0]
v = [1.0 1.0]
b = [1.0, 1.0]
c = [1.0]

# Unterjacobians (sub arrays) einzeln berechnen
#cfg1 = ForwardDiff.JacobianConfig((dv, v) -> pseudo_diffusionvertex!(dv, v, e, p, t), dv, b)
ForwardDiff.jacobian!(pseudo_jac_v, (dv, v) -> pseudo_diffusionvertex!(dv, v, e, p, t), dv, b)

#cfg11 = ForwardDiff.JacobianConfig((dv, e) -> pseudo_diffusionvertex!(dv, v, e, p, t), dv, c)
ForwardDiff.jacobian!(pseudo_jac_e, (dv, e) -> pseudo_diffusionvertex!(dv, v, e, p, t), dv, c)




# falscher vertex jacobian
cfg2 = ForwardDiff.JacobianConfig((dv, v) -> pseudo_diffusionvertex!(dv, v[1], v[2], p, t), dv, d)
ForwardDiff.jacobian!(pseudo_jac_all, (dv, x) -> pseudo_diffusionvertex!(dv, x[1], x[2], p, t), dv, d, cfg2, Val{false}())



# trying other functions

function pseudo_diffusionvertex2!(dv, v, e, p, t)
    dv[1] = e[1] + v[1]
    dv[2] = e[1] + v[1]
    nothing
end

cfg4 = ForwardDiff.JacobianConfig((dv, x) -> pseudo_diffusionvertex2!(dv, x[1:2], x[3], p, t), dv, d)
ForwardDiff.jacobian!(pseudo_jac_all, (dv, x) -> pseudo_diffusionvertex2!(dv, x[1:2], x[3], p, t), dv, d, cfg4, Val{false}())

function pseudo_diffusionvertex3!(dv, v, e, p, t)
    dv[1] = e[1] + v[1]
    dv[2] = e[1] + v[2]
    nothing
end

cfg5 = ForwardDiff.JacobianConfig((dv, x) -> pseudo_diffusionvertex3!(dv, x[1:2], x[3], p, t), dv, d)
ForwardDiff.jacobian!(pseudo_jac_all, (dv, x) -> pseudo_diffusionvertex3!(dv, x[1:2], x[3], p, t), dv, d, cfg5, Val{false}())


function pseudo_diffusionvertex4!(dv, v, e, p, t)
    dv[1] = e[1] + v[2] + v[1]
#    dv[2] = e[1] + v[2]
    nothing
end

pseudo_jac_all55 = zeros(1, 3)
dv55 = Array{Float64,2}(undef, 1, 1)

cfg55 = ForwardDiff.JacobianConfig((dv55, x) -> pseudo_diffusionvertex4!(dv55, x[1:2], x[3], p, t), dv55, d)
ForwardDiff.jacobian!(pseudo_jac_all55, (dv55, x) -> pseudo_diffusionvertex4!(dv55, x[1:2], x[3], p, t), dv55, d, cfg55, Val{false}())

















# outsourced

function attempt(dx, x, p, t)
  dx[1] = x[1]^2+x[2]^3-1
  dx[2] = x[1]^4 - x[2]^4 + x[1]*x[2]
  nothing
end


a = [2.0, 1.0] # vermutlich die aktuellen x-werte bzw. die Anfangsbedingungen für den jeweiligen Schritt


dx = Array{Float64,2}(undef, 2, 1)
results = zeros(2, 2)

p = 1.0
t = 0.0


cfg = ForwardDiff.JacobianConfig((y, x) -> attempt(y, x, p, t), dx, a)

ForwardDiff.jacobian!(results, (y, x) -> attempt(y, x, p, t), dx, a, cfg, Val{false}())

results
