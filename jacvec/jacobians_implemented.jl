using DifferentialEquations
using Reexport
using LinearAlgebra
using NetworkDynamics
using LightGraphs
using DiffEqBase
using BenchmarkTools
using ForwardDiff


N = 4

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


function pseudo_diffusionvertex!(dv, v, e, p, t)
#    dv[1] = v[1]*v[1]*e[1]
#    dv[2] = v[2]*v[2]*e[1]
    dv[1] = e[1]
    dv[2] = e[1]
    nothing
end


@Base.propagate_inbounds function jac_vertex!(J_v::AbstractMatrix, J_e::AbstractMatrix, dv_dim, v, e, p, t)

    dv = Array{Float64,2}(undef, dv_dim, 1)
    ForwardDiff.jacobian!(J_v, (dv, v) -> pseudo_diffusionvertex!(dv, v, e, p, t), dv, v)
    ForwardDiff.jacobian!(J_e, (dv, e) -> pseudo_diffusionvertex!(dv, v, e, p, t), dv, e)
end


function jac_edge!(J_s::AbstractMatrix, J_d::AbstractMatrix, e_dim, v_s, v_d, p, t)

    e = Array{Float64,2}(undef, e_dim, 1)
    ForwardDiff.jacobian!(J_s, (e, v_s) -> diffusionedge!(e, v_s, v_d, p, t), e, v_s)
    ForwardDiff.jacobian!(J_d, (e, v_d) -> diffusionedge!(e, v_s, v_d, p, t), e, v_d)
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

NDJacVecOp.graph_data.dst_edges

sum(get_dst_edges(NDJacVecOp.graph_data, 1))

view(x, NDJacVecOp.graph_structure.v_idx[1])
get_vertex(NDJacVecOp.graph_data, 1)
get_edge(NDJacVecOp.graph_data, 1)

NDJacVecOp.vertices!.dim
NDJacVecOp.edges!.f!.dim


function update_coefficients1!(Jac::NDJacVecOperator1, x, p, t)

    gs = Jac.graph_structure
    checkbounds_p(p, gs.num_v, gs.num_e)
    gd = prep_gd(x, x, Jac.graph_data, Jac.graph_structure)
    jgd = Jac.jac_graph_data

    for i in 1:gs.num_v
        @inbounds maybe_idx(Jac.vertices!, i).vertex_jacobian!(
          get_internal_vertex_jacobian(jgd, i),
          get_aggregation_vertex_jacobian(jgd, i),
          Jac.vertices!.dim, # dimension for empty array dv for jacobian calculation
          get_vertex(gd, i),
          #view(x, gs.v_idx[i]),
          sum(get_dst_edges(gd, i)),
          p_v_idx(p, i),
          t)
          println(sum(get_dst_edges(gd, i)))
    end

    for i in 1:gs.num_e
    @inbounds maybe_idx(Jac.edges!, i).edge_jacobian!(
              get_src_edge_jacobian(jgd, i),
              get_dst_edge_jacobian(jgd, i),
              #view(x, gs.s_e_idx[i]),
              #view(x, gs.d_e_idx[i]),
              Jac.edges!.f!.dim, # dimension for empty array e for jacobian calculation
              get_dst_vertex(gd, i),
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
#dx_test = similar(x_test2)

fill!(dx_test, zero(dx_test[1]))

p_test = nothing
t_test = 0.0
parallel = false

update_coefficients1!(NDJacVecOp, x_test2, p_test, t_test)

# testing single components of the new jac_vec_prod!

gs = NDJacVecOp.graph_structure
jgd = NDJacVecOp.jac_graph_data



function test(a, b, c)
    a = 1
    b = 2
    c = a+b
    return nothing
end

function test2(test)
    test = @inline (a, b, c) -> begin
        d = c
        println(d)
    end
    test
    return nothing
end

test2(test)
