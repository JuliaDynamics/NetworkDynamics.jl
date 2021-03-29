using DifferentialEquations
using Reexport
using LinearAlgebra
using NetworkDynamics
using LightGraphs
using DiffEqBase
using DiffEqOperators

#include("/Users/lalalulu/Desktop/PIK/github/NetworkDynamics.jl/Jacobians/jac_structs.jl")
#@reexport using .jac_structs

N = 4
k = 2
g = barabasi_albert(N, k)

function diffusionedge!(e, v_s, v_d, p, t)
    # usually e, v_s, v_d are arrays, hence we use the broadcasting operator .
    #e .= v_s .- v_d
    e[1] = v_s[1] - v_d[1]
    nothing
end

function diffusionvertex!(dv, v, edges, p, t)
    # usually dv, v, edges are arrays, hence we use the broadcasting operator .
    dv[1] = 0.
    #dv .= 0.
    dv[2] = 1.
    for e in edges
        dv[1] += e[1]
    end
    nothing
end

nd_diffusion_vertex = ODEVertex(f! = diffusionvertex!, dim = 2)

nd_diffusion_edge = StaticEdge(f! = diffusionedge!, dim = 1)

nd = network_dynamics(nd_diffusion_vertex, nd_diffusion_edge, g)

x0 = randn(2N) # random initial conditions
ode_prob = ODEProblem(nd, x0, (0., 4.))
sol = solve(ode_prob, Tsit5());

# collecting the graph infos

nd.f
graph_structure_ = nd.f.graph_structure
graph_data_ = nd.f.graph_data
graph_ = nd.f.graph

# build the jacobian graph data

struct JacGraphData
    v_jac_array::Array{Array{Float64, 2}, 1}
    e_jac_array::Array{Array{Array{Float64, 2}, 1}, 1}
    e_jac_product::Array{Float64, 2}
end

v_dims = nd.f.graph_structure.v_dims
e_dims = nd.f.graph_structure.e_dims
num_e = nd.f.graph_structure.num_e

v_jac_array = [Array{Float64,2}(undef, dim, dim) for dim in v_dims]
e_jac_array = [[zeros(dim, srcdim), zeros(dim, dstdim)] for (dim, srcdim, dstdim) in zip(e_dims, v_dims, v_dims)] # homogene Netzwerke: v_src_dim = v_dst_dim = v_dim
e_jac_product =  zeros(e_dims[1], num_e) # Annahme: homogene edges

jac_graph_data_object = JacGraphData(v_jac_array, e_jac_array, e_jac_product)

mutable struct NDJacVecOperator{uType, P, tType, G, GD, JGD, T} <: DiffEqBase.AbstractDiffEqLinearOperator{T} # mutable da x, p, t geupdated werden
    x::uType
    p
    t::tType
    graph::G
    graph_structure::GraphStruct
    graph_data::GD
    jac_graph_data::JGD
    parallel::Bool
end

x = similar(zeros(1), sum(v_dims))
p = nothing
t = 0.0
parallel = false

NDJacVecOperator_object = NDJacVecOperator(x, p, t, g, graph_structure_, graph_data_, jac_graph_data_object, parallel)


#### modified mutable struct

mutable struct NDJacVecOperator6{uType, tType, G, GS, GD, JGD} # mutable da x, p, t geupdated werden
    x::uType
    p
    t::tType
    graph::G
    graph_structure::GS
    graph_data::GD
    jac_graph_data::JGD
    parallel::Bool
end

NDJacVecOperator6_object = NDJacVecOperator6(x, p, t, g, graph_structure_, graph_data_, jac_graph_data_object, parallel)


mutable struct NDJacVecOperator7{uType, tType, G, GS, GD, JGD, T} <: DiffEqBase.AbstractDiffEqLinearOperator{T} # mutable da x, p, t geupdated werden
    x::uType
    p
    t::tType
    graph::G
    graph_structure::GS
    graph_data::GD
    jac_graph_data::JGD
    parallel::Bool
end

NDJacVecOperator7_object = NDJacVecOperator7(x, p, t, g, graph_structure_, graph_data_, jac_graph_data_object, parallel)

mutable struct NDJacVecOperator8{T, uType, tType, G, GS, GD, JGD} <: DiffEqBase.AbstractDiffEqLinearOperator{T} # mutable da x, p, t geupdated werden
    x::uType
    p
    t::tType
    graph::G
    graph_structure::GS
    graph_data::GD
    jac_graph_data::JGD
    parallel::Bool
end

NDJacVecOperator8_object = NDJacVecOperator8(x, p, t, g, graph_structure_, graph_data_, jac_graph_data_object, parallel)

mutable struct NDJacVecOperator9{T, uType, tType, G, GD, JGD} <: DiffEqBase.AbstractDiffEqLinearOperator{T} # mutable da x, p, t geupdated werden
    x::uType
    p
    t::tType
    graph::G
    graph_structure::GraphStruct
    graph_data::GD
    jac_graph_data::JGD
    parallel::Bool
end

NDJacVecOperator9_object = NDJacVecOperator9(x, p, t, g, graph_structure_, graph_data_, jac_graph_data_object, parallel)

mutable struct NDJacVecOperator10{T, uType, tType, G, GD, JGD} <: DiffEqBase.AbstractDiffEqLinearOperator{T} # mutable da x, p, t geupdated werden
    x::uType
    p
    t::tType
    graph::G
    graph_structure::GraphStruct
    graph_data::GD
    jac_graph_data::JGD
    parallel::Bool

    function NDJacVecOperator10{T}(x, p, t, graph, graph_structure, graph_data, jac_graph_data, parallel) where T
        new{T,typeof(x),typeof(t),typeof(graph),typeof(graph_data),typeof(jac_graph_data)}(x, p, t, graph, graph_structure, graph_data, jac_graph_data, parallel)
    end

    function NDJacVecOperator10(x, p, t, graph, graph_structure, graph_data, jac_graph_data, parallel)
        NDJacVecOperator10{eltype(x)}(x, p, t, graph, graph_structure, graph_data, jac_graph_data, parallel)
    end
end

NDJacVecOperator10_object = NDJacVecOperator10(x, p, t, g, graph_structure_, graph_data_, jac_graph_data_object, parallel)


### get functions for update_coefficients

@inline get_src_edge_jacobian(jgd::JacGraphData, i::Int) = jgd.e_jac_array[i][1]
@inline get_dst_edge_jacobian(jgd::JacGraphData, i::Int) = jgd.e_jac_array[i][2]
@inline get_vertex_jacobian(jgd::JacGraphData, i::Int) = jgd.v_jac_array[i]

### Vertex, Edge functions for update_coefficients

function VertexJacobian!(J::AbstractMatrix, v, p, t)
    #J = internal_jacobian(J, v, p, t)
    J[1, 1] = 1.0
    J[1, 2] = 0.0
    J[2, 1] = 1.0
    J[2, 2] = 0.0
end

function EdgeJacobian!(J_s::AbstractMatrix, J_d::AbstractMatrix, v_s, v_d, p, t)
   #J_s = source_jacobian(v_s, v_d, p, t)
   #J_d = dest_jacobian(v_s, v_d, p, t)
   J_s[1, 1] = 1.0
   J_s[1, 2] = 0.0
   J_d[1, 1] = 1.0
   J_d[1, 2] = 0.0
end



function update_coefficients!(Jac::NDJacVecOperator, x, p, t)

    gs = Jac.graph_structure
    checkbounds_p(p, gs.num_v, gs.num_e)
    gd = prep_gd(x, x, Jac.graph_data, Jac.graph_structure)
    jgd = Jac.jac_graph_data

    for i in 1:gs.num_v
        vertices![i].VertexJacobian!(
          get_vertex_jacobian(jgd, i),
          get_vertex(gd, i),
          p_v_idx(p, i),
          t)
    end

    for i in 1:gs.num_e
          edges![i].EdgeJacobian!(
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

function (Jac::NDJacVecOperator)(x, p, t) # auch Number bei t?
    update_coefficients!(Jac, x, p, t)
    Jac*x
end

function (Jac::NDJacVecOperator)(dx, x, p, t::Number)
    update_coefficients!(Jac, x, p, t)
    mul!(dx, Jac, x)
end

Base.:*(Jac::NDJacVecOperator, z::AbstractVector) = jac_vec_prod(Jac, z)

function LinearAlgebra.mul!(dx::AbstractVector, Jac::NDJacVecOperator, z::AbstractVector)
    jac_vec_prod!(dx, Jac, z)
end

function jac_vec_prod(Jac::NDJacVecOperator, z)

    gs = Jac.graph_structure
    checkbounds_p(p, gs.num_v, gs.num_e)
    gd = prep_gd(x, x, Jac.graph_data, Jac.graph_structure)
    jgd = Jac.jac_graph_data

    for i in 1:gs.num_e
        e_jac_product[i] .= get_src_edge_jacobian(jgd, i) * view(z, get_src_indices(i)) + get_dst_edge_jacobian(jgd, i) * view(z, get_dst_indices(i))
    end
    ## neues array wird erstellt und returnt

    dx = zeros(gs.v_dims, gs.num_v)

    for i in 1:gs.num_v
        view(dx, gs.v_idx[i]) .= get_vertex_jacobian(gd, i) * view(z, gs.v_idx[i])
        dx .+= sum(e_jac_product[i])
    end
    return dx
end

function jac_vec_prod!(dx, Jac::NDJacVecOperator, z)

    gs = Jac.graph_structure
    checkbounds_p(p, gs.num_v, gs.num_e)
    gd = prep_gd(x, x, Jac.graph_data, Jac.graph_structure)

    for i in 1:gs.num_e
        e_jac_product[i] .= get_src_edge_jacobian(gd, i) * view(z, get_src_indices(i)) + get_dst_edge_jacobian(gd, i) * view(z, get_dst_indices(i))
    end

    for i in 1:gs.num_v
        view(dx, get_vertex_indices(i)) .= get_vertex_jacobian(gd, i) * view(z, get_vertex_indices(i))
        dx .+= sum(e_jac_product[i])
    end
end



### test stuff

#v_jac_array = Array{Array{Float64, 2}, 1}(undef, 2)
#v_jac_array[1] = [1.0 2.0]
#v_jac_array[2] = [3.0 4.0]
#println(v_jac_array)

#e_jac_array = Array{Array{Array{Float64, 2}, 1}, 1}(undef, 1)
#e_jac_array[1] = [[1.0 2.0]]
#e_jac_array[2] = [[3.0 4.0]]
#print(e_jac_array)

#e_jac_product = Array{Float64, 2.0}
#e_jac_product[1] = [1.0]

#=mutable struct NDJacVecOperator1{uType, tType, JGD} # mutable da x, p, t geupdated werden
    x::uType
    p
    t::tType
    jac_graph_data::JGD
    parallel::Bool
end

NDJacVecOperator1_object = NDJacVecOperator1(x, p, t, jac_graph_data_object, parallel)


mutable struct NDJacVecOperator2{uType, tType, G, JGD} # mutable da x, p, t geupdated werden
    x::uType
    p
    t::tType
    graph::G
    jac_graph_data::JGD
    parallel::Bool
end

NDJacVecOperator2_object = NDJacVecOperator2(x, p, t, graph_, jac_graph_data_object, parallel)

mutable struct NDJacVecOperator3{uType, tType, JGD, T} <: DiffEqBase.AbstractDiffEqLinearOperator{T} # mutable da x, p, t geupdated werden
    x::uType
    p
    t::tType
    #graph::G
    #graph_structure::GraphStruct
    #graph_data::GD
    jac_graph_data::JGD
    parallel::Bool
end

NDJacVecOperator3_object = NDJacVecOperator3(x, p, t, jac_graph_data_object, parallel)

mutable struct NDJacVecOperator4{uType, tType, GD, JGD} # mutable da x, p, t geupdated werden
    x::uType
    p
    t::tType
    #graph::G
    #graph_structure::GraphStruct
    graph_data::GD
    jac_graph_data::JGD
    parallel::Bool
end

NDJacVecOperator4_object = NDJacVecOperator4(x, p, t, graph_data_, jac_graph_data_object, parallel)

mutable struct NDJacVecOperator5{uType, tType, GS, GD, JGD} # mutable da x, p, t geupdated werden
    x::uType
    p
    t::tType
    #graph::G
    graph_structure::GS
    graph_data::GD
    jac_graph_data::JGD
    parallel::Bool
end

NDJacVecOperator5_object = NDJacVecOperator5(x, p, t, graph_data_, graph_structure_, jac_graph_data_object, parallel)=#
