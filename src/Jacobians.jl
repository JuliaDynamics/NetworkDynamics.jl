

module Jacobians

using ..NetworkStructures
using ..Utilities
using ..nd_ODE_Static_mod

using DifferentialEquations
using Reexport
using LinearAlgebra
using DiffEqBase
import DiffEqBase.update_coefficients!
export JacGraphData, NDJacVecOperator

struct JacGraphData
    v_jac_array::Array{Array{Float64, 2}, 1} # contains the jacobians for each vertex
    e_jac_array::Array{Array{Array{Float64, 2}, 1}, 1} # contains the jacobians for each edge
    e_jac_product::Array{Array{Float64, 1}, 1} # is needed later in jac_vec_prod(!) as a storage for the products of edge jacobians and vectors z
end


function JacGraphData(v_jac_array, e_jac_array, e_jac_product_array, gs::GraphStruct)
    v_jac_array = [Array{Float64,2}(undef, dim, dim) for dim in gs.v_dims]
    e_jac_array = [[zeros(dim, srcdim), zeros(dim, dstdim)] for (dim, srcdim, dstdim) in zip(gs.e_dims, gs.v_dims, gs.v_dims)] # homogene Netzwerke: v_src_dim = v_dst_dim = v_dim
    e_jac_product = [zeros(gs.e_dims[1]) for i in 1:gs.num_e]
    JacGraphData(v_jac_array, e_jac_array, e_jac_product)
end


@inline get_src_edge_jacobian(jgd::JacGraphData, i::Int) = jgd.e_jac_array[i][1]

@inline get_dst_edge_jacobian(jgd::JacGraphData, i::Int) = jgd.e_jac_array[i][2]

@inline get_vertex_jacobian(jgd::JacGraphData, i::Int) = jgd.v_jac_array[i]


mutable struct NDJacVecOperator{T, uType, tType, T1, T2, G, GD, JGD} <: DiffEqBase.AbstractDiffEqLinearOperator{T} # mutable da x, p, t geupdated werden
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

    function NDJacVecOperator{T}(x, p, t, vertices!, edges!, graph, graph_structure, graph_data, jac_graph_data, parallel) where T
        new{T, typeof(x), typeof(t), typeof(vertices!), typeof(edges!),typeof(graph), typeof(graph_data), typeof(jac_graph_data)}(x, p, t, vertices!, edges!, graph, graph_structure, graph_data, jac_graph_data, parallel)
    end

    function NDJacVecOperator(x, p, t, vertices!, edges!, graph, graph_structure, graph_data, jac_graph_data, parallel)
        NDJacVecOperator{eltype(x)}(x, p, t, vertices!, edges!, graph, graph_structure, graph_data, jac_graph_data, parallel)
    end
end

# for this we will need prep_gc and inbounds_p -> reexport all nd files (nd_DDE_static, nd_ODE_ODE, nd_ODE_Static)

function update_coefficients!(Jac::NDJacVecOperator, x, p, t)

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

# functions for NDJacVecOperator: both syntaxes must be taken into account: Jac, z and dx, Jac, z

function jac_vec_prod(Jac::NDJacVecOperator, z)

    gs = Jac.graph_structure
    p = Jac.p
    x = Jac.x
    checkbounds_p(p, gs.num_v, gs.num_e)
    gd = prep_gd(x, x, Jac.graph_data, Jac.graph_structure)
    jgd = Jac.jac_graph_data

    # first for loop that considers the mutliplication of each edge jacobians with the corresponding component of z
    for i in 1:gs.num_e
        jgd.e_jac_product[i] .= get_src_edge_jacobian(jgd, i) * view(z, gs.s_e_idx[i]) + get_dst_edge_jacobian(jgd, i) * view(z, gs.d_e_idx[i])
    end

    # in this function there is no dx in which the Jacobian can be stored, so an extra array must be created and returned
    dx = zeros(gs.v_dims[1], gs.num_v)

    # second for loop in which the multiplication of vertex jacobian and the corresponding component of z is done with addition of the e_jac_product to dx
    for i in 1:gs.num_v
        view(dx, gs.v_idx[i]) .= get_vertex_jacobian(jgd, i) * view(z, gs.v_idx[i])
        dx .+= sum(jgd.e_jac_product[i])
    end
    return dx
end

function jac_vec_prod!(dx, Jac::NDJacVecOperator, z)

    gs = Jac.graph_structure
    p = Jac.p
    x = Jac.x
    checkbounds_p(p, gs.num_v, gs.num_e)
    gd = prep_gd(x, x, Jac.graph_data, Jac.graph_structure)
    jgd = Jac.jac_graph_data

    for i in 1:gs.num_e
        jgd.e_jac_product[i] .= get_src_edge_jacobian(jgd, i) * view(z, gs.s_e_idx[i]) + get_dst_edge_jacobian(jgd, i) * view(z, gs.d_e_idx[i])
    end

    for i in 1:gs.num_v
        view(dx, gs.v_idx[i]) .= get_vertex_jacobian(jgd, i) * view(z, gs.v_idx[i])
        dx .+= sum(jgd.e_jac_product[i])
    end
end


# functions for callable structs at the end of this module

Base.:*(Jac::NDJacVecOperator, z::AbstractVector) = jac_vec_prod(Jac, z)

function LinearAlgebra.mul!(dx::Matrix{Float64}, Jac::NDJacVecOperator, z::AbstractVector)
    jac_vec_prod!(dx, Jac, z)
end

# callable structs

function (Jac::NDJacVecOperator)(x, p, t) # auch Number bei t?
    update_coefficients!(Jac, x, p, t)
    Jac*x
end

function (Jac::NDJacVecOperator)(dx, x, p, t::Number)
    update_coefficients!(Jac, x, p, t)
    mul!(dx, Jac, x)
end


end # module
