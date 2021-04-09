

module Jacobians

using ..NetworkStructures
using ..Utilities

using DifferentialEquations
using Reexport
using LinearAlgebra
using DiffEqBase
import DiffEqBase.update_coefficients!
export JacGraphData, NDJacVecOperator
#export maybe_idx, p_v_idx, p_e_idx, @nd_threads, warn_parallel, checkbounds_p, prep_gd


mutable struct JacGraphDataBuffer{Tvj, Tej, Tep}
    v_Jac_array::Tvj
    e_Jac_array::Tej
    e_Jac_product_array::Tep
end

struct JacGraphData{JGDB}
    jgdb::JGDB
    v_jac_array::Array{Array{Float64, 2}, 1}
    e_jac_array::Array{Array{Array{Float64, 2}, 1}, 1}
    e_jac_product::Array{Float64, 2}
end

function JacGraphData(v_Jac_array::Tvj, e_Jac_array::Tej, e_Jac_product_array::Tep, gs::GraphStruct) where {Tvj, Tej, Tep}
    jgdb = JacGraphDataBuffer{Tvj, Tej, Tep}(v_Jac_array, e_Jac_array, e_Jac_product_array)
    JGDB = typeof(jgdb)

    v_jac = [Array{Float64,2}(undef, dim, dim) for dim in gs.v_dims]
    e_jac = [[zeros(dim, srcdim), zeros(dim, dstdim)] for (dim, srcdim, dstdim) in zip(gs.e_dims, gs.v_dims, gs.v_dims)] # homogene Netzwerke: v_src_dim = v_dst_dim = v_dim
    e_jac_product =  zeros(gs.e_dims[1], gs.num_e) # Annahme: homogene edges

    JacGraphData{JGDB}(jgdb, v_jac, e_jac, e_jac_product)
end

@inline get_src_edge_jacobian(jgd::JacGraphData, i::Int) = jgd.e_jac_array[i][1]

@inline get_dst_edge_jacobian(jgd::JacGraphData, i::Int) = jgd.e_jac_array[i][2]

@inline get_vertex_jacobian(jgd::JacGraphData, i::Int) = jgd.v_jac_array[i]

# schon in network_dynamics Funktion vorhanden

#v_jac_array = [Array{Float64,2}(undef, dim, dim) for dim in v_dims]
#e_jac_array = [[zeros(dim, srcdim), zeros(dim, dstdim)] for (dim, srcdim, dstdim) in zip(e_dims, v_dims, v_dims)] # homogene Netzwerke: v_src_dim = v_dst_dim = v_dim
#e_jac_product =  zeros(num_e, e_dims[1]) # Annahme: homogene edges


mutable struct NDJacVecOperator{T, uType, tType, G, GD, JGD} <: DiffEqBase.AbstractDiffEqLinearOperator{T} # mutable da x, p, t geupdated werden
    x::uType
    p
    t::tType
    graph::G
    graph_structure::GraphStruct
    graph_data::GD
    jac_graph_data::JGD
    parallel::Bool

    function NDJacVecOperator{T}(x, p, t, graph, graph_structure, graph_data, jac_graph_data, parallel) where T
        new{T,typeof(x),typeof(t),typeof(graph),typeof(graph_data),typeof(jac_graph_data)}(x, p, t, graph, graph_structure, graph_data, jac_graph_data, parallel)
    end

    function NDJacVecOperator(x, p, t, graph, graph_structure, graph_data, jac_graph_data, parallel)
        NDJacVecOperator{eltype(x)}(x, p, t, graph, graph_structure, graph_data, jac_graph_data, parallel)
    end
end

# for this we will need prep_gc and inbounds_p -> reexport all nd files (nd_DDE_static, nd_ODE_ODE, nd_ODE_Static)

function update_coefficients!(Jac::NDJacVecOperator, x, p, t)

    gs = Jac.graph_structure
    checkbounds_p(p, gs.num_v, gs.num_e)
    gd = prep_gd(x, x, Jac.graph_data, Jac.graph_structure)
    jgd = Jac.jac_graph_data

    for i in 1:gs.num_v
        #vertices![i].vertex_jacobian!(
        maybe_idx(vertices!, i).vertex_jacobian!(
          get_vertex_jacobian(jgd, i),
          get_vertex(gd, i),
          p_v_idx(p, i),
          t)
    end

    for i in 1:gs.num_e
          #edges![i].edge_jacobian!(
          maybe_idx(edges!, i).edge_jacobian!(
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

# final functions for NDJacVecOperator

function jac_vec_prod(Jac::NDJacVecOperator, z)

    gs = Jac.graph_structure
    checkbounds_p(p, gs.num_v, gs.num_e)
    gd = prep_gd(x, x, Jac.graph_data, Jac.graph_structure)
    jgd = Jac.jac_graph_data

    for i in 1:gs.num_e
        #e_jac_product[i] .= get_src_edge_jacobian(gd, i) * view(z, get_src_indices(i)) + get_dst_edge_jacobian(gd, i) * view(z, get_dst_indices(i))
        e_jac_product[i, :] .= get_src_edge_jacobian(jgd, i) * view(z, gs.s_e_idx[i]) + get_dst_edge_jacobian(jgd, i) * view(z, gs.d_e_idx[i])
        #println(e_jac_product[i, :])
    end

    ## neues array wird erstellt und returned
    dx = zeros(gs.v_dims[1], gs.num_v)

    for i in 1:gs.num_v
        #view(dx, get_vertex_indices(i)) .= get_vertex_jacobian(gd, i) * view(z, get_vertex_indices(i))
        view(dx, gs.v_idx[i]) .= get_vertex_jacobian(jgd, i) * view(z, gs.v_idx[i])
        #dx .+= sum(e_jac_product[i])
        dx .+= sum(e_jac_product[i, :])
    end
    return dx
end

function jac_vec_prod!(dx, Jac::NDJacVecOperator, z)

    gs = Jac.graph_structure
    checkbounds_p(p, gs.num_v, gs.num_e)
    gd = prep_gd(x, x, Jac.graph_data, Jac.graph_structure)
    jgd = Jac.jac_graph_data

    for i in 1:gs.num_e
        #e_jac_product[i] .= get_src_edge_jacobian(gd, i) * view(z, get_src_indices(i)) + get_dst_edge_jacobian(gd, i) * view(z, get_dst_indices(i))
        e_jac_product[i, :] .= get_src_edge_jacobian(jgd, i) * view(z, gs.s_e_idx[i]) + get_dst_edge_jacobian(jgd, i) * view(z, gs.d_e_idx[i])
        #println(e_jac_product[i, :])
    end

    for i in 1:gs.num_v
        #view(dx, get_vertex_indices(i)) .= get_vertex_jacobian(gd, i) * view(z, get_vertex_indices(i))
        view(dx, gs.v_idx[i]) .= get_vertex_jacobian(jgd, i) * view(z, gs.v_idx[i])
        #dx .+= sum(e_jac_product[i])
        dx .+= sum(e_jac_product[i, :])
    end
end


# functions for callable structs

Base.:*(Jac::NDJacVecOperator, z::AbstractVector) = jac_vec_prod(Jac, z)

function LinearAlgebra.mul!(dx::AbstractVector, Jac::NDJacVecOperator, z::AbstractVector)
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
