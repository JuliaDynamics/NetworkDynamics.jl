using DifferentialEquations
using Reexport
using LinearAlgebra

include("/Users/lalalulu/Desktop/PIK/github/NetworkDynamics.jl/Jacobians/jac_structs.jl")
@reexport using .jac_structs

mutable struct NDJacVecOperator1{uType, P, tType, G, GD, JGD, T} <: DiffEqBase.AbstractDiffEqLinearOperator{T} # mutable da x, p, t geupdated werden
    x::uType
    p::P
    t::tType
    graph::G
    graph_structure::GraphStruct
    graph_data::GD
    jac_graph_data::JGD # - dann muss oben noch JGD hin
    parallel::Bool
end


function update_coefficients!(Jac::NDJacVecOperator1, x, p, t)

    gs = Jac.graph_structure
    checkbounds_p(p, gs.num_v, gs.num_e)
    gd = prep_gd(x, x, Jac.graph_data, Jac.graph_structure)

    for i in 1:gs.num_v
        vertices![i].VertexJacobian!(
          get_vertex_jacobian(gd, i),
          get_vertex(gd, i),
          p_v_idx(p, i),
          t)
    end

    for i in 1:gs.num_e
          edges![i].EdgeJacobian!(
              get_src_edge_jacobian(gd, i),
              get_dst_edge_jacobian(gd, i),
              get_src_vertex(gd, i),
              get_dst_vertex(gd, i),
              p_e_idx(p, i),
              t)
      end

    Jac.x = x
    Jac.p = p
    Jac.t = t
end

function (Jac::NDJacVecOperator1)(x, p, t) # auch Number bei t?
    update_coefficients!(Jac, x, p, t)
    Jac*x
end

function (Jac::NDJacVecOperator1)(dx, x, p, t::Number)
    update_coefficients!(Jac, x, p, t)
    mul!(dx, Jac, x)
end

Base.:*(Jac::NDJacVecOperator1, z::AbstractVector) = jac_vec_prod(Jac, z)

function LinearAlgebra.mul!(dx::AbstractVector, Jac::NDJacVecOperator1, z::AbstractVector)
    jac_vec_prod!(dx, Jac, z)
end

function jac_vec_prod(Jac::NDJacVecOperator1, z)

    gs = Jac.graph_structure
    checkbounds_p(p, gs.num_v, gs.num_e)
    gd = prep_gd(x, x, Jac.graph_data, Jac.graph_structure)

    for i in 1:gs.num_e
        e_jac_product[i] .= get_src_edge_jacobian(gd, i) * view(z, get_src_indices(i)) + get_dst_edge_jacobian(gd, i) * view(z, get_dst_indices(i))
    end
    ## neues array wird erstellt und returnt

    dx = zeros(gs.v_dims, gs.num_v) # ?

    for i in 1:gs.num_v
        view(dx, get_vertex_indices(i)) = get_vertex_jacobian(gd, i) * view(z, get_vertex_indices(i))
        dx .+= sum(e_jac_product[i])
    end
    return dx
end

function jac_vec_prod!(dx, Jac::NDJacVecOperator1, z)

    gs = Jac.graph_structure
    checkbounds_p(p, gs.num_v, gs.num_e)
    gd = prep_gd(x, x, Jac.graph_data, Jac.graph_structure)

    for i in 1:gs.num_e
        e_jac_product[i] .= get_src_edge_jacobian(gd, i) * view(z, get_src_indices(i)) + get_dst_edge_jacobian(gd, i) * view(z, get_dst_indices(i))
    end

    for i in 1:gs.num_v
        view(dx, get_vertex_indices(i)) = get_vertex_jacobian(gd, i) * view(z, get_vertex_indices(i))
        dx .+= sum(e_jac_product[i])
    end
end
