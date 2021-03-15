#nd_jac = network_dynamics(nd_jac_vertex, nd_jac_edge, g)

using ForwardDiff
#using DiffEqOperators
using NetworkDynamics


### TO DO: wo werden diese Fkt integriert?
# VertexFunction and EdgeFunction will have new fields for Jacobians.
function VertexJacobian!(J::AbstractMatrix, v, p, t)
    J = internal_jacobian(J, v, p, t)
end

function EdgeJacobian!(J_s::AbstractMatrix, J_d::AbstractMatrix, v_s, v_d, p, t)
   J_s = source_jacobian(v_s, v_d, p, t)
   J_d = dest_jacobian(v_s, v_d, p, t)
end

struct NDJacVecOperator{uType, P, tType,G,GraphStruct,GD,Bool}# <: DiffEqBase.AbstractDiffEqLinearOperator{T}
    x::uType # Punkt u, für welchen das jvp berechnet wird
    p::P
    t::tType
    graph::G
    graph_structure::GraphStruct
    graph_data::GD
    parallel::Bool
end

### TO DO: in NetworkStructures.jl
# field names von graph data

v_jac_array = [Array{Float64,2}(undef, dim, dim) for dim in v_dims]
# For edges there is another array layer in the buffer since each edge has two Jacobians
e_jac_array = [[zeros(dim, srcdim), zeros(dim, dstdim)] for (dim, srcdim, dstdim) in zip(e_dims, v_src_dims, v_dst_dims)]

e_jac_product = [zeros(e_dims) , zeros(e_dims)]

### TO DO: in NetworkStructures.jl


@inline function get_src_edge_jacobian(gd::GraphData, i::Int) = gd.e_jac_array[i, 1] ## Micha: gd.e_jac_array[i][1]

@inline function get_dst_edge_jacobian(gd::GraphData, i::Int) = gd.e_jac_array[i, 2]

@inline function get_vertex_jacobian(gd::GraphData, i::Int) = gd.v_jac_array[i]


function update_coefficients!(Jac::NDJacVecOperator, x, p, t)

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

function jac_vec_prod!(dx,Jac::NDJacVecOperator, z)

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


#nd_jac(u) = JacVecProductOperator am Punkt u

#nd_jac(u)(z) = JacVecProduct Jac(u) * z
