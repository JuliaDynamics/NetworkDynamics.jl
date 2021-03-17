# using NetworkDynamics

using LightGraphs
using Reexport

include("/Users/lalalulu/Desktop/PIK/github/NetworkDynamics.jl/src/Utilities.jl")
@reexport using .Utilities

include("/Users/lalalulu/Desktop/PIK/github/NetworkDynamics.jl/Jacobians/jac_structs.jl")
@reexport using .jac_structs

include("/Users/lalalulu/Desktop/PIK/github/NetworkDynamics.jl/Jacobians/jac_functions.jl")
@reexport using .jac_functions


function diffusionedge!(e, v_s, v_d, p, t)
    # usually e, v_s, v_d are arrays, hence we use the broadcasting operator .
    e .= v_s .- v_d
    nothing
end

function diffusionvertex!(dv, v, edges, p, t)
    # usually dv, v, edges are arrays, hence we use the broadcasting operator .
    dv .= 0.
    for e in edges
        dv .+= e
    end
    nothing
end

# Jacobian stuff 1

function VertexJacobian!(J::AbstractMatrix, v, p, t)
    J[1] = 0
end

function EdgeJacobian!(J_s::AbstractMatrix, J_d::AbstractMatrix, v_s, v_d, p, t)
   J_s[1] = 0
   J_s[2] = 0

   J_d[1] = 0
   J_d[2] = 0
end

nd_diffusion_vertex = ODEVertex(f! = diffusionvertex!, dim = 1)

nd_diffusion_edge = StaticEdge(f! = diffusionedge!, dim = 1)

nd = network_dynamics(nd_diffusion_vertex, nd_diffusion_edge, g, jac = true)

### jacobian stuff 2

function update_coefficients!(Jac::NDJacVecOperator1, x, p, t)
#function update_coefficients!(Jac::NDJacVecOperator1{uType, P, tType, G, GD, JGD, T}, x, p, t)
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





### Simulation

x0 = randn(N) # random initial conditions
ode_prob = ODEProblem(nd, x0, (0., 4.))
sol = solve(ode_prob, Tsit5());
