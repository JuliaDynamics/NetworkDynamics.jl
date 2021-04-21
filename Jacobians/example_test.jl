using DifferentialEquations
using Reexport
using LinearAlgebra
using NetworkDynamics
using LightGraphs
using DiffEqBase
import DiffEqBase.update_coefficients!
using Plots


#include("/Users/lalalulu/Desktop/PIK/github/NetworkDynamics.jl/Jacobians/jac_structs.jl")
#@reexport using .jac_structs

#include("/Users/lalalulu/Desktop/PIK/github/NetworkDynamics.jl/src/nd_ODE_Static.jl")
#@reexport using .nd_ODE_Static_mod

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
    #dv[2] = 1.
    for e in edges
        dv[1] += e[1]
    end
    nothing
end

function jac_vertex!(J::AbstractMatrix, v, p, t)
    #J = internal_jacobian(J, v, p, t)
    J[1, 1] = 0.0
end

function jac_edge!(J_s::AbstractMatrix, J_d::AbstractMatrix, v_s, v_d, p, t)
   #J_s = source_jacobian(v_s, v_d, p, t)
   #J_d = dest_jacobian(v_s, v_d, p, t)
   J_s[1, 1] = 1.0
   #J_s[1, 2] = 0.0
   J_d[1, 1] = -1.0
   #J_d[1, 2] = 0.0
end

nd_diffusion_vertex = ODEVertex(f! = diffusionvertex!, dim = 1, vertex_jacobian! = jac_vertex!)

nd_diffusion_edge = StaticEdge(f! = diffusionedge!, dim = 1, edge_jacobian! = jac_edge!)

nd = network_dynamics(nd_diffusion_vertex, nd_diffusion_edge, g, jac = true)
nd2 = network_dynamics(nd_diffusion_vertex, nd_diffusion_edge, g, jac = false)

x0 = randn(N) # random initial conditions
ode_prob = ODEProblem(nd, x0, (0., 4.))
sol = solve(ode_prob, Rodas5());

ode_prob2 = ODEProblem(nd2, x0, (0., 4.))
sol2 = solve(ode_prob2, Rodas5());

plot(sol)
plot(sol2)

@time solve(ode_prob, Rodas5())
@time solve(ode_prob2, Rodas5())
