using Pkg
Pkg.activate(@__DIR__)
using OrdinaryDiffEq
using LinearAlgebra
using NetworkDynamics
using LightGraphs
using Plots
using BenchmarkTools
using Test


N = 4

g = watts_strogatz(N^2,N,0.)


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

@Base.propagate_inbounds function jac_vertex!(J_v::AbstractMatrix, J_e::AbstractMatrix, v, p, t)
    # J_v = internal_jacobian(J, v, p, t)
    J_v[1, 1] = 0.0
    J_v[2, 1] = 0.0

    J_e[1, 1] = 1.
    J_e[2, 1] = 1.
    nothing
end


nd_jac_vertex = ODEVertex(f! = diffusionvertex!, dim = 2, vertex_jacobian! = jac_vertex!)
nd_jac_edge = StaticEdge(f! = diffusionedge!, dim = 1, coupling = :undirected, edge_jacobian! = true)

# test jacobians of diffusionedge!

J_s = zeros(1, 2)
J_d = zeros(1, 2)
v_s = [1, 2]
v_d = [1, 2]
p = 0
t = 1.0

nd_jac_edge.edge_jacobian!(J_s, J_d, v_s, v_d, p, t)

@test J_s == [1.0 0.0]
@test J_d == [-1.0 0.0]



nd_jac = network_dynamics(nd_jac_vertex, nd_jac_edge, g, jac = true)


nd = network_dynamics(nd_jac_vertex, nd_jac_edge, g, jac = false)

nd_jac.jac_prototype.jac_graph_data


x0 = randn(2N^2)

ode_prob_jac = ODEProblem(nd_jac, x0, (0.0, 5.0))
ode_prob = ODEProblem(nd, x0, (0.0, 5.0))

sol_jac = solve(ode_prob_jac, TRBDF2(linsolve=LinSolveGMRES()));
sol = solve(ode_prob, TRBDF2(linsolve=LinSolveGMRES()));
