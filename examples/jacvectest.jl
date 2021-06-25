using Pkg
Pkg.activate(@__DIR__)
using OrdinaryDiffEq
using LinearAlgebra
using NetworkDynamics
using LightGraphs
using Plots
using BenchmarkTools


N = 4


g = SimpleGraph(N)

add_edge!(g, 1, 2)
add_edge!(g, 2, 3)
add_edge!(g, 3, 4)
add_edge!(g, 1, 3)


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

@Base.propagate_inbounds function jac_vertex!(J::AbstractMatrix, v, p, t)
    #J = internal_jacobian(J, v, p, t)
    J[1, 1] = 0.0
    J[2, 1] = 0.0
    nothing
end

@Base.propagate_inbounds function jac_edge!(J_s::AbstractMatrix, J_d::AbstractMatrix, v_s, v_d, p, t)
   J_s[1, 1] = 1.0
#   J_s[2, 1] = 0.0

   J_d[1, 1] = -1.0
#   J_d[2, 1] = 0.0
    nothing
end

nd_jac_vertex = ODEVertex(f! = diffusionvertex!, dim = 2, vertex_jacobian! = jac_vertex!)

nd_jac_edge = StaticEdge(f! = diffusionedge!, dim = 1, coupling = :undirected, edge_jacobian! = jac_edge!)

nd_jac = network_dynamics(nd_jac_vertex, nd_jac_edge, g, jac = true)

nd = network_dynamics(nd_jac_vertex, nd_jac_edge, g, jac = false)

x0 = [1.0, 2.0, 3.0, 4.0, 1.0, 2.0, 3.0, 4.0]

ode_prob_jac = ODEProblem(nd_jac, x0, (0.0, 5.0))
ode_prob = ODEProblem(nd, x0, (0.0, 5.0))

sol_jac = solve(ode_prob_jac,
    TRBDF2(linsolve=LinSolveGMRES()));
sol = solve(ode_prob, TRBDF2(linsolve=LinSolveGMRES()));

print(sol.destats)
print(sol_jac.destats)



dx = randn(2N)
z  = randn(2N)
mul!(dx, nd_jac.jac_prototype, z)
@allocated mul!(dx, nd_jac.jac_prototype, z)



update_coefficients!(nd_jac.jac_prototype, x0, nothing, 0.)
@allocated(update_coefficients!(nd_jac.jac_prototype, z, nothing, 0.))




plot_with_jac = plot(sol_jac, color = :black)
plot!(plot_with_jac, sol, color = :red)#, linestyle = :dash)

@btime solve(ode_prob_jac, TRBDF2(linsolve=LinSolveGMRES())); # 5933 allocations, 1.194 ms
@btime solve(ode_prob, TRBDF2(linsolve=LinSolveGMRES())); # 1721 allocations, 574.314 micro sec = 0.574314 ms
