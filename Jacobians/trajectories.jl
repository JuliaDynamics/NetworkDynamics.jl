using DifferentialEquations
using Reexport
using LinearAlgebra
using NetworkDynamics
using LightGraphs
using DiffEqBase
import DiffEqBase.update_coefficients!
using Plots
#using GraphPlot


N = 4
#k = 2
#g = barabasi_albert(N, k)

# building simple ring network consisting of 4 nodes
g = SimpleGraph(N)
add_edge!(g, 1, 2)
add_edge!(g, 2, 3)
add_edge!(g, 3, 4)
add_edge!(g, 4, 1)

function diffusionedge!(e, v_s, v_d, p, t)
    # usually e, v_s, v_d are arrays, hence we use the broadcasting operator .
    #e .= v_s .- v_d
    e[1] = v_s[1] - v_d[1]
    nothing
end

function diffusionvertex!(dv, v, edges, p, t)
    # usually dv, v, edges are arrays, hence we use the broadcasting operator .
    dv[1] = 0.
    for e in edges
        dv[1] += e[1]
    end
    nothing
end

function jac_vertex!(J::AbstractMatrix, v, p, t)
    J[1, 1] = 0.0
end

function jac_edge!(J_s::AbstractMatrix, J_d::AbstractMatrix, v_s, v_d, p, t)
   J_s[1, 1] = 1.0
   J_s[2, 1] = 0.0

   J_d[1, 1] = -1.0
   J_d[2, 1] = 0.0
end

nd_diffusion_vertex = ODEVertex(f! = diffusionvertex!, dim = 1, vertex_jacobian! = jac_vertex!)

nd_diffusion_edge = StaticEdge(f! = diffusionedge!, dim = 1, edge_jacobian! = jac_edge!)

nd_jac = network_dynamics(nd_diffusion_vertex, nd_diffusion_edge, g, jac = true)
nd = network_dynamics(nd_diffusion_vertex, nd_diffusion_edge, g, jac = false)

#x0 = randn(N) # random initial conditions
x0 = [1.0, 2.0, 3.0, 4.0]

ode_prob_jac = ODEProblem(nd_jac, x0, (0.0, 5.0))
ode_prob = ODEProblem(nd, x0, (0.0, 5.0))

sol_jac = solve(ode_prob_jac, Rodas5());
sol = solve(ode_prob, Rodas5());

plot_with_jac = plot(sol_jac, color = [:black])
plot!(plot_with_jac, sol, color = [:red])


### testing different solvers

sol_jac = solve(ode_prob_jac)
sol = solve(ode_prob)

plot_with_jac = plot(sol_jac, color = [:black], tspan=(1e-2,1e5), xscale=:log10)
plot!(plot_with_jac, sol, color = [:red], tspan=(1e-2,1e5), xscale=:log10)

# more efficient
sol_jac = solve(ode_prob_jac, Rodas5());
sol = solve(ode_prob, Rodas5());

plot_with_jac = plot(sol_jac, color = [:black], tspan=(1e-2,1e5), xscale=:log10)
plot!(plot_with_jac, sol, color = [:red], tspan=(1e-2,1e5), xscale=:log10)
# more reliable
sol = solve(ode_prob, Rodas4P());
sol_jac = solve(ode_prob_jac, Rodas4P());

plot_with_jac = plot(sol_jac, color = [:black], tspan=(1e-2,1e5), xscale=:log10)
plot!(plot_with_jac, sol, color = [:red], tspan=(1e-2,1e5), xscale=:log10)


@time solve(ode_prob_jac, Rodas5())
@time solve(ode_prob, Rodas5())

@allocated solve(ode_prob_jac, Rodas5())
@allocated solve(ode_prob, Rodas5())
