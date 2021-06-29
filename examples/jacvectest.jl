using Pkg
Pkg.activate(@__DIR__)
using OrdinaryDiffEq
using LinearAlgebra
using NetworkDynamics
using LightGraphs
using Plots
using BenchmarkTools


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

@Base.propagate_inbounds function jac_edge!(J_s::AbstractMatrix, J_d::AbstractMatrix, v_s, v_d, p, t)
   J_s[1, 1] = 1.0
   #J_s[1, 2] = 0.

   J_d[1, 1] = -1.0
   #J_d[1, 2] = 0.
    nothing
end

nd_jac_vertex = ODEVertex(f! = diffusionvertex!, dim = 2, vertex_jacobian! = jac_vertex!)
nd_jac_edge = StaticEdge(f! = diffusionedge!, dim = 1, coupling = :undirected, edge_jacobian! = jac_edge!)

nd_jac = network_dynamics(nd_jac_vertex, nd_jac_edge, g, jac = true)

nd = network_dynamics(nd_jac_vertex, nd_jac_edge, g, jac = false)

nd_jac.jac_prototype.jac_graph_data

x0 = randn(2N^2)

ode_prob_jac = ODEProblem(nd_jac, x0, (0.0, 5.0))
ode_prob = ODEProblem(nd, x0, (0.0, 5.0))

sol_jac = solve(ode_prob_jac, TRBDF2(linsolve=LinSolveGMRES()));
sol = solve(ode_prob, TRBDF2(linsolve=LinSolveGMRES()));

print(sol.destats)
print(sol_jac.destats)


M=N^2
dx = ones(2M)
J = zeros(2M,2M)
z  = randn(2M)

using ForwardDiff

z = [1; zeros(2M-1)]

ForwardDiff.jacobian((out, x) -> nd_jac(out, x, 0., 0.), dx, x0) #*  shuffle!(z)



update_coefficients!(nd_jac.jac_prototype, x0, 0., 0.)
begin
    dx = zeros(2M, 2M)
    for i = 1:2M
        z = zeros(2M)
        z[i] = 1
        mul!(view(dx,:,i), nd_jac.jac_prototype, z)
    end
    dx
end


update_coefficients!(nd_jac.jac_prototype, x0, nothing, 0.)
@allocated(update_coefficients!(nd_jac.jac_prototype, z, nothing, 0.))

dx = ones(2M)
@allocated mul!(dx, nd_jac.jac_prototype, z)

@btime solve(ode_prob_jac, KenCarp4(linsolve=LinSolveGMRES()));
@btime solve(ode_prob, KenCarp4(linsolve=LinSolveGMRES()));


plot_with_jac = plot(sol_jac, color = :black, vars=1:8)
plot!(plot_with_jac, sol, color = :red, linestyle = :dash, vars=1:8)


@btime solve(ode_prob, Tsit5());

## Large

M=501
g = watts_strogatz(M, 2,0)

nd_jac = network_dynamics(nd_jac_vertex, nd_jac_edge, g, jac = true)

nd = network_dynamics(nd_jac_vertex, nd_jac_edge, g, jac = false)

x0 = rand(2M)

ode_prob_jac = ODEProblem(nd_jac, x0, (0.0, 5.0))
ode_prob = ODEProblem(nd, x0, (0.0, 5.0))

@btime solve(ode_prob_jac,  TRBDF2(linsolve=LinSolveGMRES()));
@btime solve(ode_prob, TRBDF2());
@btime solve(ode_prob, Tsit5());

sol_jac = solve(ode_prob_jac,
     TRBDF2(linsolve=LinSolveGMRES()));
sol = solve(ode_prob,  TRBDF2(linsolve=LinSolveGMRES()));

plot_with_jac = plot(sol_jac, vars=[1,2,3,4], color = :black)
plot!(plot_with_jac, sol,  vars=[1,2,3,4], color = :red)#, linestyle = :dash)


print(sol.destats)
print(sol_jac.destats)

### HUGE

M= 500000
g = watts_strogatz(M, 2,0)

nd_jac = network_dynamics(nd_jac_vertex, nd_jac_edge, g, jac = true)

nd = network_dynamics(nd_jac_vertex, nd_jac_edge, g, jac = false)

x0 = rand(2M)

ode_prob_jac = ODEProblem(nd_jac, x0, (0.0, 5.0))
ode_prob = ODEProblem(nd, x0, (0.0, 5.0))

@btime solve(ode_prob_jac,  TRBDF2(linsolve=LinSolveGMRES()));
@time solve(ode_prob, TRBDF2());
@btime solve(ode_prob,Tsit5());

sol_jac = solve(ode_prob_jac,
     TRBDF2(linsolve=LinSolveGMRES()));
sol = solve(ode_prob,  KenCarp4(linsolve=LinSolveGMRES()));

print(sol.destats)
print(sol_jac.destats)
plot_with_jac = plot(sol_jac, vars=[10, 20, 30, 40],color = :black)
