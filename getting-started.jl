include("src/NetworkDynamics.jl")

using LightGraphs


A = barabasi_albert(40,6)
dnd = NetworkDynamics.diffusive_network_dynamics(A, (dx, x, p, t) -> dx = -x)

x0 = ones(40) + rand(40)
dx0 = zeros(40)

dnd(dx0, x0, nothing, 0.)

using DifferentialEquations

dnd_prob = ODEProblem(dnd, x0, (0., 1.))
sol = solve(dnd_prob)
using Plots
plot(sol, legend=false)
println(dnd.L)

#My first edit
