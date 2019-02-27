include("src/NetworkDynamics.jl")
include("src/StaticLines.jl")
include("src/DynamicLines.jl")

using LightGraphs
using LinearAlgebra
using DifferentialEquations
using Plots

g = barabasi_albert(10,5)
line(x_s,x_t,p,t) = x_s - x_t
lines = [line for e in edges(g)]
dline(l,x_s,x_t,p,t) = x_s -x_t - l
dlines = [dline for e in edges(g)]
node(x,l_s,l_t,p,t) = x-sum(l_s) + sum(l_t)
nodes = [node for n in vertices(g)]

dnd = NetworkDynamics.diffusive_network_dynamics(g, (dx, x, p, t) -> dx = -x^2)
ssl = StaticLines.scalar_static_line(nodes,lines,g)
sdl = DynamicLines.scalar_dynamic_line(nodes,dlines,g)

x0 = ones(10)+rand(10)
dx0 = ones(10)
#Because we need initial values for the 10 nodes and! the 25 lines
v0 = rand(35)

dnd_prob = ODEProblem(dnd, x0, (0., 2.))
sol = solve(dnd_prob)
plot(sol, legend=false)
println(dnd.L)

ssl_prob= ODEProblem(ssl, x0 , (0., 2.))
sol2 = solve(ssl_prob)
plot(sol2, legend=false)

sdl_prob = ODEProblem(sdl, v0, (0.,2.))
sol3= solve(sdl_prob)
#to separate lines and nodes, add something like vars=(1:10) in the plot function
plot(sol3, legend=false,vars=(11:35))
