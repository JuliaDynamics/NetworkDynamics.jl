include("src/NetworkDynamics.jl")
include("src/StaticLines.jl")
using LightGraphs
using LinearAlgebra



A = barabasi_albert(10,5)
dnd = NetworkDynamics.diffusive_network_dynamics(A, (dx, x, p, t) -> dx = -x^2)

x0 = ones(10)+rand(10)
dx0 = ones(10)

#dnd(dx0, x0, nothing, 0.)

using DifferentialEquations

dnd_prob = ODEProblem(dnd, x0, (0., 2.))
sol = solve(dnd_prob)
using Plots
plot(sol, legend=false)
println(dnd.L)

#My first edit
#My second edit

g = barabasi_albert(10,3)

line! = (l,x_s,x_t,p,t) -> l .= x_s .- x_t
lines! = [line! for e in edges(g)]

# l_t and l_s is an array of arrays, we have to sum twice to sum over all elements.

function ssum(a)
    if a == []
        0
    else
        sum(a)
    end
end
node! = (dx,x,l_s,l_t,p,t) -> dx = -x + ssum(l_s) - ssum(l_t)
nodes! = [node! for n in vertices(g)]

sl_nd = NetworkDynamics.static_line_network_dynamics(nodes!, lines!, g)

x0 = rand(10)
dx0 = ones(10)

sl_nd(dx0, x0, nothing, 0.)

sl_nd_prob = ODEProblem(sl_nd, x0, (0., 2.))

sl_nd_sol=solve(sl_nd_prob)

plot(sl_nd_sol, legend=false)


#My third edit:
onl = NetworkDynamics.network_dynamics_on_the_line(A, (dx,x,p,t) -> dx = -x^2, )''

onl(dx0,x0,nothing,0.)

onl_prob = ODEProblem(onl, x0, (0., 2.))

sol2=solve(onl_prob)

plot(sol2, legend=false)

#My fourth edit:

nd= NetworkDynamics.network_dynamics(A, (dx,x,p,t) -> dx= -x^2, (dl,l,p,t) -> dl= -l)

v0=ones(35)+rand(35)


nd_prob= ODEProblem(nd,v0, (0.,2.))

sol3=solve(nd_prob)

plot(sol3,legend=false)

a=ones(10)

a[1:5]
