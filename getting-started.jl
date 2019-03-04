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

dnd = NetworkDynamics.diffusive_network_dynamics(g, (dx, x, p, t) -> dx = 0)
ssl = StaticLines.scalar_static_line(nodes,lines,g)
sdl = DynamicLines.scalar_dynamic_line(nodes,dlines,g)

x0 = ones(10)+rand(10)
dx0 = ones(10)
#Because we need initial values for the 10 nodes and! the 25 lines
v0 = rand(35)

dnd_prob = ODEProblem(dnd, x0, (0., 2.))
sol = solve(dnd_prob)
plot(sol, legend=false)

ssl_prob= ODEProblem(ssl, x0 , (0., 2.))
sol2 = solve(ssl_prob)
plot(sol2, legend=false)

sdl_prob = ODEProblem(sdl, v0, (0.,50.))
sol3= solve(sdl_prob)
#to separate lines and nodes, add something like vars=(1:10) in the plot function
plot(sol3, legend=false,vars=(11:35))

line! = (l,x_s,x_t,p,t) -> l .= x_s .- x_t
lines! = [line! for e in edges(g)]
dline! = (dl,l,x_s,x_t,p,t) -> dl .= x_s .- x_t .- l
dlines! =[dline! for e in edges(g)]
node! = (dx,x,l_s,l_t,p,t) -> dx .= x .- sum(l_s) .+ sum(l_t)
nodes! = [node! for n in vertices(g)]

a= NetworkDynamics.static_line_network_dynamics(nodes!,lines!,g)
b= NetworkDynamics.dynamic_line_network_dynamics(nodes!,dlines!,g)

x0=rand(10)
dx0=rand(10)
v0=rand(35)

a(dx0,x0,nothing,0.)
b(v0,v0,nothing,0.)

a_prob=ODEProblem(a,x0,(0.,2.))
sol=solve(a_prob)

plot(sol,legend=false)

b_prob=ODEProblem(b,v0,(0.,5.))
sol2=solve(b_prob)
plot(sol2,legend=false,vars=(11:35))


mline! = [(l,x_s,x_t,p,t) -> l .= x_s .- x_t,(l,x_s,x_t,p,t) -> l .= x_s .- x_t]
mlines! = [mline! for e in edges(g)]
mnode! = [(dx,x,l_s,l_t,p,t) -> dx .= .- sum(l_s) .+ sum(l_t),(dx,x,l_s,l_t,p,t) -> dx .= .- sum(l_s) .+ sum(l_t)]
mnodes! = [mnode! for n in vertices(g)]
mdline! = [(dl,l,x_s,x_t,p,t) -> dl .= .- l[1] .+ x_s .- x_t,(dl,l,x_s,x_t,p,t) -> dl .= .-l[2] .+ x_s .- x_t]
mdlines! = [mdline! for e in edges(g)]

c = NetworkDynamics.multi_static(mnodes!,mlines!,g)

x0=rand(20)

c(x0,x0,0,0)

c_prob=ODEProblem(c,x0,(0.,5.))
sol3=solve(c_prob)
plot(sol3,legend=false)

d = NetworkDynamics.multi_dynamic(mnodes!,mdlines!,g)

v0=rand(70)

d(v0,v0,0,0)

d_prob = ODEProblem(d,v0,(0.,50.))
sol4 = solve(d_prob)
plot(sol4,legend=false,vars=(21:70))
