#=
TODO
Clean up.
=#

include("src/NetworkDynamics.jl")
using .NetworkDynamics
using LightGraphs
using LinearAlgebra
using DifferentialEquations
using Plots

g = barabasi_albert(10,5)

line! = (l,x_s,x_t,p,t) -> l .= p*(x_s - x_t)
lines! = [line! for e in edges(g)]
dline! = (dl,l,x_s,x_t,p,t) -> dl .= x_s - x_t - l
dlines! =[dline! for e in edges(g)]
node! = (dx,x,l_s,l_t,p,t) -> dx .= - sum(l_s) + sum(l_t)
nodes! = [node! for n in vertices(g)]

a= nd_ODE_Static_scalar_mod.nd_ODE_Static_scalar(nodes!,lines!,g)
b= nd_ODE_ODE_scalar(nodes!,dlines!,g)

x0=rand(10)
dx0=rand(10)
z0=rand(35)

a(dx0,x0,ones(35),0.)
b(z0,z0,ones(35),0.)

a_prob=ODEProblem(a,x0,(0.,2.),ones(35)*0.5)
sol=solve(a_prob)

plot(sol,legend=false)
b_prob=ODEProblem(b,z0,(0.,5.))
sol2=solve(b_prob)
plot(sol2,legend=false,vars=(1:10))

function diffusion_vertex!(dv, v, e_ss, e_ds, p, t)
    # Note that e_ss and e_ds might be empty, the code needs to be able to deal
    # with this situation.
    dv .= 0
    for e_s in e_ss
        dv .-= e_s
    end
    for e_d in e_ds
        dv .+= e_d
    end
    nothing
end

vertices! = [diffusion_vertex! for vertex in vertices(g)]
edges! = [(l,x_s,x_t,p,t) -> l .= x_s - x_t for edge in edges(g)]
dedges! = [(dl,l,x_s,x_t,p,t) -> dl .= x_s - x_t - l for edge in edges(g)]

dim_v = 2 * ones(Int32, nv(g))
dim_e = 2* ones(Int32, ne(g))

c = nd_ODE_Static(vertices!,edges!,g,dim_v,dim_e)

x0=ones(20) + rand(20)

c(x0,x0,nothing,0)

c_prob=ODEProblem(c,x0,(0.,5.))
sol3=solve(c_prob)
plot(sol3,legend=false)

d = nd_ODE_ODE(vertices!,dedges!,g,dim_v,dim_e)

v0=rand(70)

d(v0,v0,nothing,0)

d_prob = ODEProblem(d,v0,(0.,50.))
sol4 = solve(d_prob)
plot(sol4,legend=false,vars=(21:70))
