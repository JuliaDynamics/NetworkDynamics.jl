using Pkg
Pkg.activate(@__DIR__)
using Revise

using NetworkDynamics
using LightGraphs
using LinearAlgebra
using SparseArrays
using OrdinaryDiffEq
using Plots


struct SolAnalytic
    L
end
function (sa::SolAnalytic)(x0, p, t)
    exp(- t * sa.L) * x0
end

g = barabasi_albert(10,5)

L = laplacian_matrix(g)

sol_analytic = SolAnalytic(Symmetric(Array(L)))

diff_network_L = ODEFunction((dx, x, p, t) -> dx .= - L * x, analytic=sol_analytic)

x0 = rand(10)

prob_L = ODEProblem(diff_network_L,x0,(0.,5.))
sol_L = solve(prob_L, Tsit5())

println(sol_L.errors)

plot(sol_L)

sol_ana = zeros(length(sol_L.t), length(x0))
for i in 1:length(sol_L.t)
    sol_ana[i, :] .= sol_analytic(x0, nothing, sol_L.t[i])
end
plot(sol_L.t, sol_ana)

#Now for NetworkDynamics

@inline function diffusion_edge!(e,v_s,v_d,p,t)
    e .= v_s .- v_d
    nothing
end

@inline function diffusion_vertex!(dv, v, e_s, e_d, p, t)
    dv .= 0.
    oriented_edge_sum!(dv, e_s, e_d) # Oriented sum of the incoming and outgoing edges
    nothing
end

odevertex = ODEVertex(f! = diffusion_vertex!,dim = 1)
staticedge = StaticEdge(f! = diffusion_edge!, dim = 1)
odeedge = ODEEdge(staticedge) # We promote the static edge to an ODEEdge artifically

diff_network_st = network_dynamics(odevertex, staticedge, g)

x0 = rand(10)
dx_L = similar(x0)
dx_st = similar(x0)

prob_st = ODEProblem(diff_network_st,x0,(0.,5.))
sol_st = solve(prob_st, Tsit5())
plot(sol_st)

diff_network_L(dx_L, x0, nothing, 0.)
diff_network_st(dx_st, x0, nothing, 0.)

isapprox(dx_L, dx_st)

# ODE Edges

diff_network_ode = network_dynamics(odevertex,odeedge,g)

x0_ode = find_valid_ic(diff_network_ode, randn(10 + 25))
dx0_ode = similar(x0_ode)

prob_st2 = ODEProblem(diff_network_st,x0_ode[1:10],(0.,5.))
prob_ode = ODEProblem(diff_network_ode,x0_ode,(0.,5.))

sol_st2 = solve(prob_st2, Tsit5())
sol_ode = solve(prob_ode, Rodas4())

plot(sol_ode)

vertex_syms = syms_containing(diff_network_ode, "v")

plot(sol_ode, vars=vertex_syms)
plot(sol_st2, vars=1:10)

sol = solve(prob_st2, Rodas4())

using BenchmarkTools

@btime solve(prob_st2, Rodas4(autodiff=false))
@btime solve(prob_st2, Rodas4())
