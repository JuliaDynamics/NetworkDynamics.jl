using Pkg
Pkg.activate(@__DIR__)
using Revise

using NetworkDynamics
using LightGraphs
using LinearAlgebra
using SparseArrays
using OrdinaryDiffEq
using Plots

### Solving the problem without NetworkDynamics

struct SolAnalytic
    L
end

# callable struct of type SolAnaltic
function (sa::SolAnalytic)(x0, p, t)
    exp(- t * sa.L) * x0
end

### Defining a graph

# random graph, nodes N=10, nodedegree k=5
g = barabasi_albert(10,5)
L = laplacian_matrix(g)

sol_analytic = SolAnalytic(Symmetric(Array(L)))

### Solving the DDEProblem

# function f is declared with '->' and is callable struct of ODEFunction
diff_network_L = ODEFunction((dx, x, p, t) -> dx .= - L * x, analytic=sol_analytic)

x0 = rand(10) # random initial conditions of length N=10
prob_L = ODEProblem(diff_network_L,x0,(0.,5.))

# method: Tsit5 - non-stiff solver
sol_L = solve(prob_L, Tsit5())

println(sol_L.errors)
plot(sol_L)

# create multidimensional solution array
sol_ana = zeros(length(sol_L.t), length(x0))

# for each row, calculate callable struct of SolAnalytic sa(x0, p, t)
for i in 1:length(sol_L.t)
    sol_ana[i, :] .= sol_analytic(x0, nothing, sol_L.t[i])
end

### Plotting

plot(sol_L.t, sol_ana)


### Now for NetworkDynamics

### Functions for edges and vertices

@inline function diffusion_edge!(e,v_s,v_d,p,t)
    # usually e, v_s, v_d are arrays, hence we use the broadcasting operator .
    e .= v_s .- v_d
    nothing
end

@inline function diffusion_vertex!(dv, v, e_s, e_d, p, t)
    dv .= 0.
    oriented_edge_sum!(dv, e_s, e_d) # Oriented sum of the incoming and outgoing edges
    nothing
end

### Constructing the network dynamics

# constructors of NetworkDynamics
odevertex = ODEVertex(f! = diffusion_vertex!,dim = 1)
staticedge = StaticEdge(f! = diffusion_edge!, dim = 1)
odeedge = ODEEdge(staticedge) # we promote the static edge to an ODEEdge artifically

# setting up the key constructor network_dynamics with static edges
diff_network_st = network_dynamics(odevertex, staticedge, g)

### Simulation

x0 = rand(10) # random initial conditions length N=10
dx_L = similar(x0)
dx_st = similar(x0)

# we use the same methods as above (in the case without NetworkDynamics) to solve and plot
prob_st = ODEProblem(diff_network_st,x0,(0.,5.))
sol_st = solve(prob_st, Tsit5())
plot(sol_st)

### Comparision of the solutions

# callable structs of ODEFunction
# without usage of NetworkDynamics
diff_network_L(dx_L, x0, nothing, 0.)
# with usage of NetworkDynamics
diff_network_st(dx_st, x0, nothing, 0.)

isapprox(dx_L, dx_st)

### Case with ODEEdges

# network dynamics with a StaticEdge that got promoted to an ODEEdge
diff_network_ode = network_dynamics(odevertex,odeedge,g)

### Simulation
# ODEEdge's (with internal dynamics) need initial conditions for the edges
# we first test valid initial conditions for the 10 nodes + 25 edges
x0_ode = find_valid_ic(diff_network_ode, randn(10 + 25))
dx0_ode = similar(x0_ode)

# compare the treatment of the same problem, with StaticEdges and ODEEdges
prob_st2 = ODEProblem(diff_network_st,x0_ode[1:10],(0.,5.))
prob_ode = ODEProblem(diff_network_ode,x0_ode,(0.,5.))

# Rodas4 - stiff solver capable of handling mass matrices
sol_st2 = solve(prob_st2, Tsit5())
sol_ode = solve(prob_ode, Rodas4())

### Plotting

plot(sol_ode)

# extract the indices of the network, containing the 'v' (vertex) -symbol
vertex_syms = syms_containing(diff_network_ode, "v")

plot(sol_st2, vars=1:10)
plot(sol_ode, vars=vertex_syms) # only plotting the vertices, not the edges

# also solving the static problem with the solver Rodas4 for better comparision
sol = solve(prob_st2, Rodas4())

### Benchmarking

using BenchmarkTools

@btime solve(prob_st2, Rodas4(autodiff=false))
@btime solve(prob_st2, Rodas4())
