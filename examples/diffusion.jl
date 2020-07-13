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
# takes x0=initial conditions, p=parameters, t=time
function (sa::SolAnalytic)(x0, p, t)
    exp(- t * sa.L) * x0
end

### Defining a graph

# constructing a Barabasi-Albert graph with N=10 nodes and nodedegree k=5
# signature: barabasi_albert(N,k), method of LightGraphs
g = barabasi_albert(10,5)
L = laplacian_matrix(g)

# first some broadcasting operations are made:
# Array(L) converts 2 dimensional- to multidimensional array
# Symmetric (Array(L)) ensures that object is symmetric
# then we assign the symmetric multidemsnional array to the constructor SolAnalytic
sol_analytic = SolAnalytic(Symmetric(Array(L)))

### Solving the DDEProblem

# ODEFunction is struct of OrdinaryDiffEq
# signature: ODEFunction(f, analytic=nothing, mass_matrix=I, jac=...)
# function f is declared with '->' and is callable struct of ODEFunction
diff_network_L = ODEFunction((dx, x, p, t) -> dx .= - L * x, analytic=sol_analytic)

x0 = rand(10) # random initial conditions of length N=10

# ODEProblem is a struct of OrdinaryDiffEq
# Signature: ODEProblem(f,u0,tspan, ...)
# by handing over the object diff_network_L, the callable struct f of the ODEFunction object is delivered to ODEProblem
prob_L = ODEProblem(diff_network_L,x0,(0.,5.))

# solve is a struct of OrdinaryDiffEq
# signature: solve(prob::ODEProblem,alg;kwargs)
# choose solver method with the alg keyword in solve
# method: Tsit5 - non-stiff solver
sol_L = solve(prob_L, Tsit5())

println(sol_L.errors)
plot(sol_L)

# create multidimensional solution array with zeros
# dimension (timesteps, nodes)
sol_ana = zeros(length(sol_L.t), length(x0))

# for each row, calculate callable struct of SolAnalytic sa(x0, p, t)
for i in 1:length(sol_L.t)
    sol_ana[i, :] .= sol_analytic(x0, nothing, sol_L.t[i])
end

### Plotting

plot(sol_L.t, sol_ana)


### Now for NetworkDynamics

### Functions for edges and vertices

# edge and vertex functions of NetworkDynmaics are modifying functions, so we use '!'
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

# constructors of NetworkDynamics with signatures:
# ODEVertex{T}(f!,dim, sym, mass_matrix)
# StaticEdge{T}(f!, dim, sym)
# ODEEdge{T}(f!, dim, sym, mass_matrix); here also the edges can have internal dynamics
odevertex = ODEVertex(f! = diffusion_vertex!,dim = 1)
staticedge = StaticEdge(f! = diffusion_edge!, dim = 1)
odeedge = ODEEdge(staticedge) # we promote the static edge to an ODEEdge artifically

# setting up the key constructor network_dynamics with static edges
diff_network_st = network_dynamics(odevertex, staticedge, g)

### Simulation

x0 = rand(10) # random initial conditions length N=10
# create empty arrays of same size as x0
dx_L = similar(x0)
dx_st = similar(x0)

# we use same methods as above (in the case without NetworkDynamics) to solve and plot
prob_st = ODEProblem(diff_network_st,x0,(0.,5.))
sol_st = solve(prob_st, Tsit5())
plot(sol_st)

### Comparision of the solutions

# without usage of NetworkDynamics: diff_network_L
# with usage of NetworkDynamics: diff_network_st

# diff_network_L is object of type ODEFunction
# when we hand over parameters (dx,x0,p,t) callable struct is called
diff_network_L(dx_L, x0, nothing, 0.)

# diff_network_st is object of type NetworkDynamics
# callable struct calls ODEFunction
diff_network_st(dx_st, x0, nothing, 0.)

isapprox(dx_L, dx_st) # inexact equality comparision

### Case with ODEEdges

# now we create the network dynamics object with the artificially as ODEEdge promoted static edges (line 91)
diff_network_ode = network_dynamics(odevertex,odeedge,g)

### Simulation
# as the edges are declared as ODEEdge's now (with internal dynamics) we also need initial conditions for the edges
# we first test valid initial conditions for the 10 nodes + 25 edges
x0_ode = find_valid_ic(diff_network_ode, randn(10 + 25))
dx0_ode = similar(x0_ode)

# compare the treatment of the same problem, with StaticEdges and ODEEdges
prob_st2 = ODEProblem(diff_network_st,x0_ode[1:10],(0.,5.))
prob_ode = ODEProblem(diff_network_ode,x0_ode,(0.,5.))

# Rodas4 - stiff solver
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
