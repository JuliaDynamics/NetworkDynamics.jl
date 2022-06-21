# Introduction to stochastic modeling based on getting-started.jl
# For introductory tutorials to NetworkDynamics visit the docs at
# https://pik-icone.github.io/NetworkDynamics.jl/dev/

using NetworkDynamics
using Graphs
using OrdinaryDiffEq
using StochasticDiffEq
using Plots

### Defining a graph

N = 20 # number of nodes
k = 4  # average degree
g = barabasi_albert(N, k) # a little more exciting than a bare random graph


### Functions for edges and vertices

function diffusionedge!(e, v_s, v_d, p, t)
    # usually e, v_s, v_d are arrays, hence we use the broadcasting operator .
    e .= v_s - v_d
    nothing
end

function diffusionvertex!(dv, v, edges, p, t)
    # usually v, e_s, e_d are arrays, hence we use the broadcasting operator .
    dv .= 0.0
    # edges for which v is the source
    for e in edges
        dv .+= e
    end
    nothing
end


### Constructing the network dynamics

nd_diffusion_vertex = ODEVertex(; f=diffusionvertex!, dim=1)
nd_diffusion_edge = StaticEdge(; f=diffusionedge!, dim=1)

nd = network_dynamics(nd_diffusion_vertex, nd_diffusion_edge, g)

### Simulation & Plotting

x0 = randn(N) # random initial conditions
tspan = (0.0, 4.0)
ode_prob = ODEProblem(nd, x0, tspan)
sol = solve(ode_prob, Tsit5());

plot(sol; vars=syms_containing(nd, "v"))

### The noise can be thought of as an extra layer to the network
# In this case we have independent additive noise at the nodes
# Therefor the extra layer doesn't need any edges
# More complicated noise structures are possible but have to be implemented with care

h = SimpleGraph(N, 0)

function noisevertex!(dv, v, edges, p, t)
    # the noise is written as an ODE, here dv_noise is a square wave
    # when added in the SDE it seems to affect v directly and not dv
    dv[1] = round(Int, cos(t * 0.5 * pi)) * p
end

nd_noisevertex = ODEVertex(; f=noisevertex!, dim=1)
nd_noise = network_dynamics(nd_noisevertex, nd_diffusion_edge, h)


### Simulation of noise & plotting

# we are using the tuple parameter syntax, i.e. p[i] is passed to the corresponing node
p = (randn(N), nothing)
tspan = (0.0, 15.0)
noise_prob = ODEProblem(nd_noise, x0, tspan, p)
sol_noise = solve(noise_prob, Tsit5())

# The plot may be a bit misleading, since we are seeing v and not dv
plot(sol_noise; vars=syms_containing(nd, "v"))


### Setting up the SDEProblem is straightforward

# SDEProblems can be constructed by inserting ODEFunctions for `f` and `g`
# even the symbols stay the same
sde_prob = SDEProblem(nd, nd_noise, x0, tspan, p)
# SOSRA is recommended for additive noise
sde_sol = solve(sde_prob, SOSRA())

plot(sde_sol; vars=syms_containing(nd, "v"), legend=false)
