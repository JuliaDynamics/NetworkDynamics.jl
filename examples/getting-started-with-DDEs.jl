# Experimental DDE tutorial

using NetworkDynamics
using LightGraphs
using OrdinaryDiffEq
using DelayDiffEq
using Plots

### Defining a graph

N = 20 # number of nodes
k = 8  # average degree
g = watts_strogatz(N, k, 0.) # a little more exciting than a bare random graph

# the signature of the edge and vertex functions differs from the ODE signature
function diffusionedge!(e, v_s, v_d, h_v_s, h_v_d, p, t)
    # usually e, v_s, v_d are arrays, hence we use the broadcasting operator .
    e .= v_s - v_d
    nothing
end

function diffusionvertex!(dv, v, e_s, e_d, h_v, p, t)
    # usually v, e_s, e_d are arrays, hence we use the broadcasting operator .
    dv .= -h_v
    # edges for which v is the source
    for e in e_s
        dv .-= e
    end
    # edges for which v is the destination
    for e in e_d
        dv .+= e
    end
    nothing
end

### Constructing the network dynamics

# New Vertex type DDEVertex - DDEEdge might be part of a future release
nd_diffusion_vertex = DDEVertex(f! = diffusionvertex!, dim = 1)
nd_diffusion_edge = StaticEdge(f! = diffusionedge!, dim = 1)

nd = network_dynamics(nd_diffusion_vertex, nd_diffusion_edge, g)

### Simulation

x0 = randn(N) # random initial conditions
# history function defaults to all 1. and is in-place to save allocations
h(out, p, t) = (out .= .5)
tspan = (0., 4.)

# i extended the tuple syntax to pass the delay τ as a parameter
# the first argument should be an array (or other object) containing the vertex parameters
# the second argument holds the edge parameters and the third specifies the delay time τ
# p = (vertexparameters, edgeparameters, delaytime)
p = (nothing, nothing, 1.)

dde_prob = DDEProblem(nd, x0, h, tspan, p)
sol = solve(dde_prob, MethodOfSteps(Tsit5()))

### Plotting

plot(sol, vars = syms_containing(nd, "v"), legend=false)


### Bonus: Two independet diffusions with fancy symbols


# We will have two independent diffusions on the network, hence dim = 2
nd_diffusion_vertex_2 = DDEVertex(f! = diffusionvertex!, dim = 2, sym = [:x, :ϕ])
nd_diffusion_edge_2 = StaticEdge(f! = diffusionedge!, dim = 2)
nd_2 = network_dynamics(nd_diffusion_vertex_2, nd_diffusion_edge_2, g)

# at the moment there are issues with higher dimensional arrays in the DDE solve
# these will be patched
# for now we have to use flat arrays that contain the initial conditions in the right order
# x_0_2 = (x₀_1, ϕ₀_1, x₀_2, ϕ₀_2, x₀_3, ϕ₀_3  ...)

x0_2 = Array{Float64,1}(vec([randn(N).-10 randn(N).^2]')) # x ~ N(0,1); ϕ ~ N(0,1)^2

p = (nothing, nothing, .5) # p = (vertexparameters, edgeparameters, delaytime)
dde_prob_2 = DDEProblem(nd_2, x0_2, h, tspan, p)

sol_2 = solve(dde_prob_2, MethodOfSteps(Tsit5()));


plot(sol_2, legend=false)

using BenchmarkTools

@btime nd($x0_2, $x0_2, $h, 0., 0.)

### Kuramoto model

function kuramoto_delay_edge!(e, v_s, v_d, h_v_s, h_v_d, p, t)
    # usually e, v_s, v_d are arrays, hence we use the broadcasting operator .
    @. e = p * sin(h_v_s - v_d)
    nothing
end

function kuramoto_delay_vertex!(dv, v, e_s, e_d, h_v, p, t)
    # usually v, e_s, e_d are arrays, hence we use the broadcasting operator .
    dv .= p
    # edges for which v is the source
    for e in e_s
        dv .-= e
    end
    # edges for which v is the destination
    for e in e_d
        dv .+= e
    end
    nothing
end

kdedge! = StaticEdge(f! = kuramoto_delay_edge!, dim=1)
kdvertex! = DDEVertex(f! = kuramoto_delay_vertex!, dim = 1)

nd! = network_dynamics(kdvertex!, kdedge!, g)

x0 = randn(N) # random initial conditions
# history function defaults to all 1. and is in-place to save allocations
h(out, p, t) = (out .= pi/2)
# p = (vertexparameters, edgeparameters, delaytime)
p = (randn(N), 10., 1.)
tspan = (0.,10.)
dde_prob = DDEProblem(nd!, x0, h, tspan, p)

sol = solve(dde_prob, MethodOfSteps(Tsit5()));

### Plotting

plot(sol, vars = syms_containing(nd, "v"), legend=false)

@btime nd!($x0, $x0, $h, p, 0.)
