# You will find a step-by-step guide to this example in the docs and the
# corresponding jupyter notebook on our github repository.

using NetworkDynamics
using LightGraphs
using OrdinaryDiffEq
using Plots

### Defining a graph

N = 20 # number of nodes
k = 4  # average degree
g = barabasi_albert(N, k) # a little more exciting than a bare random graph


### Functions for edges and vertices

function diffusionedge!(e, v_s, v_d, p, t)
    # usually e, v_s, v_d are arrays, hence we use the broadcasting operator .
    e .= v_s .- v_d
    nothing
end

function diffusionvertex!(dv, v, e_s, e_d, p, t)
    # usually dv, v, e_s, e_d are arrays, hence we use the broadcasting operator .
    dv .= 0.
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

nd_diffusion_vertex = ODEVertex(f! = diffusionvertex!, dim = 1)

nd_diffusion_edge = StaticEdge(f! = diffusionedge!, dim = 1)

nd = network_dynamics(nd_diffusion_vertex, nd_diffusion_edge, g)

### Simulation

x0 = randn(N)
ode_prob = ODEProblem(nd, x0, (0., 4.))
sol = solve(ode_prob, Tsit5());

### Plotting

plot(sol, vars = syms_containing(nd, "v"))

# accessing edge values via helper function GetGD
gd_nd = nd(sol(1.0), 1.0, p, GetGD) # exposes underlying graph data struct
e_values_1 = gd_nd.e_array

plot(e_values_1, vars = syms_containing(nd, "e"))
print(e_values_1[:])

# accessing edge values using SavingCallback
