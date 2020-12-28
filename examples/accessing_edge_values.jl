using NetworkDynamics
using LightGraphs
using OrdinaryDiffEq
using Plots

### Defining a graph

N = 20 # number of nodes
k = 4  # average degree
g = barabasi_albert(N, k)


### Functions for edges and vertices

function diffusionedge!(e, v_s, v_d, p, t)
    e .= v_s .- v_d
    nothing
end

function diffusionvertex!(dv, v, edges, p, t)
    dv .= 0.
    for e in edges
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
gd_nd = nd(sol(1.0), 1.0, nothing, GetGD) # exposes underlying graph data struct
e_values = gd_nd.gdb.e_array


# accessing edge values using SavingCallback

using DiffEqCallbacks: SavingCallback, SavedValues

saved_values = SavedValues(Float64, Vector{Float64})
function saving_func(u, t, integrator)
    edgevals = Float64[]
    for i in 1:integrator.f.f.graph_structure.num_e
        push!(edgevals, integrator.f.f.graph_data.e[i]...)
    end
    edgevals
end

cb = SavingCallback(saving_func, saved_values)
sol = solve(ode_prob,Tsit5(),callback=cb)

saved_values # stores e_values
