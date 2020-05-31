# You will find a step-by-step guide to this example in the docs and the
# corresponding jupyter notebook on our github repository.

using NetworkDynamics
using LightGraphs
using OrdinaryDiffEq
using DelayDiffEq
using Plots

### Defining a graph

N = 20 # number of nodes
k = 4  # average degree
g = barabasi_albert(N, k) # a little more exciting than a bare random graph


function diffusionedge!(e, v_s, v_d, p, t)
    # usually e, v_s, v_d are arrays, hence we use the broadcasting operator .
    e .= v_s - v_d
    nothing
end

function diffusionvertex!(dv, v, e_s, e_d, h, p, t)
    # usually v, e_s, e_d are arrays, hence we use the broadcasting operator .
    dv .= h
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

nd_diffusion_vertex = DDEVertex(f! = diffusionvertex!, dim = 1)
nd_diffusion_edge = StaticEdge(f! = diffusionedge!, dim = 1)

nd = network_dynamics(nd_diffusion_vertex, nd_diffusion_edge, g)

### Simulation

x0 = randn(N) # random initial conditions
h(out, p, t) = (out .= 1.) #constant history

ode_prob = DDEProblem(nd, x0, h, (0., 5.))
sol = solve(ode_prob, MethodOfSteps(Tsit5()));

### Plotting

plot(sol, vars = syms_containing(nd, "v"))


### Bonus: Two independet diffusions with fancy symbols


# We will have two indepent diffusions on the network, hence dim = 2
nd_diffusion_vertex_2 = DDEVertex(f! = diffusionvertex!, dim = 2, sym = [:x, :ϕ])
nd_diffusion_edge_2 = StaticEdge(f! = diffusionedge!, dim = 2)
nd_2 = network_dynamics(nd_diffusion_vertex_2, nd_diffusion_edge_2, g)

x0_2 = Array{Float64,1}(vec([10*randn(N) randn(N).^2]')) # x ~ N(0,1); ϕ ~ x^2
ode_prob_2 = DDEProblem(nd_2, x0_2, h,(0., 4.))
sol_2 = solve(ode_prob_2, MethodOfSteps(Tsit5()));


# Try plotting the variables ϕ_i yourself. [Type \phi and press TAB]
plot(sol_2)
