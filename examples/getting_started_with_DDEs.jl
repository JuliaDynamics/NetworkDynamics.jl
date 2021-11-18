# Corresponds to the DDE tutorial in the docs

using NetworkDynamics
using Graphs
using OrdinaryDiffEq
using DelayDiffEq
using Plots

### Defining a graph

N = 20 # number of nodes
k = 8  # average degree
g = watts_strogatz(N, k, 0.0) # ring network

# the signature of the edge and vertex functions differs from the ODE signature
function diffusionedge!(e, v_s, v_d, p, t)
    e .= 0.1 * (v_s - v_d)
    nothing
end

function diffusionvertex!(dv, v, edges, h_v, p, t)
    dv .= -h_v
    sum_coupling!(dv, edges)
    nothing
end

### Constructing the network dynamics

# DDEVertex and StaticDelayEdge  both have access to the VERTEX history
# DDEVertex is expected to have call signature (dv, v, edges, h_v, p, t)
# StaticDelayEdge is expected to have  signature (e, v_s, v_d, h_v_s, h_v_d, p, t)
# StaticEdges get promoted to StaticDelayEdges [then their signature changes]
nd_diffusion_vertex = DDEVertex(; f=diffusionvertex!, dim=1)
nd_diffusion_edge = StaticEdge(; f=diffusionedge!, dim=1)

nd = network_dynamics(nd_diffusion_vertex, nd_diffusion_edge, g)

### Simulation

x0_1 = randn(N) # random initial conditions
# constant history function, is in-place to save allocations
h(out, p, t) = (out .= x0_1)
tspan = (0.0, 10.0)

# fo DDES the tuple parameter syntax is extended to hold the delay time τ as a parameter
# the first argument should be an array (or other object) containing the vertex parameters
# the second argument holds the edge parameters and the third specifies the delay time τ
# p = (vertexparameters, edgeparameters, delaytime)
p = (nothing, nothing, 1.0)

# you might want to use the keyword constant_lags here. check the docs of DelayDiffEq for details
dde_prob = DDEProblem(nd, x0_1, h, tspan, p)
sol = solve(dde_prob, MethodOfSteps(Tsit5()))

### Plotting

plot(sol; vars=syms_containing(nd, "v"), legend=false)


### Bonus: Two independent diffusions with fancy symbols


# We will have two independent diffusions on the network, hence dim = 2
nd_diffusion_vertex_2 = DDEVertex(; f=diffusionvertex!, dim=2, sym=[:x, :ϕ])
nd_diffusion_edge_2 = StaticEdge(; f=diffusionedge!, dim=2)
nd_2 = network_dynamics(nd_diffusion_vertex_2, nd_diffusion_edge_2, g)

# at the moment there are issues with higher dimensional arrays in the DDE solve
# these will be patched
# for now we have to use flat arrays that contain the initial conditions in the right order
# x_0_2 = (x₀_1, ϕ₀_1, x₀_2, ϕ₀_2, x₀_3, ϕ₀_3  ...)

const x0_2 = Vector{Float64}(vec([randn(N) .- 10 randn(N) .^ 2]')) # x ~ N(0,1); ϕ ~ N(0,1)^2
h(out, p, t) = (out .= x0_2)
p = (nothing, nothing, 1.0) # p = (vertexparameters, edgeparameters, delaytime)
dde_prob_2 = DDEProblem(nd_2, x0_2, h, tspan, p)



sol_2 = solve(dde_prob_2, MethodOfSteps(Tsit5()));


plot(sol_2; legend=false)

### Kuramoto model

function kuramoto_delay_edge!(e, v_s, v_d, h_v_s, h_v_d, p, t)
    e[1] = p * sin(v_s[1] - h_v_d[1])
    nothing
end

function kuramoto_vertex!(dv, v, edges, p, t)
    dv[1] = p
    for e in edges
        dv .+= e
    end
    nothing
end

kdedge! = StaticDelayEdge(; f=kuramoto_delay_edge!, dim=1)
kdvertex! = ODEVertex(; f=kuramoto_vertex!, dim=1)

nd! = network_dynamics(kdvertex!, kdedge!, g)

const x0_3 = randn(N) # random initial conditions
h(out, p, t) = (out .= x0_3)
# p = (vertexparameters, edgeparameters, delaytime)
ω = randn(N)
ω .-= sum(ω) / N
p = (ω, 2.0, 1.0)
tspan = (0.0, 15.0)
dde_prob = DDEProblem(nd!, x0_3, h, tspan, p)

sol = solve(dde_prob, MethodOfSteps(Tsit5()));

### Plotting

plot(sol; vars=syms_containing(nd, "v"), legend=false)


## Add in a static Vertex


svertex! = StaticVertex(; f=(x, edges, p, t) -> x .= 1, dim=1)

# convert to list of VertexFunctions, that may contain different types of vertices
vlist = Array{VertexFunction}([kdvertex! for i in 1:nv(g)])
vlist[1] = svertex!
# at the moment if either edges or vertices is a list, the other has to be a list as well
elist = [kdedge! for i in 1:ne(g)]

snd! = network_dynamics(vlist, elist, g)

# adjust initial condition to meet the constraint
x0_1[1] = 1.0


ddae_prob = DDEProblem(snd!, x0_1, h, tspan, p)

# autodiff fails for this problem, tuple params don't seem to be the problem, but it might have to do with https://docs.sciml.ai/stable/basics/faq/#I-get-Dual-number-errors-when-I-solve-my-ODE-with-Rosenbrock-or-SDIRK-methods-1

# until we fix this turn off autodiff
sol = solve(ddae_prob, MethodOfSteps(Rosenbrock23(; autodiff=false)));
plot(sol; vars=syms_containing(nd, "v"), legend=false)
