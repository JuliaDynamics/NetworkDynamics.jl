#=
# Cascading Failure
This script was written by Hans Würfel. It reimplements the minimal
example of a dynamic cascading failure described in Schäfer et al. (2018) [1].
It is an example how to use callback functions on single nodes or edges.

[1] Schäfer, B., Witthaut, D., Timme, M., & Latora, V. (2018).
Dynamically induced cascading failures in power grids.
Nature communications, 9(1), 1-13.
https://www.nature.com/articles/s41467-018-04287-5

The system is modeled using swing equation and active power edges. The nodes are
characterized by the voltage angle `δ`, the active power on each line is symmetric
and a function of the difference between source and destination angle `δ_src - δ_dst`.
=#

using NetworkDynamics
using Graphs
using OrdinaryDiffEq
using DiffEqCallbacks
using Plots

#=
For the nodes we define the swing equation. State `v[1] = δ`, `v[2] = ω`.
The swing equation has three parameters: `p = (P_ref, I, γ)` where `P_ref`
is the power setpopint, `I` is the inertia and `γ` is the droop or damping coeficcient.
=#
function swing_equation(dv, v, edges, p,t)
    P, I, γ = p
    dv[1] = v[2]
    dv[2] = P - γ * v[2] + flow_sum(edges)
    dv[2] = dv[2] / I
    nothing
end

#=
As an auxilliary function we need to define the flowsum, which adds up all incomming active
power flows to the nodes.
=#
function flow_sum(edges)
    sum = 0.0
    for e in edges
        sum -= e[1]
    end
    return sum
end

#=
Lets define a simple purely active power line whose active power flow is
completlye determined by the connected voltage angles and the coupling constant
`K`.
=#
function simple_edge(e, v_s, v_d, K, t)
    e[1] = K * sin(v_d[1] - v_s[1])
end

#=
The following function returns a `VectorContinoursCallback` which compares the `max_flow` of all
edges in the network an shuts them down once the flow is to high.
=#
function watch_line_limit_callback(max_flow)
    num_e = length(max_flow)

    condition = function(out, u, t, integrator)
        ## get current edge values from integrator
        gd = integrator.f.f(u, integrator.p, integrator.t, GetGD)
        ## collect edge values (the power on each line)
        edge_values = [get_edge(gd, i)[1] for i in 1:ne(g)]
        for i in 1:num_e
            out[i] = max_flow[i] - abs(edge_values[i])
        end
    end

    return VectorContinuousCallback(condition, trip_line_affect!, num_e)
end

#=
This affect will set the coupling strength `K` of a given line to zero. This is equivalent to
removing the edge from the network.
=#
function trip_line_affect!(integrator, idx)
    integrator.p[2][idx] = 0
    auto_dt_reset!(integrator)
end

#=
Axilliary function to plot the line loads at the end.
=#
function plot_flow(nd, saved_edge_data)
    t = saved_edge_data.t
    vals = saved_edge_data.saveval

    data = zeros(length(t), length(vals[1]))
    for i in 1:length(vals)
        # devide by 1.63 to normalize
        data[i, :] = abs.(vals[i]) / 1.63
    end

    g = nd.f.graph
    sym = Array{String}(undef, (1, ne(g)))
    for (i, e) in enumerate(edges(g))
        sym[1, i] = string(e.src) * "->" * string(e.dst)
    end

    plot(t, data; label=sym, size=(800, 600))
end

#=
We can now define the graph topology as in the paper of Schäfter.
=#
g = SimpleGraph([0 1 1 0 1;
                 1 0 1 1 0;
                 1 1 0 1 0;
                 0 1 1 0 1;
                 1 0 0 1 0])

# Definition of nodes and edges according to schäfer18
I = 1.0
γ = 0.1

node_p = [(-1.0, I, γ),
          ( 1.5, I, γ),
          (-1.0, I, γ),
          (-1.0, I, γ),
          ( 1.5, I, γ)]

edge_p = [1.63 for i in 1:ne(g)]

p = (node_p, edge_p)

limits = [0.6*1.63 for i in 1:ne(g)]

# Define nodes/edges and network
odevertex  = ODEVertex(; f=swing_equation, dim=2, sym=[:δ, :ω])
staticedge = StaticEdge(; f=simple_edge, dim=1, sym=[:P], coupling=:antisymmetric)
swing_network = network_dynamics(odevertex, staticedge, g)

# u0 determined from static solution
u0 = [-0.12692637482862684, -1.3649456633810975e-6, 0.14641121510104085, 4.2191082676726005e-7, -0.24376507587890778,
      1.567589744768255e-6, -0.12692637482862684, -1.3649456633810975e-6, 0.35120661043511864, 7.403907552948938e-7]

# define callback to save edge values (for plotting)
function save_edges(u, t, integrator)
    gd = integrator.f.f(u, integrator.p, integrator.t, GetGD)
    return [get_edge(gd, i)[1] for i in 1:ne(g)]
end

saved_edgevalues = SavedValues(Float64, Vector{Float64})
save_callback = SavingCallback(save_edges, saved_edgevalues)

# define problem
prob = ODEProblem(swing_network, u0, (0.0, 6), p)

# define trip_first callback, at t=1.0s remove edge(2,4)
first_trip_idx = findfirst(e->e.src==2 && e.dst==4, collect(edges(g)))
trip_first = PresetTimeCallback(1.0, integrator -> trip_line_affect!(integrator, first_trip_idx))

# define load_watch_callback which will trip lines which exceed the limit
load_watch_callback = watch_line_limit_callback(limits)

# solve the system
sol = solve(prob, Tsit5(); callback=CallbackSet(save_callback, trip_first, load_watch_callback), dtmax=0.001);

# plot solution
plot_flow(swing_network, saved_edgevalues)
