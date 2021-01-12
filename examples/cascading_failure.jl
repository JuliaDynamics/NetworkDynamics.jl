println("This script works with NetworkDynamics.jl 0.4 but is broken on 0.5 due to syntax changes.")


# This script was written by Hans Würfel. It reimplements the minimal
# example of a dynamic cascading failure described in Schäfer et al. (2018) [1].
# It is an example how to use callback functions on single nodes or edges.
#
# [1] Schäfer, B., Witthaut, D., Timme, M., & Latora, V. (2018).
# Dynamically induced cascading failures in power grids.
# Nature communications, 9(1), 1-13.
# https://www.nature.com/articles/s41467-018-04287-5




using NetworkDynamics
using LightGraphs
using OrdinaryDiffEq
using DiffEqCallbacks
using Plots

struct SwingVertex
    P::Float64 # power of node
    I::Float64 # inertia of node
    γ::Float64 # damping of node
end
function (sv::SwingVertex)(dv, v, edges, p, t)
    # v[1] -> δ, dv[1] = dγ = ω = v[2]
    # v[2] -> ω, dv[2] -> dω
    dv[1] = v[2]
    dv[2] = sv.P - sv.γ*v[2] + flow_sum(edges)
    dv[2] = dv[2]/sv.I
    nothing
end

function flow_sum(edges)
    sum = 0.0
    for e in edges
        sum -= e[1]
    end
    return sum
end

mutable struct SwingEdge
    from::Int8
    to::Int8
    K::Float64
    max_flow::Float64
end

function (se::SwingEdge)(e, v_s, v_d, p, t)
    e[1] = se.K * sin(v_d[1] - v_s[1])
end

"""
Function returns callable object, which can be used as 'affect!' of a
callback. The affect! will set the coupling of edge idx to zero.
"""
function kill_line_affect(idx)
    function affect!(integrator)
        @info "Line $idx killed at t=$(integrator.t)."
        integrator.f.f.edges![idx].f!.K = 0
        u_modified!(integrator, false)
    end
    return affect!
end


"""
Function returns callable object, which can be used as 'affect!' of a
callback. The affect! will set the coupling of edge (source, destination)
in network to zero.
"""
function kill_line_affect(network, source, destination)
    gs = network.f.graph_structure
    s_ind = findall(x->x==source, gs.s_e)
    d_ind = findall(x->x==destination, gs.d_e)
    ind = intersect(s_ind, d_ind)
    @assert length(ind) == 1
    return kill_line_affect(ind[1])
end


"""
Function returns a VectorContinoursCallback which compares the max_flow of all
edges in the network an shuts them down once the flow is to high.
"""
function watch_line_limit_callback(network)
    num_e = network.f.graph_structure.num_e
    edges = network.f.edges!

    function condition(out, u, t, integrator)
        # get current edge values from integrator
        graph_data = integrator.f.f(u, integrator.p, integrator.t, GetGD)
        edge_values = graph_data.gdb.e_array
        for i in 1:num_e
            out[i] = edges[i].f!.max_flow - abs(edge_values[i])
        end
    end

    function affect!(integrator, idx)
        kill_line_affect(idx)(integrator)
    end
    return VectorContinuousCallback(condition, affect!, num_e)
end

function plot_flow(saved_edge_data, network)
    t = saved_edge_data.t
    vals = saved_edge_data.saveval

    data = zeros(length(t), length(vals[1]))
    for i in 1:length(vals)
        # devide by 1.63 to normalize
        data[i, :] = abs.(vals[i])/1.63
    end

    gs = network.f.graph_structure
    sym = Array{String}(undef, (1,gs.num_e))
    for i in 1:gs.num_e
        sym[1,i] = string(gs.s_e[i]) * "->" * string(gs.d_e[i])
    end

    plot(t, data, label=sym, size=(800,600))
end

# Definition of nodes and edges according to schäfer18
I = 1.0
γ = 0.1
verticies = [
    SwingVertex(-1.0, I, γ), #1
    SwingVertex( 1.5, I, γ), #2
    SwingVertex(-1.0, I, γ), #3
    SwingVertex(-1.0, I, γ), #4
    SwingVertex( 1.5, I, γ), #5
]

swingedges = [
    SwingEdge(1, 2, 1.63, 0.6*1.63), #1
    SwingEdge(1, 3, 1.63, 0.6*1.63), #2
    SwingEdge(1, 5, 1.63, 0.6*1.63), #3
    SwingEdge(2, 3, 1.63, 0.6*1.63), #4
    SwingEdge(2, 4, 1.63, 0.6*1.63), #5
    SwingEdge(3, 4, 1.63, 0.6*1.63), #6
    SwingEdge(4, 5, 1.63, 0.6*1.63)  #7
]

# BEWARE: ordering of edges ≠ ordering of add_edge calls!
# In this case I've ordered them manualy for this to match.
g = SimpleGraph(length(verticies))
for e in swingedges
    add_edge!(g, e.from, e.to)
end

# Define nodes/edges and network
odeverts = [ODEVertex(f! = vert, dim=2, sym=[:d, :w]) for vert in verticies]
staticedges = [StaticEdge(f! = edge, dim=1, sym=[:F]) for edge in swingedges]
swing_network! = network_dynamics(odeverts, staticedges, g)

# x0 determined from static solution
# x0 = [0.0 for i in 1:2*nv(g)]
x0 = [-0.12692637482862684, -1.3649456633810975e-6, 0.14641121510104085, 4.2191082676726005e-7, -0.24376507587890778, 1.567589744768255e-6, -0.12692637482862684, -1.3649456633810975e-6, 0.35120661043511864, 7.403907552948938e-7]

# define callback to save edge values (for plotting)
function save_edges(u, t, integrator)
    graph_data = integrator.f.f(u, integrator.p, integrator.t, GetGD)
    return copy(graph_data.gdb.e_array)
end

saved_edgevalues = SavedValues(Float64, Array{Float64, 1})
save_callback = SavingCallback(save_edges, saved_edgevalues)

# define problem
prob = ODEProblem(swing_network!, x0, (0., 6))

# define kill_first callback, at t=1.0s remove edge(2,4)
kill_first = PresetTimeCallback( 1.0, kill_line_affect(swing_network!, 2, 4))

# define load_watch_callback
load_watch_callback = watch_line_limit_callback(swing_network!)

sol = solve(prob, Tsit5(), callback=CallbackSet(save_callback, kill_first, load_watch_callback), dtmax=0.01)

# plot(sol, vars = idx_containing(swing_network!, :w))
# plot(sol, vars = idx_containing(swing_network!, :d))
plot_flow(saved_edgevalues, swing_network!)
