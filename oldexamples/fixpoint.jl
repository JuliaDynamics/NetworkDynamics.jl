using OrdinaryDiffEq
using NetworkDynamics
using NLsolve
using Graphs
using Plots
using GraphPlot

g = star_graph(4)
gplot(g)

function swing_eq!(dv, v, edges, P, t)
    dv[1] = v[2]
    dv[2] = P - 0.1 * v[2]
    for e in edges
        dv[2] += e[1]
    end
end

swing_vertex = ODEVertex(; f=swing_eq!, dim=2, sym=[:θ, :ω])

function load_eq!(dv, v, edges, P, t)
    dv[1] = P
    sum_coupling!(dv, edges)
    nothing
end

load_vertex = ODEVertex(; f=load_eq!, dim=1, mass_matrix=0, sym=[:θ])

function powerflow_eq!(e, v_s, v_d, K, t)
    e[1] = K * sin(v_s[1] - v_d[1])
end

powerflow_edge = StaticEdge(; f=powerflow_eq!, dim=1)

vertex_array = [swing_vertex, swing_vertex, load_vertex, load_vertex]
edge_array = [powerflow_edge for e in edges(g)]
nd = network_dynamics(vertex_array, edge_array, g)

K = 6.0
P = [1.0, 1.0, -1.0, -1.0]
p = (P, K)

u0 = find_fixpoint(nd, p, zeros(6))
tspan = (0.0, 100.0)
ode_prob = ODEProblem(nd, u0, tspan, p)
ode_sol = solve(ode_prob, Rosenbrock23())
plot(ode_sol; vars=syms_containing(nd, "θ"), legend=false)
