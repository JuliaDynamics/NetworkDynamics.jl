using DelimitedFiles
using SimpleWeightedGraphs, LightGraphs
using NetworkDynamics
using OrdinaryDiffEq
using DelayDiffEq
using Plots
using GraphPlot
using BenchmarkTools


### Constructing the graph
const N = 500
const k=340 # node degree = 2*R
const β=0
const ϵ=0.05
const sigma_inter=0.01
const sigma_intra=0.2
const a = .5
const ϕ=π/2 -0.1


ws = watts_strogatz(N,k,β)
A_ws= adjacency_matrix(ws)
A_0 = zeros(N, N)
#print(A_ws)
#gplot(ws,nodelabel=vertices(ws))




# combining the adjacency matrices for the outer wats-strogatz layers
L1_0 = vcat(A_ws,A_0)
L3_0 = vcat(A_0,A_ws)
L1_L3 = hcat(L1_0, L3_0)

g1 = SimpleWeightedDiGraph(L1_L3)

# adding the hub
add_vertex!(g1)

# adding the edges, where hub is connected to all other nodes
for i in 1:2N
    add_edge!(g1, 2N+1, i)
end
for i in 1:2N
    add_edge!(g1, i, 2N+1)
end

adjacency_matrix_complete = adjacency_matrix(g1)

#print(adjacency_matrix_complete)

#gplot(g1,nodelabel=vertices(g1))

print(ne(g1))
#collect(edges(g1))

edge_weights = getfield.(collect(edges(g1)), :weight);

# vertex function

@inline Base.@propagate_inbounds function fhn_electrical_vertex!(dv, v, e_s, e_d, p, t) #local term F(x)
    dv[1] = v[1] - v[1]^3 / 3 - v[2]
    dv[2] = (v[1] - a) * ϵ
    @inbounds for e in e_s
        dv[1] -= e[1]
        dv[2] -= e[2]
    end
    @inbounds for e in e_d
        dv[1] += e[1]
        dv[2] += e[2]
    end
    nothing
end


# edge functions

const H = [(1/ϵ)*cos(ϕ) (1/ϵ)*sin(ϕ); -sin(ϕ) cos(ϕ)] # rotation matrix
const G=[1/ϵ 0.0; 0.0 0.0]

@inline Base.@propagate_inbounds function electrical_edge_intra!(e, v_s, v_d, p, t)
    #H = [(1/ϵ)*cos(ϕ) (1/ϵ)*sin(ϕ); -sin(ϕ) cos(ϕ)]

    e[1] = sigma_intra * (v_s[1] - v_d[1])
    e[2] = sigma_intra * (v_s[2] - v_d[2])

    e[1] = H[1] * e[1] + H[2] * e[2]
    e[2] = H[3] * e[1] + H[4] * e[2]
    nothing
end

@inline Base.@propagate_inbounds function electrical_edge_inter_hub!(e, v_s, v_d, h_v_s, h_v_d, p, t)
    #G=[1/ϵ 0.0; 0.0 0.0]

    e[1] = p[1]/(2*N) * (h_v_s[1] - v_d[1])
    e[2] = p[2]/(2*N) * (h_v_s[2] - v_d[2])

    e[1]=G[1]*e[1] + G[2]*e[2]
    e[2]=G[3]*e[1] + G[4]*e[2]
    nothing
end

@inline Base.@propagate_inbounds function electrical_edge_inter!(e, v_s, v_d, h_v_s, h_v_d, p, t)
    #G=[1/ϵ 0.0; 0.0 0.0]

    e[1] = p[1] * (h_v_s[1] - v_d[1])
    e[2] = p[2] * (h_v_s[2] - v_d[2])

    e[1]=G[1]*e[1] + G[2]*e[2]
    e[2]=G[3]*e[1] + G[4]*e[2]
    nothing
end


edge_intra = StaticEdge(f! = electrical_edge_intra!, dim = 2, sym=[:u, :v])
edge_inter_outer = StaticDelayEdge(f! = electrical_edge_inter!, dim = 2, sym=[:u, :v])
edge_inter_hub = StaticDelayEdge(f! = electrical_edge_inter_hub!, dim = 2, sym=[:u, :v])

local_vertex = [ODEVertex(f! = fhn_electrical_vertex!, dim = 2, sym=[:u, :v]) for v in 1:(2N + 1)]

number_edges = ne(g1)

edge_list_combined_1 = [e % (k+1) == 0 ? edge_inter_outer : edge_intra for e in 1:(number_edges - (2*N))]
edge_list_combined_2 = [edge_inter_hub for e in ((number_edges - 2*N) + 1):number_edges]

edge_list_combined = vcat(edge_list_combined_1, edge_list_combined_2)

fhn_network! = network_dynamics(local_vertex, edge_list_combined, g1)

p = (nothing, sigma_inter * edge_weights, 0.0)
x0 = Array{Float64, 1}(vec([cos.(rand(2N+1)) sin.(rand(2N+1))]))
h(out, p, t) = (out .= 1.)
tspan = (0., 0.001)
prob2  = DDEProblem(fhn_network!, x0, h, tspan, p)
sol   = solve(prob2, MethodOfSteps(AutoTsit5(TRBDF2())))

### Benchmarking

dx0=similar(x0)
@benchmark(fhn_network!(dx0,x0,h,p,0))


### Plotting

plot(sol, vars = idx_containing(fhn_network!, :u), legend = false, ylim=(-3, 3))
plot(sol, vars = idx_containing(fhn_network!, :v), legend = false, ylim=(-5, 5))
