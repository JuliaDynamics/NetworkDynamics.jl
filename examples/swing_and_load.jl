using Graphs
using OrdinaryDiffEq
using NetworkDynamics
using Plots
using Random
using Distributions

function kuramoto_edge!(edge, v_s, v_d, p, t)
    # coupling strength σ = 9.0
    edge[1] = 9.0 * sin(v_s[1] - v_d[1])
    return nothing
end

function kuramoto_inertia!(dv, v, edges, P, t)
    # inertia α = 0.1
    dv[1] = v[2]
    dv[2] = P - 0.1 * v[2]
    for edge in edges
        dv[2] += edge[1]
    end
    return nothing
end

# A load constraint on the active power
function load_vertex(dv, v, edges, P, t)
    dv[1] = P
    for edge in edges
        dv[1] += edge[1]
    end
    return nothing
end

static = ODEVertex(; f=load_vertex, mass_matrix=0.0,
                   dim=1, sym=[:θ])

inertia = ODEVertex(; f=kuramoto_inertia!, dim=2, sym=[:θ, :ω])

f_edge = StaticEdge(; f=kuramoto_edge!, dim=1)

N = 8
# Simple ring topology
g = watts_strogatz(N, 2, 0.0)

v_arr = shuffle!([repeat([inertia], 4); repeat([static], 4)])
nd    = network_dynamics(v_arr, f_edge, g)
P     = map(x -> typeof(x) == typeof(inertia) ? 1.0 : -1.0, v_arr)
p     = (P, nothing)
tspan = (0.0, 150.0)

# Project the initial condition on the constraints 
x0 = find_valid_ic(nd, randn(N * 3 ÷ 2); p=p)
#... or try random ic
# x0 = rand((N * 3 ÷ 2)

prob = ODEProblem(nd, x0, tspan, p)
sol = solve(prob, Rodas4())

# Define node colors
begin
    nodefillc = map(x -> typeof(x) == typeof(inertia) ? colorant"lightseagreen" : colorant"orange", v_arr)
    nodefillc = reshape(nodefillc, 1, N)
end

vars_θ = syms_containing(nd, :θ)
plot(sol; vars=vars_θ, lc=nodefillc)
vars_ω = syms_containing(nd, :ω)
plot!(sol; vars=vars_ω, lc=colorant"darkred")


# Another way to construct the fixed point 
# ndfp = network_dynamics(inertia, f_edge, g)
# probfp = ODEProblem(ndfp, zeros(2N), tspan, p)
# solfp = solve(probfp, Tsit5(); reltol=1e-6, abstol=0.0)
# plot(solfp; vars=syms_containing(ndfp, :ω))
# fp = zeros(N * 3 ÷ 2)
# fp[idx_containing(nd, :θ)] .= solfp[end,idx_containing(ndfp, :θ)]
# prob0 = ODEProblem(nd, fp, tspan, p)
# sol0 = solve(prob0, Rodas4(); reltol=1e-6, abstol=0.0)
# plot(sol0)