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

# Find an inital condition

# ic_guess = rand(N * 3 ÷ 2)
# idx_ω = idx_containing(nd, :ω)
# ic_guess[idx_ω] .= 0.0
# x0 = find_valid_ic(nd, rand(N * 3 ÷ 2); p=p)

#fp = find_fixpoint(ndfp, p, rand(2N))

ndfp = network_dynamics(inertia, f_edge, g)
probfp = ODEProblem(ndfp, zeros(2N), tspan, p)
solfp = solve(probfp, Tsit5())
plot(solfp; vars=syms_containing(ndfp, :ω))
ic = find_fixpoint(ndfp, p, zeros(2N))
x0 = zeros(12)
x0[idx_containg]



# or try random
x0 = find_valid_ic(nd, randn(N * 3 ÷ 2); p=p)

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

