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

static = StaticVertex(; f=(θ, edges, c, t) -> θ[1] = c,
                      dim=1, sym=[:θ])

inertia = ODEVertex(; f=kuramoto_inertia!, dim=2, sym=[:θ, :ω])

f_edge = StaticEdge(; f=kuramoto_edge!, dim=1)



N = 9
# Simple ring topology
g = watts_strogatz(N, 2, 0.0)

v_arr = Array{VertexFunction}([inertia for _ in 1:N-1])
P     = shuffle!([ones((N - 1) ÷ 2); -ones((N - 1) ÷ 2)])
x0    = vcat([[rand(Uniform(-pi, pi)), rand(Uniform(-15, 15))] for _ in 1:N-1]...)

# first vertex is grid connection point
insert!(v_arr, 1, static)
# the microgrid is actually balanced
insert!(x0, 1, 0.0)

tspan = (0.0, 150.0)
p     = ([0.0; P], nothing)

nd = network_dynamics(v_arr, f_edge, g)
# i multiplied the initial conditions with 0.1 to get convergence, for a power grid model this would be important
# actually i don't really understand why it ended up in a desynchronized state, despite the grid connection :D
prob = ODEProblem(nd, x0 .* 0.1, tspan, p)
sol = solve(prob, Rodas4())


# Define node colors
begin
    nodefillc = [colorant"lightseagreen" for _ in 1:N-1]
    insert!(nodefillc, 1, colorant"orange")
    nodefillc = reshape(nodefillc, 1, N)
end

vars_θ = syms_containing(nd, :θ)
plot(sol; vars=vars_θ, lc=nodefillc)
vars_ω = syms_containing(nd, :ω)
plot!(sol; vars=vars_ω, lc=colorant"darkred")

