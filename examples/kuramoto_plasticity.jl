using Graphs
using OrdinaryDiffEq
using NetworkDynamics
using Plots

const N_plastic = 10 # number of nodes
k = 4  # average degree
g = barabasi_albert(N_plastic, k)

#=
  Berner, Rico, Eckehard Schöll, and Serhiy Yanchuk.
  "Multiclusters in Networks of Adaptively Coupled Phase Oscillators."
  SIAM Journal on Applied Dynamical Systems 18.4 (2019): 2227-2266.
=#


@inline function kuramoto_plastic_edge!(de, e, v_s, v_d, p, t)
    # The coupling function is modeled by a differential algebraic equation with mass matrix
    # 0 * de[1] = e[2] * sin(v_s[1] - v_d[1] + α) / N - e[1] is equivalent to
    # e[1] = e[2] * sin(v_s[1] - v_d[1] + α) / N

    de[1] = e[2] * sin(v_s[1] - v_d[1] + α) / N_plastic - e[1]
    de[2] = -ϵ * (sin(v_s[1] - v_d[1] + β) + e[2])

    nothing
end

@inline function kuramoto_plastic_vertex!(dv, v, edges, p, t)
    dv .= 0
    for e in edges
        dv .-= e[1]
    end
end


# Global parameters need to be const for type stability
const ϵ = 0.1
const α = 0.2π
const β = -.95π

# NetworkDynamics Setup
plasticvertex = ODEVertex(; f=kuramoto_plastic_vertex!, dim=1)
mass_matrix_plasticedge = zeros(2, 2)
mass_matrix_plasticedge[2, 2] = 1.0 # First variables is set to 0

plasticedge = ODEEdge(; f=kuramoto_plastic_edge!, dim=2, sym=[:e, :de], coupling=:undirected,
                      mass_matrix=mass_matrix_plasticedge);
kuramoto_plastic! = network_dynamics(plasticvertex, plasticedge, g)

# ODE Setup & Solution
x0_plastic     = vcat(randn(N_plastic), ones(4ne(g)))
tspan_plastic  = (0.0, 100.0)
params_plastic = (nothing, nothing)
prob_plastic   = ODEProblem(kuramoto_plastic!, x0_plastic, tspan_plastic, params_plastic)

sol_plastic = solve(prob_plastic, Rosenbrock23(); reltol=1e-6)

# Plotting
v_idx = idx_containing(kuramoto_plastic!, :v)

plot(sol_plastic; vars=v_idx, legend=false, ylabel="θ")

# Shows Coupling terms
e_idx = idx_containing(kuramoto_plastic!, :e)
#plot!(sol_plastic, vars=e_idx, legend=false, color=:black, linewidth=0.1)
