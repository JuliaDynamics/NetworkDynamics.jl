using Graphs
using OrdinaryDiffEq
using Plots

N = 10
g = barabasi_albert(N, 5)

B = incidence_matrix(g; oriented=true)
L = laplacian_matrix(g)

P = randn(N)

P .-= sum(P) / N

# Kuramoto Dynamics, read up on incidence matrices to understand the code below.

struct kur_dyn{T,T2}
    B::T
    B_trans::T2
end
function (dd::kur_dyn)(dx, x, p, t)
    dx .= p .- dd.B * sin.(dd.B_trans * x)
    nothing
end

kn = kur_dyn(B, transpose(B))
kur_prob = ODEProblem(kn, randn(N), (0.0, 10.0), P)

kur_sol = solve(kur_prob, Tsit5())

plot(kur_sol)

# Linear diffusion dynamics, read up on Graph Laplacian matrices

struct diff_dyn{T}
    L::T
end
function (dd::diff_dyn)(dx, x, p, t)
    dx .= p .- dd.L * x
    nothing
end

dd = diff_dyn(L)
dd_prob = ODEProblem(dd, randn(N), (0.0, 10.0), P)

dd_sol = solve(dd_prob, Tsit5())

plot(dd_sol)

# diffusion dynamics with a x_0 node coupling in from the outside.

struct diff_dyn_x0{T}
    L::T
end
function (dd::diff_dyn_x0)(dx, x, par, t)
    P, x_0 = par
    dx .= P .- dd.L * x
    dx[1] += x_0(t) - x[1]
    nothing
end

x_0(t) = sin(t)# + 0.4 * sin(2t + 0.1) - 0.7 * sin(3t + 0.2)

ddx0 = diff_dyn_x0(L)
ddx0_prob = ODEProblem(ddx0, randn(N), (0.0, 50), (P, x_0))


ddx0_sol = solve(ddx0_prob, Tsit5())

plot(ddx0_sol)
