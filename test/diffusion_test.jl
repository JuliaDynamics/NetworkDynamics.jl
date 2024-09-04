using Test
using NetworkDynamics
using Graphs
using LinearAlgebra
using SparseArrays
using OrdinaryDiffEq
using Random
Random.seed!(42)

printstyled("--- Diffusion dynamics --- \n"; bold=true, color=:white)

# We will work on a random graph:
g = barabasi_albert(10, 5)
L = laplacian_matrix(g)

N = nv(g)

# The test checks the diffusion dynamics on this graph. The analytic solution
# to this is given by the following function.

struct SolAnalytic
    L::LinearAlgebra.Symmetric{Int,Array{Int,2}}
end
function (sa::SolAnalytic)(x0, p, t)
    exp(-t * sa.L) * x0
end

sol_analytic = SolAnalytic(LinearAlgebra.Symmetric(Array(L)))

# We now define an ODEFunction for differential equations
diff_network_L = ODEFunction((dx, x, p, t) -> dx .= -L * x)

# We can integrate this function using differential equations for some initial
# conditions:

x0 = rand(nv(g))
prob_L = ODEProblem(diff_network_L, x0, (0.0, 5.0))
sol_L = solve(prob_L, Tsit5())

# This checks for the accuracy of the integration, but does not test our own
# code yet, hence is not listed as a test.
@assert all([all(sol_L(t) .- sol_analytic(x0, nothing, t) .< 10^-2) for t in sol_L.t])

# Now for NetworkDynamics


#Now for full complexity:

println("Building Static Network Dynamics...")

# We define the fundamental edge and vertex functions:

# This tests that NetworkDynamics correectly reproduces the dynamics above.

@inline function diffusion_edge!(e, v_s, v_d, p, t)
    e .= v_s .- v_d
    nothing
end

@inline function diffusion_vertex!(dv, v, esum, p, t)
    dv .= esum
    nothing
end


odevertex = ODEVertex(; f=diffusion_vertex!, dim=1, pdim=0)
staticedge = StaticEdge(; f=diffusion_edge!, dim=1, pdim=0, coupling=AntiSymmetric())

vertex_list = [odevertex for v in vertices(g)]
edge_list = [staticedge for e in edges(g)]

diff_network_st = Network(g, vertex_list, edge_list)

x0 = rand(nv(g))
dx_L = similar(x0)
dx_st = similar(x0)

@testset "Network dynamics static accuracy" begin
    for i in 1:30
        x0 = randn(nv(g))

        diff_network_L(dx_L, x0, nothing, 0.0)
        diff_network_st(dx_st, x0, nothing, 0.0)

        @test isapprox(dx_L, dx_st)

    end
end

println("Building Static Network Dynamics with artifical ODE Edges...")
# This promotes the static edges to dynamic edges, the resulting ODEFunction
# has a mass matrix that enforces the edge equations.

@inline function real_ode_edge!(de, e, v_s, v_d, p, t)
    de[1] = v_s[1] - v_d[1] - e[1]
    de[2] = v_d[1] - v_s[1] - e[2]
    nothing
end
odeedge = ODEEdge(; f=real_ode_edge!, dim=2, pdim=0, coupling=Fiducial(), mass_matrix=0)

ode_edge_list = [odeedge for e in edges(g)]

diff_network_ode = Network(g, vertex_list, ode_edge_list)

x0_ode = randn(nv(g) + 2 * ne(g))
dx0_ode = similar(x0_ode)

diff_network_ode(dx0_ode, x0_ode, nothing, 0.0)
@test @allocated(diff_network_ode(dx0_ode, x0_ode, nothing, 0.0)) == 0

prob_L = ODEProblem(diff_network_L, x0_ode[1:N], (0.0, 5.0))
prob_st = ODEProblem(diff_network_st, x0_ode[1:N], (0.0, 5.0))
prob_ode = ODEProblem(diff_network_ode, x0_ode, (0.0, 5.0))
Main.test_execution_styles(prob_st) # testing all ex styles #src
Main.test_execution_styles(prob_ode) # testing all ex styles #src

sol_L = solve(prob_L, Tsit5(); reltol=1e-5)
sol_st = solve(prob_st, Tsit5(); reltol=1e-5)
sol_ode = solve(prob_ode, Rodas5())

max_L = [maximum(abs.(sol_L(t) .- sol_analytic(x0_ode[1:N], nothing, t))) for t in sol_L.t] |> maximum
max_st = [maximum(abs.(sol_st(t) .- sol_analytic(x0_ode[1:N], nothing, t))) for t in sol_L.t] |> maximum
max_ode = [maximum(abs.(sol_ode(t)[1:N] .- sol_analytic(x0_ode[1:N], nothing, t))) for t in sol_L.t] |> maximum

@test max_L < 1e-6
@test max_st < 1e-6
@test max_ode < 1e-6

println("Maximum difference between analytic solution and...")
println("\t * and solution with explicit Laplacian: $max_L")
println("\t * and solution with static ND: $max_st")
println("\t * and solution with promoted static edges ND: $max_ode")
