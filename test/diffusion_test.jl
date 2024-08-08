using Test
using NetworkDynamics
using Graphs
using LinearAlgebra
using SparseArrays
using OrdinaryDiffEq

printstyled("--- Diffusion dynamics --- \n"; bold=true, color=:white)

# We will work on a random graph:
g = barabasi_albert(10, 5)
L = laplacian_matrix(g)

N = nv(g)

# The test checks the diffusion dynamics on this graph. The analytic solution
# to this is given by the following function.

struct SolAnalytic
    L::Symmetric{Int,Array{Int,2}}
end
function (sa::SolAnalytic)(x0, p, t)
    exp(-t * sa.L) * x0
end

sol_analytic = SolAnalytic(Symmetric(Array(L)))

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

@inline function diffusion_vertex!(dv, v, edges, p, t)
    dv .= 0.0
    sum_coupling!(dv, edges) # Oriented sum of the incoming and outgoing edges
    nothing
end


odevertex = ODEVertex(; f=diffusion_vertex!, dim=1)
staticedge = StaticEdge(; f=diffusion_edge!, dim=1, coupling=:antisymmetric)

vertex_list = [odevertex for v in vertices(g)]
edge_list = [staticedge for e in edges(g)]

diff_network_st = network_dynamics(vertex_list, edge_list, g)

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

@testset "Promotion rules for static edges" begin
    @test_throws ErrorException odeedge = ODEEdge(staticedge) # We promote the static edge to an ODEEdge artifically

    @inline function promotable_diffusion_edge!(e, v_s, v_d, p, t)
        e[1] = v_s[1] - v_d[1]
        nothing
    end

    promotable_staticedge = StaticEdge(; f=promotable_diffusion_edge!,
                                       dim=1, coupling=:undirected)
    odeedge = ODEEdge(promotable_staticedge)

    ode_edge_list = [odeedge for e in edges(g)]

    diff_network_ode = network_dynamics(vertex_list, ode_edge_list, g)

    x0_ode = find_valid_ic(diff_network_ode, randn(nv(g) + 2 * ne(g)))
    dx0_ode = similar(x0_ode)
    diff_network_ode(dx0_ode, x0_ode, nothing, 0.0)
    @test NetworkDynamics.allocations(diff_network_ode, nothing) == 0
end

@inline function real_ode_edge!(de, e, v_s, v_d, p, t)
    de[1] = v_s[1] - v_d[1] - e[1]
    de[2] = v_d[1] - v_s[1] - e[2]
    nothing
end
odeedge = ODEEdge(; f=real_ode_edge!, dim=2, coupling=:fiducial, mass_matrix=0.0)

ode_edge_list = [odeedge for e in edges(g)]

const diff_network_ode = network_dynamics(vertex_list, ode_edge_list, g)

const x0_ode = find_valid_ic(diff_network_ode, randn(nv(g) + 2 * ne(g)))
const dx0_ode = similar(x0_ode)

diff_network_ode(dx0_ode, x0_ode, nothing, 0.0)
@test @allocated(diff_network_ode(dx0_ode, x0_ode, nothing, 0.0)) == 0

prob_L = ODEProblem(diff_network_L, x0_ode[1:N], (0.0, 5.0))

prob_st = ODEProblem(diff_network_st, x0_ode[1:N], (0.0, 5.0))
prob_ode = ODEProblem(diff_network_ode, x0_ode, (0.0, 5.0))


sol_L = solve(prob_L, Tsit5(); reltol=1e-5)
sol_st = solve(prob_st, Tsit5(); reltol=1e-5)
sol_ode = solve(prob_ode, Rodas5())

# These two are different code paths that we want to cover. If there is a type:
# mismatch between the internal GraphData and the argument being passed in, a
# new GraphData object is allocated of the right type. As this is essentially a
# fancy view into the Graph, the allocation is relatively mild.
diff_network_ode(dx0_ode, Array{Float32}(x0_ode), nothing, 0.0)
diff_network_ode(dx0_ode, x0_ode, nothing, 0.0)


max_L = [maximum(abs.(sol_L(t) .- sol_analytic(x0_ode[1:N], nothing, t))) for t in sol_L.t] |> maximum
max_st = [maximum(abs.(sol_st(t) .- sol_analytic(x0_ode[1:N], nothing, t))) for t in sol_L.t] |> maximum
max_ode = [maximum(abs.(sol_ode(t)[1:N] .- sol_analytic(x0_ode[1:N], nothing, t))) for t in sol_L.t] |> maximum
println("Maximum difference between analytic solution and...")
println("\t * and solution with explicit Laplacian: $max_L")
println("\t * and solution with static ND: $max_st")
println("\t * and solution with promoted static edges ND: $max_ode")

# We test that these helper function extract the right symbols:
syms_v = syms_containing(diff_network_ode, "v")
idx_v = idx_containing(diff_network_ode, "v")
@test all(diff_network_ode.syms[idx_v] .== syms_v)
