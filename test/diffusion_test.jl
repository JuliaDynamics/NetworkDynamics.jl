using Test
using NetworkDynamics
using LightGraphs
using LinearAlgebra
using SparseArrays
using OrdinaryDiffEq

printstyled("--- Diffusion dynamics --- \n", bold=true, color=:white)

# We will work on a random graph:
g = barabasi_albert(10,5)
L = laplacian_matrix(g)

N = nv(g)

# The test checks the diffusion dynamics on this graph. The analytic solution
# to this is given by the following function.

struct SolAnalytic
    L::Symmetric{Int,Array{Int,2}}
end
function (sa::SolAnalytic)(x0, p, t)
    exp(- t * sa.L) * x0
end

sol_analytic = SolAnalytic(Symmetric(Array(L)))

# We now define an ODEFunction for differential equations
diff_network_L = ODEFunction((dx, x, p, t) -> dx .= - L * x)

# We can integrate this function using differential equations for some initial
# conditions:

x0 = rand(nv(g))
prob_L = ODEProblem(diff_network_L,x0,(0.,5.))
sol_L = solve(prob_L, Tsit5())

# This checks for the accuracy of the integration, but does not test our own
# code yet, hence is not listed as a test.
@assert all([all(sol_L(t) .- sol_analytic(x0, nothing, t) .< 10^-2) for t in sol_L.t])

# Now for NetworkDynamics

println("Building with simple API")

#This cas is covered by the simple api:
diff_network_simple_api = scalar_network_with_sum((x, p, t) -> 0., (x_s, x_t, p, t) -> x_t - x_s, g, nothing)

#Now for full complexity:

println("Building Static Network Dynamics")

# We define the fundamental edge and vertex functions:

# This tests that NetworkDynamics correectly reproduces the dynamics above.

@inline function diffusion_edge!(e,v_s,v_d,p,t)
    e .= v_s .- v_d
    nothing
end

@inline function diffusion_vertex!(dv, v, p, t, i)
    dv .= i
    nothing
end


odevertex = ODEVertex(f! = diffusion_vertex!, dim = 1)
staticedge = StaticEdge(f! = diffusion_edge!, dim = 1)

vertex_list = [odevertex for v in vertices(g)]
edge_list = [staticedge for e in edges(g)]

diff_network_st = network_dynamics(vertex_list, edge_list, oriented_edge_sum, g)

@test diff_network_st isa ODEFunction
@test diff_network_simple_api isa ODEFunction

x0 = rand(nv(g))
dx_L = similar(x0)
dx_st = similar(x0)
dx_sa = similar(x0)

@testset "Network dynamics static accuracy" begin
    for i in 1:30
        x0 = randn(nv(g))

        diff_network_L(dx_L, x0, nothing, 0.)
        diff_network_st(dx_st, x0, nothing, 0.)
        diff_network_simple_api(dx_sa, x0, nothing, 0.)

        @test isapprox(dx_L, dx_st)
        @test isapprox(dx_L, dx_sa)
    end
end

println("Building Static Network Dynamics with artifical ODE Edges")
# This promotes the static edges to dynamic edges, the resulting ODEFunction
# has a mass matrix that enforces the edge equations.

odeedge = ODEEdge(staticedge) # We promote the static edge to an ODEEdge artifically
ode_edge_list = [odeedge for e in edges(g)]

diff_network_ode = network_dynamics(vertex_list, ode_edge_list, oriented_edge_sum, g)

x0_ode = find_valid_ic(diff_network_ode, randn(nv(g) + ne(g)))
dx0_ode = similar(x0_ode)

diff_network_ode(dx0_ode, x0_ode, nothing, 0.)

prob_L = ODEProblem(diff_network_L,x0_ode[1:N],(0.,5.))

prob_st = ODEProblem(diff_network_st,x0_ode[1:N],(0.,5.))
prob_ode = ODEProblem(diff_network_ode,x0_ode,(0.,5.))
#
# Jv = diff_network_ode.jac_prototype
#
# Jv(dx0_ode, x0_ode, nothing, 0.)

sol_L = solve(prob_L, Tsit5())
sol_st = solve(prob_st, Tsit5())
sol_ode = solve(prob_ode, Rodas5())

 # These two are different code paths that we want to cover. If there is a type:
 # mismatch between the internal GraphData and the argument being passed in, a
 # new GraphData object is allocated of the right type. As this is essentially a
 # fancy view into the Graph, the allocation is relatively mild.
diff_network_ode(dx0_ode, Array{Float32}(x0_ode),nothing,0.)
diff_network_ode(dx0_ode, x0_ode,nothing,0.)


max_L = [maximum(abs.(sol_L(t) .- sol_analytic(x0_ode[1:N], nothing, t))) for t in sol_L.t] |> maximum
max_st = [maximum(abs.(sol_st(t) .- sol_analytic(x0_ode[1:N], nothing, t))) for t in sol_L.t] |> maximum
max_ode = [maximum(abs.(sol_ode(t)[1:N] .- sol_analytic(x0_ode[1:N], nothing, t))) for t in sol_L.t] |> maximum

println("Maximum difference to analytic solution with explicit Laplacian: $max_L")
println("Maximum difference to analytic solution with static ND: $max_st")
println("Maximum difference to analytic solution with fake ode edge ND: $max_ode")

# We test that these helper function extract the right symbols:
syms_v = syms_containing(diff_network_ode, "v")
idx_v = idx_containing(diff_network_ode, "v")
@test all( diff_network_ode.syms[idx_v] .== syms_v )
