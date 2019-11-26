using Pkg
Pkg.activate(".")

using NetworkDynamics
using LightGraphs
using OrdinaryDiffEq
# using DiffEqOperators
using BenchmarkTools

N = 100
g = barabasi_albert(N, 5)

B = incidence_matrix(g, oriented=true)

struct kuramoto_dyn{T, T2, U}
    B::T
    B_trans::T2
    ω::U
    N::Int
end
function (dd::kuramoto_dyn)(dx, x, p, t)
    dx[1:N] .= dd.ω .- x[1:N] .- 5. .* dd.B * sin.(dd.B_trans * x[N+1:2N])
    dx[N+1:2N] .= x[1:N]
    nothing
end

ω = Array{Float64}(1:N) ./ N

kn = kuramoto_dyn(B, transpose(B), ω, N)
# Jv = JacVecOperator(kn, randn(nv(g)), nothing, 0.0)
kur_network_L = ODEFunction(kn) #, jac_prototype=Jv)

#Now for NetworkDynamics

@inline Base.@propagate_inbounds function kuramoto_edge!(e,v_s,v_d,p,t)
    e[1] = 5. * sin(v_s[2] - v_d[2])
    nothing
end

@inline Base.@propagate_inbounds function kuramoto_dedge!(de, e, v_s, v_d, p, t)
    de[1] = 100. * (5. * sin(v_s[2] - v_d[2]) - e[1])
    nothing
end

# @inline Base.@propagate_inbounds

@inline Base.@propagate_inbounds function kuramoto_vertex!(dv, v, e_s, e_d, p, t)
    dv[1] = p - v[1]
    for e in e_s
        dv[1] -= e[1]
    end
    for e in e_d
        dv[1] += e[1]
    end
    dv[2] = v[1]
    nothing
end


odevertex = ODEVertex(f! = kuramoto_vertex!, dim = 2, sym=[:ω, :ϕ])
staticedge = StaticEdge(f! = kuramoto_edge!, dim = 1)

vertex_list = [odevertex for v in vertices(g)]
edge_list = [staticedge for e in edges(g)]
ode_sd_edge_list = [ODEEdge(se) for se in edge_list]
ode_edge_list = [ODEEdge(f! = kuramoto_dedge!, dim = 1) for e in edges(g)]
p = vcat(ω, zeros(ne(g)))

kur_network_nd = network_dynamics(vertex_list,edge_list, g, p)
kur_network_eode = network_dynamics(vertex_list,ode_edge_list, g, p)


x0_L = 0.1 .* Array{Float64}(1:2N)
x0_nd = similar(x0_L)
x0_nd2 = similar(x0_L)

# Different ways to set the inital conditions:

# Explicitly addressing the ordering in x0_nd
ω_idx = idx_containing(kur_network_nd, :ω)
ϕ_idx = idx_containing(kur_network_nd, :ϕ)
x0_nd[ω_idx] .= x0_L[1:N]
x0_nd[ϕ_idx] .= x0_L[1+N:2N]

# A new way to interact with the state of the system is to call the
# network rhs with the signature (x, p, t, GetGD). This will return a GraphData
# object built on x, that allows reading and manipulating the underlying data.
# GraphData can also be called in various ways to get views into the underlying
# arrays:


# Using a GraphData built on top of x0_nd we can set the initial conditions this
# way:
gd_0 = kur_network_nd(x0_nd2, p, 0., GetGD)
gd_0(ViewV, :ω) .= x0_L[1:N]
gd_0(ViewV, :ϕ) .= x0_L[1+N:2N]

@assert x0_nd == x0_nd2

# Note that now the edge values gd_0(ViewE) are out of sync with the vertex
# values. There is no automatic updating as ViewV gives direct access to the
# underlying data. If we call gd_0 again on the same data, the internal edge
# cache get's updated and contains the appropriate values:

gd_0 = kur_network_nd(x0_nd2, p, 0., GetGD)


# For the ODE version we also have to set the edge values. We can use the
# updated gd_0 object for that.

x0_ode = Array{Float64}(1:(2N+ne(g)))
gd_ode = kur_network_eode(x0_ode, p, 0., GetGD)
gd_ode(ViewV) .= gd_0(ViewV)
gd_ode(ViewE) .= gd_0(ViewE)

# Now everything should be the same:

dx_L = similar(x0_L)
dx_nd = similar(x0_nd)
dx_ode = similar(x0_ode)

kur_network_nd(dx_nd, x0_nd, p, 0.)
kur_network_eode(dx_ode, x0_ode, p, 0.)
kur_network_L(dx_L, x0_L, nothing, 0.)

println("Accuracy")

println(dx_nd[ω_idx] .- dx_L[1:N] .|> abs |> maximum)
println(dx_ode[ω_idx] .- dx_L[1:N] .|> abs |> maximum)

println("Benchmarking")

@btime kur_network_nd(dx_nd, x0_nd, p, 0.)
@btime kur_network_L(dx_L, x0_L, nothing, 0.)
@btime kur_network_eode(dx_ode, x0_ode, p, 0.)

tspan = (0., 100.)
prob_nd = ODEProblem(kur_network_nd, x0_nd, tspan, p)
prob_ode = ODEProblem(kur_network_eode, x0_ode, tspan, p)
prob_L = ODEProblem(kur_network_L, x0_L, tspan, nothing)

println("Solving")

sol_nd = solve(prob_nd)
sol_nd = solve(prob_nd)
@time sol_nd = solve(prob_nd)

sol_L = solve(prob_L)
sol_L = solve(prob_L)
@time sol_L = solve(prob_L)

sol_ode = solve(prob_ode)
sol_ode = solve(prob_ode)
@time sol_ode = solve(prob_ode)

# Let's have a look at these solutions:

using Plots
pyplot() # GR backend doesn't display unicode

ptspan = (90., 100.)

plot(sol_nd, vars=ω_idx[1:10], tspan=ptspan)
plot(sol_L, vars=collect(1:10), tspan=ptspan)
plot(sol_ode, vars=ω_idx[1:10], tspan=ptspan)
