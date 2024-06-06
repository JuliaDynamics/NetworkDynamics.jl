using Pkg
Pkg.activate(@__DIR__)
using Revise

using NetworkDynamics
using Graphs
using OrdinaryDiffEq
using BenchmarkTools

#### First without using NetworkDynamics

### Defining the graph
N = 100 # nodes
g = barabasi_albert(N, 5)
B = incidence_matrix(g; oriented=true)

struct kuramoto_dyn{T,T2,U}
    B::T
    B_trans::T2
    ω::U
    N::Int
end

# callable struct with differential equation of Kurmaoto-Oscialltor (2-dim)
# x and dx arrays containing 2 variables (x[1:N] : ω , x[N+1:2N] : ϕ)
function (dd::kuramoto_dyn)(dx, x, p, t)
    dx[1:N] .= dd.ω .- x[1:N] .- 5.0 .* dd.B * sin.(dd.B_trans * x[N+1:2N])
    dx[N+1:2N] .= x[1:N]
    nothing
end

# constructing ordered eigenfrequencies of the N vertices
ω = Array{Float64}(1:N) ./ N

kn = kuramoto_dyn(B, transpose(B), ω, N)
# Jv = JacVecOperator(kn, randn(nv(g)), nothing, 0.0)
# callable struct kn is handed over to ODEFunction
kur_network_L = ODEFunction(kn) #, jac_prototype=Jv)

### Now for NetworkDynamics

# StaticEdge fct
Base.@propagate_inbounds function kuramoto_edge!(e, v_s, v_d, p, t)
    # coupling strength K=5
    e[1] = 5.0 * sin(v_s[2] - v_d[2])
    nothing
end

# StaticEdge fct
Base.@propagate_inbounds function promotable_kuramoto_edge!(e, v_s, v_d, p, t)
    # coupling strength K=5
    e[1] = 5.0 * sin(v_s[2] - v_d[2])
    e[2] = 5.0 * sin(v_d[2] - v_s[2])
    nothing
end

# ODEdedge fct
Base.@propagate_inbounds function kuramoto_dedge!(de, e, v_s, v_d, p, t)
    de[1] = 100.0 * (5.0 * sin(v_s[2] - v_d[2]) - e[1])
    de[2] = 100.0 * (5.0 * sin(v_d[2] - v_s[2]) - e[2])
    nothing
end


# ODEVertex function
Base.@propagate_inbounds function kuramoto_vertex!(dv, v, edges, p, t)
    dv[1] = p - v[1]
    for e in edges
        dv[1] += e[1]
    end
    dv[2] = v[1]
    nothing
end

### Constructing the network dynamics

# StaticEdge case
odevertex = ODEVertex(; f=kuramoto_vertex!, dim=2, sym=[:ω, :ϕ])
staticedge = StaticEdge(; f=kuramoto_edge!, dim=1)
ode_static_edge = ODEEdge(; f=promotable_kuramoto_edge!, dim=2, coupling=:fiducial)
#odeedge = ODEEdge(f = real_kuramoto_dedge!, dim = 2, coupling = :fiducial, mass_matrix = 0.)
# we can also create lists of the vertices/edges, helpful if vertices/edges have different DE's
vertex_list = [odevertex for v in vertices(g)]
edge_list = [staticedge for e in edges(g)]

# StaticEdge with artificially promoted StaticEdges --> ODEEdges
ode_sd_edge_list = [ode_static_edge for se in edge_list]
#ode_sd_edge_list = [ODEEdge(se) for se in edge_list]
# ODEEdges with ODEEdgefct
ode_edge_list = [ODEEdge(; f=kuramoto_dedge!, dim=2, coupling=:fiducial) for e in edges(g)]
p = (ω, nothing)

# now we create the NetworkDynamics objects
# StaticEdge case
kur_network_hom = network_dynamics(odevertex, staticedge, g)
# StaticEdge case in list form
kur_network_nd = network_dynamics(vertex_list, edge_list, g)
# promotes StaticEdges --> ODEEdges
kur_network_e_static_ode = network_dynamics(vertex_list, ode_sd_edge_list, g)
# ODEEdges with ODEEdgefct
kur_network_eode = network_dynamics(vertex_list, ode_edge_list, g)

### Initial conditions

x0_L = 0.1 .* Array{Float64}(1:2N)
x0_nd = similar(x0_L)
x0_nd2 = similar(x0_L)

# different ways to set the inital conditions:

# explicitly addressing the ordering in x0_nd:
# ϕ-indices are even, ω-indices are uneven
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
gd_0 = kur_network_nd(x0_nd2, p, 0.0, GetGD)
gs_0 = kur_network_nd(GetGS)

view_v(gd_0, gs_0, :ω) .= x0_L[1:N]
view_v(gd_0, gs_0, :ϕ) .= x0_L[1+N:2N]

# we can also skip this step and call view_v with all the required arguments:

x0_nd3 = similar(x0_nd)

ω_view = view_v(kur_network_nd, x0_nd3, p, 0.0, :ω)
ϕ_view = view_v(kur_network_nd, x0_nd3, p, 0.0, :ϕ)

ω_view .= x0_L[1:N]
ϕ_view .= x0_L[1+N:2N]

# @assert cond [text]
# throw an AssertionError if cond is falsex0_nd == x0_nd2
@assert x0_nd == x0_nd3

# Note that now the edge values gd_0(ViewE) are out of sync with the vertex
# values. There is no automatic updating as ViewV gives direct access to the
# underlying data. If we call gd_0 again on the same data, the internal edge
# cache get's updated and contains the appropriate values:
gd_0 = kur_network_nd(x0_nd2, p, 0.0, GetGD)


# For the ODE version we also have to set the edge values. We can use the
# updated gd_0 object for that.
# ne(g) - number of edges
# array with initial conditions for vertices und edges
x0_ode = Array{Float64}(1:(2N+ne(g)*2))
gd_ode = kur_network_eode(x0_ode, p, 0.0, GetGD)
gs_ode = kur_network_eode(GetGS)
view_v(gd_ode, gs_ode) .= view_v(gd_0, gs_0)
view_e(gd_ode, gs_ode) .= view_e(gd_0, gs_0)

# Now everything should be the same:

dx_L = similar(x0_L)
dx_nd = similar(x0_nd)
dx_ode = similar(x0_ode)

# callable structs, calling the ODEFunction
kur_network_nd(dx_nd, x0_nd, p, 0.0)
kur_network_eode(dx_ode, x0_ode, p, 0.0)
kur_network_L(dx_L, x0_L, nothing, 0.0)

println("Accuracy")

println(dx_nd[ω_idx] .- dx_L[1:N] .|> abs |> maximum)
println(dx_ode[ω_idx] .- dx_L[1:N] .|> abs |> maximum)

println("Benchmarking")

@btime kur_network_nd(dx_nd, x0_nd, p, 0.0)
@btime kur_network_L(dx_L, x0_L, nothing, 0.0)
@btime kur_network_eode(dx_ode, x0_ode, p, 0.0)
@btime kur_network_hom(dx_nd, x0_nd, p, 0.0)

### Simulation

tspan = (0.0, 100.0)
prob_nd = ODEProblem(kur_network_nd, x0_nd, tspan, p)
prob_ode = ODEProblem(kur_network_eode, x0_ode, tspan, p)
prob_L = ODEProblem(kur_network_L, x0_L, tspan, nothing)

println("Solving")

# Tsit5 - non-stiff solver
sol_nd = solve(prob_nd, Tsit5())
sol_nd = solve(prob_nd, Tsit5())
@time sol_nd = solve(prob_nd, Tsit5())

sol_L = solve(prob_L, Tsit5())
sol_L = solve(prob_L, Tsit5())
@time sol_L = solve(prob_L, Tsit5())

sol_ode = solve(prob_ode, Tsit5())
sol_ode = solve(prob_ode, Tsit5())
@time sol_ode = solve(prob_ode, Tsit5())

### Plotting

using Plots

ptspan = (90.0, 100.0)

plot(sol_nd; vars=ω_idx[1:10], tspan=ptspan)
plot(sol_L; vars=collect(1:10), tspan=ptspan)
plot(sol_ode; vars=ω_idx[1:10], tspan=ptspan)
