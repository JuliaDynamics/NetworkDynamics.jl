using NetworkDynamics
using LightGraphs
using OrdinaryDiffEq
using Plots

N = 50
k = 5
g = barabasi_albert(N, k)

function order_parameter

@inline function kuramoto_vertex!(dv, v, e_s, e_d, p, t)
    dv .= p
    oriented_edge_sum!(dv, e_s, e_d)
    nothing
end

@inline function kuramoto_edge!(e,v_s,v_d,p,t)
    # usually e, v_s, v_d are arrays, hence we use the broadcasting operator .
    e .= p * sin.(v_s - v_d)
    nothing
end

odevertex = ODEVertex(f! = kuramoto_vertex!, dim = 1)
staticedge = StaticEdge(f! = kuramoto_edge!, dim = 1)

v_pars = [1. * randn() for v in vertices(g)]
e_pars = [1. /3. for e in edges(g)]

parameters = (v_pars, e_pars)

kuramoto_network! = network_dynamics(odevertex, staticedge, g)

x0 = randn(nv(g)) # nv(g) - number of vertices in g
dx = similar(x0)

prob = ODEProblem(kuramoto_network!, x0, (0.,15), parameters)

sol = solve(prob, Tsit5())

plot(sol)

# Inspecting the Kuramoto order parameter

function order_parameter(s)
    θ  = 0
    for i in s
        θ += exp(im * i)
    end
    θ / length(s)
end

# First we get the indices of the vertex variables...
u_idx = idx_containing(kuramoto_network!, :v)
# ... then we compute the order parameter at each time step
plot(abs.(mapslices(order_parameter, sol[u_idx, :], dims = 1))[:])
