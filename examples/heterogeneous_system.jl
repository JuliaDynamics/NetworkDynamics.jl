# # Modeling a heterogeneous system
#
# One of the main purposes of NetworkDynamics.jl is to facilitate
# modeling coupled systems with heterogenities. This means that
# components can differ in their parameters as well as in their dynamics.

# ## Heterogenous parameters
# We start by setting up a simple system of Kuramoto oscillators.

using NetworkDynamics, OrdinaryDiffEq, Plots, Graphs

N = 8
g = watts_strogatz(N, 2, 0) # ring network

function kuramoto_edge!(e, θ_s, θ_d, K, t)
    e[1] = K * sin(θ_s[1] - θ_d[1])
end

function kuramoto_vertex!(dθ, θ, edges, ω, t)
    dθ[1] = ω
    sum_coupling!(dθ, edges)
end

vertex! = ODEVertex(; f=kuramoto_vertex!, dim=1, sym=[:θ])
edge!   = StaticEdge(; f=kuramoto_edge!, dim=1)
nd! = network_dynamics(vertex!, edge!, g);

# Introducing heterogeneous parameters is as easy as defining an array.
# Here the vertex parameters are hetereogeneous, while the edges share the same coupling
# parameter K.

ω = (collect(1:N) .- sum(1:N) / N) / N
K = 3.0
p = (ω, K); # p[1] vertex parameters, p[2] edge parameters

# Integrate and plot

x0 = (collect(1:N) .- sum(1:N) / N) / N
tspan = (0.0, 10.0)
prob = ODEProblem(nd!, x0, tspan, p)
sol = solve(prob, Tsit5())
plot(sol; ylabel="θ")

# ## Heterogenous dynamics

# Two paradigmatic modifications of the node model above are static nodes and nodes with
# inertia. A static node has no internal dynamics and instead fixes the variable at a
# constant value. A Kuramoto model with inertia consits of two interal variables leading to
# more complicated (and for many applications more realistic) local dynamics.

static! = StaticVertex(; f=(θ, edges, c, t) -> θ .= c, dim=1, sym=[:θ])


function kuramoto_inertia!(dv, v, edges, P, t)
    dv[1] = v[2]
    dv[2] = P - 1.0 * v[2]
    for e in edges
        dv[2] += e[1]
    end
end

inertia! = ODEVertex(; f=kuramoto_inertia!, dim=2, sym=[:θ, :ω]);


# Since now we model a system with hetereogeneous node dynamics we can no longer
# straightforwardly pass a single VertexFunction to `network_dynamics` but instead have to
# hand over an Array.

vertex_array    = Array{VertexFunction}([vertex! for i in 1:N])
vertex_array[1] = static!
vertex_array[5] = inertia! # index should correspond to the node's index in the graph
nd_hetero! = network_dynamics(vertex_array, edge!, g);


# Now we have to take a bit more care with defining initial conditions and parameters.
# For the static! node the initial condition has to match its parameter. For simplicity
# we use the same parameters as above.

x0[1] = ω[1];

# The node with inertia is two-dimensional, hence we need to specify two initial conditions.
# For the first dimension we keep the initial condition from above and insert! another one
# into `x0` at the correct index.

inertia_ic_2 = 5
insert!(x0, 6, inertia_ic_2);

# `x0[1:4]` holds ic for nodes 1 to 4, `x0[5:6]` holds the two
# initial conditions for node 5, `x0[7:9]` holds ic for nodes 6 to 8.


prob_hetero = ODEProblem(nd_hetero!, x0, tspan, p)
sol_hetero = solve(prob_hetero, Rodas4());

# For clarity we plot only the variables refering to the oscillator's angle θ and color
# them according to their type.

membership = ones(Int64, N)
membership[1] = 2
membership[5] = 3
nodecolor = [colorant"lightseagreen", colorant"orange", colorant"darkred"];
nodefillc = reshape(nodecolor[membership], 1, N);

vars = syms_containing(nd_hetero!, :θ);
plot(sol_hetero; ylabel="θ", vars=vars, lc=nodefillc)



# ## Components with algebraic constraints
#
# If one of the network components contains an algebraic as well as dynamical component,
# then there is the option to supply a mass matrix for the given component. In general this
# will look as follows:

function edgeA!(de, e, v_s, v_d, p, t)
    de[1] = f(e, v_s, v_d, p, t) # dynamic variable
    e[2]  = g(e, v_s, v_d, p, t) # static varibale
end

M = zeros(2, 2)
M[1, 1] = 1

nd_edgeA! = ODEEdge(; f=edgeA!, dim=2, coupling=:undirected, mass_matrix=M);


# This handles the second equations as `0 = M[2,2] * de[2] = g(e, v_s, v_d, p, t) - e[2]`.
#
# See the example kuramoto_plasticity.jl and the discussion on [github](https://github.com/pik-icone/NetworkDynamics.jl/issues/45#issuecomment-659491913) for more details.
