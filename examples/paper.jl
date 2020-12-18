# Chapter 1 - Pure Julia

using LightGraphs, OrdinaryDiffEq

N = 10
g = watts_strogatz(N,  2, 0.)
const B = incidence_matrix(g, oriented=true)
const B_t = transpose(B)

function kuramoto_network!(dθ, θ, ω, t)
    dθ .= ω .- 5 .* (B * sin.(B_t * θ))
    return nothing
end

ω = (collect(1:N) .- (sum(1:N) / N) ) / N
x0 = (collect(1:N) .- (sum(1:N) / N) ) / N
tspan = (0., 4.)

prob = ODEProblem(kuramoto_network!, x0, tspan, ω)

sol = solve(prob, Tsit5())


# Chapter 2 - Homogeneous ND

function kuramoto_edge!(e, θ_s, θ_d, σ, t)
    e .=  σ .* sin.(θ_s .- θ_d)
end
function kuramoto_vertex!(dθ, θ, edges, ω, t)
    dθ .= ω
    for e in edges
        dθ[1] += e[1]
    end
end

using NetworkDynamics

vertex! = ODEVertex(f! = kuramoto_vertex!, dim = 1, sym=[:θ])
edge!   = StaticEdge(f! = kuramoto_edge!, dim = 1)
nd! = network_dynamics(vertex!, edge!, g)

vertexp = ω
edgep   = 5.
p = (vertexp, edgep)

nd_prob = ODEProblem(nd!, x0, tspan, p)
nd_sol = solve(nd_prob, Tsit5())

using Plots

#plotlyjs()
plot(nd_sol, ylabel="θ", color_palette = :seaborn_dark)

savefig("chaos_homogeneous.pdf")

# Chapter 3 - Heterogeneous ND

# Hidden
membership = Int64.(ones(N))
membership[1] = 2
membership[N ÷ 2] = 3
nodecolor = [colorant"lightseagreen", colorant"orange", colorant"darkred"];
# membership color
nodefillc = nodecolor[membership];
nodefillc = reshape(nodefillc, 1, N);

function kuramoto_inertia!(dv, v, edges, p, t)
    dv[1] = v[2]
    dv[2] = p - v[2]
    for e in edges
        dv[1] += e[1]
    end
end

inertia! = ODEVertex(f! = kuramoto_inertia!, dim = 2, sym= [:θ, :ω])

static! = StaticVertex(f! = (θ, edges, c, t) -> θ .= c, dim = 1, sym = [:θ])

function kuramoto_edge!(e, θ_s, θ_d, σ, t)
    e[1] = σ * sin(θ_s[1] - θ_d[1])
end

vertex_array    = Array{VertexFunction}( [vertex! for v in vertices(g)])
vertex_array[1] = inertia!
vertex_array[N ÷ 2] = static!
nd_hetero!      = network_dynamics(vertex_array, edge!, g);

# Parameters and inital conditions

insert!(x0, 2, 3.) # add initial condition for inertia vertex
prob_hetero = ODEProblem(nd_hetero!, x0, tspan, p);
sol_hetero = solve(prob_hetero, Rodas4());



vars = syms_containing(nd_hetero!, :θ);
plot(sol_hetero, ylabel="θ", vars=vars, lc = nodefillc)
savefig("chaos_heterogeneous.pdf")


# Chapter 4 - Fancy ND (Delays)


function kuramoto_delay_edge!(e, v_s, v_d, h_v_s, h_v_d, p, t)
    e[1] = p * sin(v_s[1] - h_v_d[1])
    nothing
end
kdedge! = StaticDelayEdge(f! = kuramoto_delay_edge!, dim=1)
nd_delay! = network_dynamics(vertex_array, kdedge!, g)

using DelayDiffEq

h(out, p, t) = (out .= x0)
τ = 0.12
tspan = (0.,4.)
p = (vertexp,  edgep, τ)
prob_delay = DDEProblem(nd_delay!, x0, h, tspan, p)

sol_delay = solve(prob_delay, MethodOfSteps(Rodas4(autodiff=false)))

plot(sol_delay, ylabel="θ", vars=vars, lc = nodefillc)
savefig("chaos_delay.pdf")
# Chapter  5 - Fancy ND II (Callbacks)

using DiffEqCallbacks


θ_idxs = idx_containing(nd_delay!, :θ)
function condition(out, u, t, integrator)
    out .= (u[θ_idxs] .- 0.5) .* (u[θ_idxs] .+ 0.5)
    nothing
end
function affect!(integrator, idx)
    stable_edges = map(e -> idx ∉ e, Pair.(edges(g)))
    integrator.p = (integrator.p[1], stable_edges .* integrator.p[2], integrator.p[3])
    nothing
end
cb = VectorContinuousCallback(condition, affect!, 10)
prob_cb = remake(prob_delay, p =(vertexp,  edgep .* ones(N), τ))

sol_cb = solve(prob_cb, MethodOfSteps(Rodas4(autodiff=false)), callback=cb)

plot(sol_cb, ylabel="θ", vars=vars, lc = nodefillc)
hline!([-.5], color = [:black], width = [1.], line=[:dot], label="")
hline!([.5], color = [:black], width = [1.], line=[:dot], label="")

savefig("chaos_callback.pdf")
