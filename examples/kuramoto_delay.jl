# # Network with time delays - Kuramoto model with delayed coupling

using DelayDiffEq
using Random
using Graphs
using NetworkDynamics
using Distributions

# A common modification of the [Kuramoto model](https://en.wikipedia.org/wiki/Kuramoto_model) is to include time-lags in the coupling function. In neuroscience this may be used to account for transmission delays along synapses connecting different neurons.
#
# In this tutorial we solve the system

# ```math
# \begin{aligned}
# \dot \theta_i(t) = \omega_i + \sum_{j=1}^N \kappa_{ji} (\theta_j(t - \tau_{ji}) - \theta_i(t))
# \end{aligned}
# ```
# where $\kappa_{ji}$ denotes the coupling strength of the edge from $j$ to $i$, $\tau_{ji}$ its time delay and $\theta_j(t - \tau_{ji})$ the past state of its source vertex
#
# For a general introduction to solving delay differential equations in Julia see this [tutorial]( https://diffeq.sciml.ai/stable/tutorials/dde_example/).
#
# To implement this in NetworkDynamics.jl a `StaticDelayEdge` has to be defined. Such a coupling function receives two additional arguments `h_v_s` and `h_v_d`.
# These are wrappers of the global history function that directly compute the history values of the local variables at the
# source and destination vertex respectively. The delay time should be
# passed like any other parameter. The `idxs` keyword argument of the history function can be used to access specific local variables.

function delay_coupling(e, θ_s, θ_d, h_θ_s, h_θ_d, p, t)
    τ, κ = p
    hist1 = h_θ_s(t - τ; idxs=1)
    e[1] = κ * sin(hist1 - θ_d[1])
    return nothing
end

edge = StaticDelayEdge(; f=delay_coupling, dim=1, coupling=:directed);

# The vertex dynamics are simple ODE vertices.

function kuramoto_vertex(dθ, θ, edges, p, t)
    ω = p
    dθ[1] = ω
    for e in edges
        dθ[1] += e[1]
    end
    return nothing
end

vertex = ODEVertex(; f=kuramoto_vertex, dim=1);

# For this example we use a complete graph. Bear in mind however that the data structures of Network Dynamics are best suited for sparse problems and might introduce some additional overhead for dense graphs.

N = 6
g = SimpleDiGraph(complete_graph(N))
nd = network_dynamics(vertex, edge, g)

Random.seed!(1)

ω = rand(nv(g)) # internal frequency
τ = rand(ne(g)) # time-lag per edge
κ = rand(ne(g)) # coupling strength per edge
p = (ω, [τ'; κ'])

θ₀ = rand(Uniform(0, 2π), N); # initial conditions

# Define random initial history function 

const past = rand(Uniform(0, 2π), N)
h(p, t; idxs=nothing) = typeof(idxs) <: Number ? past[idxs] : past;

# When constructing the DDEProblem lags may be specified.
# The solver will then track discontinuities arising from the evaluation of
# the history function and step to each of the discontinuities. 
# Since each multiplication combination of the lags may be connected to a
# discontinuity, this may be slow if many different lags are specified.
prob = DDEProblem(nd, θ₀, h, (0.0, 1.0), p; constant_lags=τ)
@time solve(prob, MethodOfSteps(BS3()); abstol=1e-10, reltol=1e-5); # ~50000 steps because of discontinuities
@time solve(prob, MethodOfSteps(BS3()); abstol=1e-10, reltol=1e-5);

# We recommend to solve with lowered absolute and relative error tolerances, since the Kuramoto system is highly multistable and the simulations may else produce very different results.
#
#
# The discontinuities arise from the initial history function and quickly get smoothed out (i.e. reduced in order) when the integration time is larger than the maximum lag. If the asymptotic behaviour is more interesting than the correct solution for a specific initial condition, it is possible to trade accuracy for computational speed by leaving the `constant_lags` undeclared.

fast_prob = DDEProblem(nd, θ₀, h, (0.0, 1.0), p)
@time solve(fast_prob, MethodOfSteps(BS3()); abstol=1e-10, reltol=1e-5); # ~200 steps
@time solve(fast_prob, MethodOfSteps(BS3()); abstol=1e-10, reltol=1e-5); 

# The `MethodOfSteps` algortihm extends an ODE solver to DDEs. For an overview of available solvers consult the manual of DifferentialEquations.jl. For example, for stiff systems, such as this one, it might be beneficial to use a stiff solver such as `TRBDF2`.

@time solve(fast_prob, MethodOfSteps(TRBDF2()); abstol=1e-10, reltol=1e-5);
@time solve(fast_prob, MethodOfSteps(TRBDF2()); abstol=1e-10, reltol=1e-5);

# Some further helpful comments for dealing within initial discontinuities in DDEs may be found in the [manual](https://jitcdde.readthedocs.io/en/stable/#dealing-with-initial-discontinuities) of the Python software JiTCDDE

using Plots
fast_prob = DDEProblem(nd, θ₀, h, (0.0, 10.0), p)
sol = solve(fast_prob, MethodOfSteps(BS3()); abstol=1e-10, reltol=1e-5, saveat=0.01)
plot(sol; fmt=:png)

