# # Cascading Failure
# This script reimplements the minimal example of a dynamic cascading failure
# described in Schäfer et al. (2018) [1].
# In is an example how to use callback functions to change network parameters. In
# this case to disable certain lines.
#
# > [1] Schäfer, B., Witthaut, D., Timme, M., & Latora, V. (2018). Dynamically induced cascading failures in power grids. Nature communications, 9(1), 1-13. https://www.nature.com/articles/s41467-018-04287-5
#
# The system is modeled using swing equation and active power edges. The nodes are
# characterized by the voltage angle `δ`, the active power on each line is symmetric
# and a function of the difference between source and destination angle `δ_src - δ_dst`.

using NetworkDynamics
using Graphs
using OrdinaryDiffEqTsit5
using DiffEqCallbacks
using Plots
using Test #hide
import SymbolicIndexingInterface as SII

# For the nodes we define the swing equation. State `v[1] = δ`, `v[2] = ω`.
# The swing equation has three parameters: `p = (P_ref, I, γ)` where `P_ref`
# is the power setpopint, `I` is the inertia and `γ` is the droop or damping coeficcient.
#
# The output of the node is just the first state. `g=1` is a shorthand for `g=StateMask(1:1)`
# which implements a trivial output function `g` which just takes the first element of the
# state vector.

function swing_equation(dv, v, esum, p,t)
    P, I, γ = p
    dv[1] = v[2]
    dv[2] = P - γ * v[2] .+ esum[1]
    dv[2] = dv[2] / I
    nothing
end
vertex = VertexModel(f=swing_equation, g=1, sym=[:δ, :ω], psym=[:P_ref, :I=>1, :γ=>0.1])

# Lets define a simple purely active power line whose active power flow is
# completlye determined by the connected voltage angles and the coupling constant
# `K`.
# We give an additonal parameter, the line limit, which we'll use later in the callback.

function simple_edge(e, v_s, v_d, (K,), t)
    e[1] = K * sin(v_s[1] - v_d[1])
end
edge = EdgeModel(;g=AntiSymmetric(simple_edge), outsym=:P, psym=[:K=>1.63, :limit=>1])

# With the definition of the graph topology we can build the `Network` object:

g = SimpleGraph([0 1 1 0 1;
                 1 0 1 1 0;
                 1 1 0 1 0;
                 0 1 1 0 1;
                 1 0 0 1 0])
swing_network = Network(g, vertex, edge)

# For the parameters, we create the `NWParameter` object prefilled with default p values

p = NWParameter(swing_network)
# vertices 1, 3 and 4 act as loads
p.v[(1,3,4), :P_ref] .= -1
# vertices 2 and 5 act as generators
p.v[(2,5), :P_ref] .= 1.5

# We can use `find_fixpoint` to find a valid initial condition of the network

u0 = find_fixpoint(swing_network, p)

# In order to implement the line failures, we need to create a `VectorContinousCallback`.
# In the callback, we compare the current flow on the line with the limit. If the limit is reached,
# the coupling `K` is set to 0.
#
# First we can define the affect function:

function affect!(integrator, idx)
    println("Line $idx tripped at t=$(integrator.t)")
    p = NWParameter(integrator) # get indexable parameter object
    p.e[idx, :K] = 0
    auto_dt_reset!(integrator)
    save_parameters!(integrator)
    nothing
end

# There is one important aspect to this function: the `save_parameters!` call.
# In the callback, we change the parameters of the network, making the parameters time
# dependent. The flow on the line is a function `P(t) = f(u(t), p(t))`. Thus we need to
# inform the integrator, that a discrete change in parameters happend. With this, the
# solution object not only tracks `u(t)` but also `p(t)` and we may extract the observable `P(t)`
# directly.
#
# The callback trigger condition is a bit more complicated. The straight forward version looks like this:

function naive_condition(out, u, t, integrator)
    # careful,  u != integrator.u
    # therefore construct nwstate with Network info from integrator but u
    s = NWState(integrator, u, integrator.p, t)
    for i in eachindex(out)
        out[i] = abs(s.e[i,:P]) - s.p.e[1,:limit] # compare flow with limit for line
    end
    nothing
end

# However, from a performacne perspectiv there are problems with this solution: on
# every call, we need to perform symbolic indexing into the `NWState` object.
# Symbolic indexing is not cheap, as it requires to gather meta data about the network.
# Luckily, the `SymbolicIndexingInterface` package which powers the symbolic indexing
# provides the lower level functions `getp` and `getu` which can be used to create and
# cache accessors to the internal states.
#
# This still isn't ideal beacuse both `getlim` and `getflow` getters will create arrays
# within the callback. But is far better then resolving the flat state indices every time.

condition = let getlim = SII.getp(swing_network, epidxs(swing_network, :, :limit)),
                getflow = SII.getu(swing_network, eidxs(swing_network, :, :P))
    function (out, u, t, integrator)
        # careful,  u != integrator.u
        # therefore construct nwstate with Network info from integrator but u
        s = NWState(integrator, u, integrator.p, t)
        out .= getlim(s) .- abs.(getflow(s))
        nothing
    end
end

# We can combine affect and condition to form the callback.

trip_cb = VectorContinuousCallback(condition, affect!, ne(g));

# However, there is another component missing. If we look at the powerflow on the
# lines in the initial steady state

u0.e[:, :P]

# We see that every flow is below the trip value 1.0. Therefor we need to add a distrubance
# to the network. We do this by manually disabeling line 5 at time 1.

trip_first_cb = PresetTimeCallback(1.0, integrator->affect!(integrator, 5));

# With those components, we can create the problem and solve it.

prob = ODEProblem(swing_network, uflat(u0), (0,6), copy(pflat(p));
                  callback=CallbackSet(trip_cb, trip_first_cb))
sol = solve(prob, Tsit5());
# we want to test the reconstruction of the observables # hide
@test all(!iszero, sol(sol.t; idxs=eidxs(sol,:,:P))[begin]) # hide
@test all(iszero, sol(sol.t; idxs=eidxs(sol,:,:P))[end][[1:5...,7]]) # hide
nothing #hide

# Through the magic of symbolic indexing we can plot the power flows on all lines:

plot(sol; idxs=eidxs(sol,:,:P))

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
