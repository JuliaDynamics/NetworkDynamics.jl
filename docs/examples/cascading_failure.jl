#=
# Cascading Failure
This script reimplements the minimal example of a dynamic cascading failure
described in Schäfer et al. (2018) [1]. It is an example how to use callback
functions to change network parameters. In this case to disable certain lines.

[1] Schäfer, B., Witthaut, D., Timme, M., & Latora, V. (2018).
Dynamically induced cascading failures in power grids.
Nature communications, 9(1), 1-13.
https://www.nature.com/articles/s41467-018-04287-5

The system is modeled using swing equation and active power edges. The nodes are
characterized by the voltage angle `δ`, the active power on each line is symmetric
and a function of the difference between source and destination angle `δ_src - δ_dst`.
=#

using NetworkDynamics
using Graphs
using OrdinaryDiffEq
using DiffEqCallbacks
using Plots
import SymbolicIndexingInterface as SII

#=
For the nodes we define the swing equation. State `v[1] = δ`, `v[2] = ω`.
The swing equation has three parameters: `p = (P_ref, I, γ)` where `P_ref`
is the power setpopint, `I` is the inertia and `γ` is the droop or damping coeficcient.
=#
function swing_equation(dv, v, esum, p,t)
    P, I, γ = p
    dv[1] = v[2]
    dv[2] = P - γ * v[2] .+ esum[1]
    dv[2] = dv[2] / I
    nothing
end
odevertex  = ODEVertex(swing_equation; sym=[:δ, :ω], psym=[:P_ref, :I=>1, :γ=>0.1], depth=1)

#=
Lets define a simple purely active power line whose active power flow is
completlye determined by the connected voltage angles and the coupling constant
`K`.
We give an additonal parameter, the line limit, which we'll use later in the callback.
=#
function simple_edge(e, v_s, v_d, (K,), t)
    e[1] = K * sin(v_s[1] - v_d[1])
end
staticedge = StaticEdge(simple_edge; sym=:P, psym=[:K=>1.63, :limit=>1], coupling=AntiSymmetric())

#=
With the definition of the graph topology we can build the `Network` object:
=#
g = SimpleGraph([0 1 1 0 1;
                 1 0 1 1 0;
                 1 1 0 1 0;
                 0 1 1 0 1;
                 1 0 0 1 0])
swing_network = Network(g, odevertex, staticedge)

#=
For the parameters, we create the `NWParameter` object prefilled with default p values
=#
p = NWParameter(swing_network)
## vertices 1, 3 and 4 act as loads
p.v[(1,3,4), :P_ref] .= -1
## vertices 2 and 5 act as generators
p.v[(2,5), :P_ref] .= 1.5
nothing #hide #md

#=
We can use `find_fixpoint` to find a valid initial condition of the network
=#
u0 = find_fixpoint(swing_network, p)
nothing #hide #md

#=
In order to implement the line failures, we need to create a `VectorContinousCallback`.
In the callback, we compare the current flow on the line with the limit. If the limit is reached,
the coupling `K` is set to 0.

First we can define the affect function:
=#
function affect!(integrator, idx)
    println("Line $idx tripped at t=$(integrator.t)")
    p = NWParameter(integrator) # get indexable parameter object
    p.e[idx, :K] = 0
    auto_dt_reset!(integrator)
    nothing
end
nothing #hide #md

#=
The callback trigger condition is a bit more complicated. The straight forward version looks like this:
=#
function naive_condition(out, u, t, integrator)
    ## careful,  u != integrator.u
    ## therefore construct nwstate with Network info from integrator but u
    s = NWState(integrator, u, integrator.p, t)
    for i in eachindex(out)
        out[i] = abs(s.e[i,:P]) - s.p.e[1,:limit] # compare flow with limit for line
    end
    nothing
end
nothing #hide #md

#=
However, from a performacne perspectiv there are problems with this solution: on
every call, we need to perform symbolic indexing into the `NWState` object.
Symbolic indexing is not cheap, as it requires to gather meta data about the network.
Luckily, the `SymbolicIndexingInterface` package which powers the symbolic indexing
provids the lower level functions `getp` and `getu` which can be used to create and
cache accessors to the internal states.
=#
condition = let getlim = SII.getp(swing_network, EPIndex(1:ne(g), :limit)),
                getflow = SII.getu(swing_network, EIndex(1:ne(g), :P))
    function (out, u, t, integrator)
        ## careful,  u != integrator.u
        ## therefore construct nwstate with Network info from integrator but u
        s = NWState(integrator, u, integrator.p, t)
        out .= getlim(s) .- abs.(getflow(s))
        nothing
    end
end
nothing #hide #md

#=
We can combine affect and condition to form the callback.
=#
trip_cb = VectorContinuousCallback(condition, affect!, ne(g));

#=
However, there is another component missing. If we look at the powerflow on the
lines in the initial steady state
=#
u0.e[:, :P]
#=
We see that every flow is below the trip value 1.0. Therefor we need to add a distrubance
to the network. We do this by manually disabeling line 5 at time 1.
=#
trip_first_cb = PresetTimeCallback(1.0, integrator->affect!(integrator, 5));

#=
With those components, we can create the problem and solve it.
=#

prob = ODEProblem(swing_network, uflat(u0), (0,6), copy(pflat(p));
                  callback=CallbackSet(trip_cb, trip_first_cb))
Main.test_execution_styles(prob) # testing all ex styles #src
sol = solve(prob, Tsit5());

plot(sol; idxs=eidxs(sol,:,:P))

#=
Currently, plotting of "observed" states such as static edge states does not work,
because they depend on `e(t) = f(u(t), p(t), t)`, and currently the ODESolution
does not track `p` over time.

TODO: eventually SII will allow for "parameter timeseries" which track the parameters
=#
