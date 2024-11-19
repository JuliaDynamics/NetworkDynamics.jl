#=
# Modeling a heterogeneous system

This example can be dowloaded as a normal Julia script [here](@__NAME__.jl). #md

One of the main purposes of NetworkDynamics.jl is to facilitate
modeling coupled systems with heterogenities. This means that
components can differ in their parameters as well as in their dynamics.

## Heterogenous parameters

We start by setting up a simple system of Kuramoto oscillators.
=#

using NetworkDynamics, OrdinaryDiffEqTsit5, Plots, Graphs

N = 8
g = watts_strogatz(N, 2, 0) # ring network

function kuramoto_edge!(e, θ_s, θ_d, (K,), t)
    e[1] = K * sin(θ_s[1] - θ_d[1])
    nothing
end
edge! = EdgeModel(g=AntiSymmetric(kuramoto_edge!), outdim=1, psym=[:K=>3])
#-

function kuramoto_vertex!(dθ, θ, esum, (ω0,), t)
    dθ[1] = ω0 + esum[1]
    nothing
end
vertex! = VertexModel(f=kuramoto_vertex!, g=StateMask(1:1), sym=[:θ], psym=[:ω0], name=:kuramoto)
#-

nw = Network(g, vertex!, edge!)

#=
To assign parameters, we can create a `NWParameter` object based on the `nw` definition.
This parameter object will be pre-filled with the default parameters.
=#
p = NWParameter(nw)

#=
To set the vertex parameters, we can use indexing of the `p.v` field:
=#
ω = collect(1:N) ./ N
ω .-= sum(ω) / N
p.v[:, :ω0] = ω
nothing #hide #md

#=
Here, the index pairing `:, :ω` is used to index state ω for all node indices.

The parameter object contains information about the network structure. For the
actual problem definition we need to throw away this wrapper and use the
flat-vector representation of the parameters `pflat(p)`. Note that `pflat(p)`

Similarily, we could use `NWState(nw)` to create an indexable wrapper of the
initial state. However in this case we can also fill create the flat state
array manually:
=#
x0 = collect(1:N) ./ N
x0 .-= sum(x0) ./ N
tspan = (0.0, 10.0)
prob = ODEProblem(nw, x0, tspan, pflat(p))
Main.test_execution_styles(prob) # testing all ex styles #src
sol = solve(prob, Tsit5())
plot(sol; ylabel="θ", fmt=:png)

#=
## Heterogeneous dynamics

Two paradigmatic modifications of the node model above are static nodes and nodes with
inertia.
A static node has no internal states and instead fixes the variable at a
constant value.
=#
function static_g(out, u, p, t)
    out[1] = p[1]
    nothing
end
static! = VertexModel(g=static_g, outsym=[:θ], psym=[:θfix => ω[1]], name=:static)

#=
But wait! NetworkDynamics classified this as [`PureFeedForward`](@ref), because it cannot
distinguish between the function signatures
```
g(out, u, p, t)    # PureFeedForward
g(out, ins, p, t)  # NoFeedForward
```
and since `dim(u)=0` it wrongfully assumes that the latter is meant.
We can overwrite the classification by passing the ff keyword:
=#
static! = VertexModel(g=static_g, outsym=[:θ], psym=[:θfix => ω[1]], ff=NoFeedForward(), name=:static)

#=
A Kuramoto model with inertia consists of two internal variables leading to
more complicated (and for many applications more realistic) local dynamics.
=#
function kuramoto_inertia!(dv, v, esum, (ω0,), t)
    dv[1] = v[2]
    dv[2] = ω0 - 1.0 * v[2] + esum[1]
    nothing
end

inertia! = VertexModel(f=kuramoto_inertia!, g=1:1, sym=[:θ, :ω], psym=[:ω0], name=:inertia)

#=
Since now we model a system with heterogeneous node dynamics we can no longer
straightforwardly pass a single VertexModel to the `Network` constructor but
instead have to hand over an Array.
=#

vertex_array    = VertexModel[vertex! for i in 1:N]
vertex_array[1] = static!
vertex_array[5] = inertia! # index should correspond to the node's index in the graph
nw_hetero! = Network(g, vertex_array, edge!)

#=
Now we have to take a bit more care with defining initial conditions and parameters.

First, we can generate a `NWState` object based on the `nw_hetero!` object which
will be populated with the default values.
=#

state = NWState(nw_hetero!)

#=
The node with inertia is two-dimensional, hence we need to specify two initial conditions.
For the first dimension we keep the initial conditions from above and insert! another one into `x0` at
the correct index.

For the θ states we will use the same initial conditins as before:
=#
state.v[2:8,:θ] = x0[2:8]
nothing #hide #md

#=
We're still missing one initial condition: the second variable ω of the 5th vertex.
=#
state.v[5,:ω] = 5
nothing #hide #md

#=
The `NWState` object also contains a parameter object accessible via `state.p`.
The edge parameters are already filled with default values.
The vertex parameters can be copied from our old parmeter object `p`.
=#
state.p.v[2:8, :ω0] = p.v[2:8, :ω0]
nothing #hide #md

#=
For the problem construction, we need to convert the nested stuctures to flat arrays using the [`uflat`](@ref) and [`pflat`](@ref) methods.
=#
prob_hetero = ODEProblem(nw_hetero!, uflat(state), tspan, pflat(state))
Main.test_execution_styles(prob_hetero) # testing all ex styles #src
sol_hetero = solve(prob_hetero, Tsit5());
nothing #hide #md
plot(sol_hetero)

#=
For clarity we plot only the variables referring to the oscillator's angle θ and color
them according to their type.
=#

colors = map(vertex_array) do vertexf
    if vertexf.name == :kuramoto
        colorant"lightseagreen"
    elseif vertexf.name == :static
        colorant"orange"
    elseif vertexf.name == :inertia
        colorant"darkred"
    end
end

plot(sol_hetero; ylabel="θ", idxs=vidxs(1:8,:θ), lc=colors', fmt=:png)
