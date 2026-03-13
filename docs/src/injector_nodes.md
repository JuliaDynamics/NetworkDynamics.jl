# [Injector Nodes](@id injector-nodes)

In large network models, vertices often contain multiple internal components (e.g., generators, loads, storage devices). While these can be modeled as a single monolithic vertex model, splitting them into separate "injector nodes" connected via special loopback edges can offer performance and modularity advantages. This page explains the injector node pattern and demonstrates how to use `LoopbackConnection` edges.

## Concept

Let's consider the following subset of a Network:

```asciiart
   ⋮
  ┏┷━━┓      input
⋯─┨vₖ ┠───╮  aggr.  ┏━━━━━━━━┓
  ┗━━━┛   ╰─────╮   ┃        ┃
   ┏━━━┓       (+)──┨ Vertex ┃
⋯──┨vⱼ ┠────────╯   ┃        ┃
   ┗━┯━┛            ┗━━━━━━━━┛
     ⋮
```
we have a vertex of interest which is connected to two other vertices in the network via edges.
Generally, we follow the interface of having **potential like** outputs at the nodes and **flow like** outputs at the edges. I.e. the *potential on the nodes* depends on the *sum of flows through the edges* while the *flows through the edges* depend on the *potential on the adjacent nodes*.

The input-output structure of this system looks something like this:
```asciiart
                                 more edges
                                     △
n ⋯───╮             ╭────────────────┼────────────────╮             ╭───⋯ n
e     │             │      potential │ φ out          │             │     e
x  ┏━━▽━━━━━━━━━━━━━▽━━┓   ╔═════════△═════════╗   ┏━━▽━━━━━━━━━━━━━▽━━┓  x
t  ┃ EdgeModel         ┃   ║ VertexModel       ║   ┃ EdgeModel         ┃  t
   ┃ ẋ = f(x, φ, p, t) ┃   ║ ẋ = f(x, Φ, p, t) ║   ┃ ẋ = f(x, φ, p, t) ┃
n  ┃ Φ = g(x, φ, p, t) ┃   ║ φ = g(x, p, t)    ║   ┃ Φ = g(x, φ, p, t) ┃  n
o  ┗━━▽━━━━━━━━━━━━━▽━━┛   ╚═════════△═════════╝   ┗━━▽━━━━━━━━━━━━━▽━━┛  o
d     │        flow │ Φ out        ╭─┴─╮         flow │ Φ out       │     d
e ⋯───╯             ╰──────────────▷ + ◁──────────────╯             ╰───⋯ e
                                   ╰─△─╯
                                     │
                                 more edges
```
where, notably, only the edge models support feed forward behavior.

In typical NetworkDynamics modeling, the entire nodal dynamic is contained within a single VertexModel.
However, vertices often have modular internal structure consisting of multiple components that inject or draw flows.
For example, an electrical bus might have generators, loads, and storage devices all connected to it.

```asciiart
   ⋮                ┏━━━━━━━━━━━━━━━━━━┓
  ┏┷━━┓      input  ┃Vertex            ┃
⋯─┨vₖ ┠───╮  aggr.  ┃     ╭──────────╮ ┃
  ┗━━━┛   ╰─────╮   ┃  ╭──┤Injector A│ ┃
   ┏━━━┓       (+)──╂─(+) ╰──────────╯ ┃
⋯──┨vⱼ ┠────────╯   ┃  │  ╭──────────╮ ┃
   ┗━┯━┛            ┃  ╰──┤Injector B│ ┃
     ⋮              ┃     ╰──────────╯ ┃
                    ┗━━━━━━━━━━━━━━━━━━┛
```

While not strictly necessary, splitting these vertex models into "clusters" can improve performance and code organization. A cluster consists of a hub vertex and several injector vertices that connect to it. This approach can be particularly beneficial because:
- NetworkDynamics performs best when there are many identical components. Splitting components into smaller parts makes it more likely to have repeated, identical components.
- The model structure matches the physical modularity of the system.
- For ModelingToolkit models, large monolithic components can lead to higher compilation and symbolic simplification times compared to multiple smaller models.
Notably, injector models have a flipped input-output scheme compared to normal vertices: they take the hub's potential as a direct input and output a flow.
```asciiart
                Hub    Loopback  Satelites
              ╭──────╮╭────────╮╭──────────╮
   ⋮
  ┏┷━━┓                         ┏━━━━━━━━━━┓
⋯─┨vₖ ┠───╮   ┏━━━━━━┓    ╭─────┨Injector A┃
  ┗━━━┛   ╰───┨      ┠────╯     ┗━━━━━━━━━━┛
   ┏━━━┓      ┃ Σi=0 ┃
⋯──┨vⱼ ┠──────┨      ┠────╮     ┏━━━━━━━━━━┓
   ┗━┯━┛      ┗━━━━━━┛    ╰─────┨Injector B┃
     ⋮                          ┗━━━━━━━━━━┛

              ╰────────────────────────────╯
                     Vertex Cluster
```

To connect this kind of injector nodes, we use the special EdgeModel `LoopbackConnection`.
See the docstring below for a detailed explaination of the interfaces.
```@docs; canonical=false
LoopbackConnection
```

## Example
Following the other examples we'll showcase the feature on the a small electrical network.
The network to model looks like this:
```
                v1   Resistor   v2
                ●─←────███────→─●
                │            ╭──┼──╮
ideal v source (↗)           ┴  █  ⚕
                │            ┬  █  ⚕  C + R + L
                │            ╰──┼──╯
                ⏚               ⏚
```
For demonstration purposes we'll model second vertex in two ways: as a single model enclosing all three components and as separate injector nodes.

As always, this is mainly a pedagogical example. For such a simple system, it is probably much cleaner to model it as a single vertex.
However thats not always the case for very large networks with many complex vertex models!

### Prerequisites
The first few components building blocks are identical to the docs on [ModelingToolkit Integration](@ref).
```@example injector 
using NetworkDynamics
using ModelingToolkitBase
using ModelingToolkitBase: t_nounits as t, D_nounits as D
using SciCompDSL
using OrdinaryDiffEqTsit5
using CairoMakie

@mtkmodel NWTerminal begin
    @variables begin
        v(t), [description="Voltage at node"]
        i(t), [description="Current flowing into node"]
    end
end

@mtkmodel VoltageSource begin
    @components begin
       p = NWTerminal()
    end
    @parameters begin
        V = 1.0
    end
    @variables begin
        i(t), [description="produced current by ideal voltage source (observable)"]
    end
    @equations begin
        i ~ -p.i
        p.v ~ V
    end
end
@named vs = VoltageSource()
vs_vertex = VertexModel(vs, [:p₊i], [:p₊v])

@mtkmodel Resistor begin
    @components begin
        src = NWTerminal()
        dst = NWTerminal()
    end
    @parameters begin
        R = 1
    end
    @equations begin
        dst.i ~ (src.v - dst.v)/R
        src.i ~ -dst.i
    end
end
@named resistor = Resistor()
nothing # hide
```
## Part A: Modeling with Injector Nodes and LoopbackConnection

We'll model the circuit using separate components connected via loopback edges.
Since our capacitor has the equation
```math
\dot{u} = \frac{1}{C} i
```
it is a natural voltage source. We'll use it as the "hub" node which will be accompanied by two injector nodes for the resistor and inductor.
```@example injector
@mtkmodel Capacitor begin
    @components begin
        p = NWTerminal(;v=0)
    end
    @parameters begin
        C = 1.0
    end
    @equations begin
        D(p.v) ~ p.i / C
    end
end
@named cap = Capacitor()
hub_vertex = VertexModel(cap, [:p₊i], [:p₊v], name=:hub)
```
Next, we define the resistor as an injector node. Unlike the regular `Resistor` edge model, this takes voltage as input and outputs current:
```@example injector
@mtkmodel ResistorInjector begin
    @components begin
        p = NWTerminal()
    end
    @parameters begin
        R = 100.0
    end
    @equations begin
        p.i ~ p.v/R
    end
end
@named resistor_inj = ResistorInjector()
Rinj_vertex = VertexModel(resistor_inj, [:p₊v], [:p₊i], name=:R_injector)
```
Notice how we've flipped the interface: voltage becomes an input (`[:p₊v]`) and current becomes an output (`[:p₊i]`).
However, the model above has a constraint instead of feed-forward behavior. By default, the `VertexModel` constructor transforms feed-forward relationships into constraint states—a sensible default since most vertex models should not have feed-forward behavior. For injector nodes, we need to opt out of this transformation:
```@example injector
Rinj_vertex = VertexModel(resistor_inj, [:p₊v], [:p₊i], name=:R_injector, ff_to_constraint=false)
```

Next, we go for the Inductor injector node:
```@example injector
@mtkmodel InductorInjector begin
    @components begin
        p = NWTerminal()
    end
    @parameters begin
        L = 0.1
    end
    @equations begin
        D(p.i) ~ p.v / L
    end
end
@named inductor_inj = InductorInjector()
Linj_vertex = VertexModel(inductor_inj, [:p₊v], [:p₊i], name=:L_injector, ff_to_constraint=false)
```
Now, we need to define the connections between the hub and the injectors, and between voltage source and hub.
For the two injectors we define the special `LoopbackConnection` endges from injector to hub.
For the connection between voltage source and hub we use a regular edge model.
```@example injector
edges = [
    LoopbackConnection(potential=[:u], flow=[:i], src=:R_injector, dst=:hub, name=:R_loopback),
    LoopbackConnection(potential=[:u], flow=[:i], src=:L_injector, dst=:hub, name=:L_loopback),
    EdgeModel(resistor, [:src₊v], [:dst₊v], [:src₊i], [:dst₊i]; src=:vs, dst=:hub)
]
vertices = [vs_vertex, hub_vertex, Rinj_vertex, Linj_vertex]
nw = Network(vertices, edges; warn_order=false)
```
To simulate the system we use the default initial conditions. We can inspect the states of our network to see how the different variables span the components:
```@example injector
s0 = NWState(nw)
s0.v
```
With that knowlege, we can set the initial condition:
```@example injector
s0.v[:hub, :p₊v] = 0.0
s0.v[:L_injector, :p₊i] = 0.0
nothing #hide
```

From initial state we can simulate and plot the results:
```@example injector
prob = ODEProblem(nw, s0, (0.0, 10.0))
sol = solve(prob, Tsit5())

let
    fig = Figure()
    ax = Axis(fig[1,1])
    plot!(ax, sol; idxs=VIndex(:hub, :p₊v), label="Capacitor Voltage (Injector Nodes)", color=Cycled(1))
    plot!(ax, sol; idxs=VIndex(:L_injector, :p₊i), label="Inductor Current (Injector Nodes)", color=Cycled(2))
    axislegend(ax)
    fig
end
```

### Part B: Modeling with a Single VertexModel

For comparison, we now model the same system using a single monolithic vertex that contains all three components (capacitor, resistor, and inductor) internally.
```@example injector
@mtkmodel CRLModel begin
    @components begin
        cap = Capacitor()
        resistor = ResistorInjector()
        inductor = InductorInjector()
    end
    @variables begin
        v(t), [description="Voltage at node"]
        i(t)=0, [description="Current flowing into node"]
    end
    @equations begin
        0 ~ resistor.p.i + inductor.p.i + cap.p.i - i
        v ~ cap.p.v
        v ~ resistor.p.v
        v ~ inductor.p.v
    end
end
@named crl_model = CRLModel()
crl_vertex = VertexModel(crl_model, [:i], [:v], name=:CRL_vertex)
```
With that definition we can define the network again:
```@example injector
edges2 = [
    EdgeModel(resistor, [:src₊v], [:dst₊v], [:src₊i], [:dst₊i]; src=:vs, dst=:CRL_vertex)
]
vertices2 = [vs_vertex, crl_vertex]
nw2 = Network(vertices2, edges2; warn_order=false)

s0_2 = NWState(nw2)
s0_2[VIndex(:CRL_vertex, :v)] = 0.0
s0_2[VIndex(:CRL_vertex, :inductor₊p₊i)] = 0.0

prob2 = ODEProblem(nw2, s0_2, (0.0, 10.0))
sol2 = solve(prob2, Tsit5())
let
    fig = Figure()
    ax = Axis(fig[1,1])
    plot!(ax, sol; idxs=VIndex(:hub, :p₊v), label="Capacitor Voltage (Injector Nodes)", color=Cycled(1), alpha=0.5)
    plot!(ax, sol2; idxs=VIndex(:CRL_vertex, :v), label="Capacitor Voltage (Single Vertex)", color=Cycled(1), linestyle=:dash)
    plot!(ax, sol; idxs=VIndex(:L_injector, :p₊i), label="Inductor Current (Injector Nodes)", color=Cycled(2), alpha=0.5)
    plot!(ax, sol2; idxs=VIndex(:CRL_vertex, :inductor₊p₊i), label="Inductor Current (Single Vertex)", color=Cycled(2), linestyle=:dash)
    axislegend(ax)
    fig
end
```
As expected, we get identical results from both modeling approaches.
