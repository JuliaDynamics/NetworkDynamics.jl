# [Injector Nodes](@id injector-nodes)

Lets consider the following subset of a Network:

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
we have a vertex of interested which is connected to two other vertices in the network via edges.
Generally, we follow the interface of having **potential like** outputs at the nodes and **flow like** outputs at the edges. I.e. the *potential on the nodes* depends on the *sum of flows through the edges* while the *flows through the edges* depend on the *potential on the adjacent nodes*.

The input-output structrue of this system looks something like this:
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

In typical NetworkDynamics modeling, the entire nodal dynamic will be subsumed in the single VertexModel.
However, it is quite common that there is some modulare structure within the vertex which
itself looks like multiple "flow injectors".

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

While not strictly necessary, it can be useful from a performance standpoint to split up those kind of vertex models in to "clusters", meaning we have a hub vertex and several injector vertices connecting them.
Notably, the "injector" models have a flipped input-output scheme, i.e. they take the potential of the hub as an direct input while outputting a flow.
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
```@docs
NetworkDynamics.LoopbackConnection
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
For demonstration purposes we'll model second vertex in two ways: as a single model enclsoing all three components and as separate injector nodes.

As allways, this is mainly a pedagogical example. For such a simple system, it is probably much cleaner to model it as a single vertex.
However thats not allwasy the case for very large networks with many complex vertex models!

### Prequisites
The first few components building blocks are identical to the docs on [ModelingToolkit Integration](@ref).
```@example injector 
using NetworkDynamics
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
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
## Part A: Modeling using Injector nodes and LoopbackConnection
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
Next we go for the resistor injector node:
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
We flipped the input-output scheme here: the injector takes the potential as input and outputs a flow.
However we notice, that the model above has no feed forward but instead a constraint. This is because per default, the VertexModel constructor will transform feed forwards to constraint states. This is a sane default, sinse most of our models won't be injector models. We can opt out of this behavior:
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

### Part B: Modeling using a single VertexModel
For comparison, we now model the same system using a single vertex model for the capacitor, resistor and inductor.
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
