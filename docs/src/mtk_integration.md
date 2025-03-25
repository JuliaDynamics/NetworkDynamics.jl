# ModelingToolkit Integration

NetworkDynamics.jl is compatible with [`ModelingTookit.jl`](https://github.com/SciML/ModelingToolkit.jl) (MTK).
The general idea is to use MTK to define *component models* (i.e. edge and vertex dynamics)
which are then connected on network scale using NetworkDynamics.

The main entry point for this interop are the constructors
```julia
VertexModel(::ODESystem, inputs, outputs)
EdgeModel(::ODESystem, srcin, dstin, [srscout], dstout)
```
whose docstrings can be found in the [Component Models with MTK](@ref) section in the API.

These constructors will:
- transforming the states marked as input to parameters and `structural_simplify`ing the system,
- generating the `f` and `g` functions,
- generate code for observables,
- port all supported [Metadata](@ref) from MTK symbols to component symbols and
- output a `Vertex-`/`EdgeModel` function compatible with NetworkDynamics.jl.

The main usecase for this feature is when you want to build relatively complex
component models but interconnect them in a very homogeneous way (i.e. having the
same output/input pairings in the whole system).

In theory, you can achieve everything you want to do with plain MTK. The idea of
combining the two is, that NetworkDynamics offers *far less flexibility* when in
comes to interconnection of subsystems on the network level. This might allow ND
to exploit more knowledge of the structure without very expensive operations such
as tearing of thousands of equations.

!!! warning
    ModelingToolkit is a fast paced library with lots of functionality and ever
    growing complexity. As such the provided interface is kinda experimental.
    Some features of MTK are straight up unsupported, for example events within
    models or delay differential equations.

## RC-Circuit Example
In good [MTK tradition](https://docs.sciml.ai/ModelingToolkit/stable/tutorials/acausal_components/), this feature will be explained along a simple RC circuit example.
The [Gas Network Example](@ref gas-example) or [Initialization Tutorial](@ref init-tutorial) also showcase the MTK constructors.

The system to model is 2 node, 1 edge network. The node output states are the voltage (to ground), the edge output sates are the currents at both ends.
```

ideal v source      Resistor     Capacitor
            v1 o─←────MMM────→─o v2
               │               ┴ 
              (↗)              ┬
               │               │
               ⏚               ⏚
```
Obviously there is no need in modeling such a small system using NetworkDynamics, however
the method extends quite easily to construct large electrical networks reusing the same
fundamental building blocks.

```@example mtk
using NetworkDynamics
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using OrdinaryDiffEqTsit5
using CairoMakie
```

All our components have "terminals", which have a voltage and current. We don't use the `@connector` from MTK here because our pins mark the interface towards the network and do not follow the MTK connector semantics.

```@example mtk
@mtkmodel NWTerminal begin
    @variables begin
        v(t), [description="Voltage at node"]
        i(t), [description="Current flowing into node"]
    end
end
nothing #hide
```

An ideal voltage source is just a model which pins its output voltage to a fixed parameter.
The source ejects whatever current is necessary. We introduce another variable `i(t)`
to "capture" this current. This variable will be removed during structural simplify, but will
be available for plotting through the [Observables](@ref) mechanism.
The `VertexModel` can be generated from an `ODESystem` by providing names of the input and output states:

```@example mtk
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
vs_vertex = VertexModel(vs, [:p₊i], [:p₊v]; vidx=1)
```

A capacitor is a slightly more complicated model. Its voltage is defined as an differential
equation based on the inflowing current.

```@example mtk
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
cap_vertex = VertexModel(cap, [:p₊i], [:p₊v], vidx=2)
```

For the resistor we need two pins, one for the `src` and one for the `dst` side.
The equations are straight forward.

```@example mtk
@mtkmodel Resistor begin
    @components begin
        src = NWTerminal()
        dst = NWTerminal()
    end
    @parameters begin
        R = 1.0
    end
    @equations begin
        dst.i ~ (src.v - dst.v)/R
        src.i ~ -dst.i
    end
end
@named resistor = Resistor()
resistor_edge = EdgeModel(resistor, [:src₊v], [:dst₊v], [:src₊i], [:dst₊i]; src=:vs, dst=:cap)
```

Having all those components defined, we can build the network. We don't need to provide a graph
object here because we specified the placement in the graph on a per component basis.
 
```@example mtk
nw = Network([vs_vertex, cap_vertex], [resistor_edge])
```

We can see, that NetworkDynamics internally is able to reduce all of the "output" states. We end up with a plain ODE of a single state.

Now we can simulate the system. For that we generate the `u0` object. Since the metadata (such as default values) was automatically transferred, we can straight away construct the `ODEProblem`
and solve the system.
 
```@example mtk
u0 = NWState(nw) # generate state based on default values
prob = ODEProblem(nw, uflat(u0), (0, 10.0), pflat(u0))
sol = solve(prob, Tsit5())

# plot the solution
fig, ax1, p = plot(sol; idxs=VIndex(:cap, :p₊v), label="Capacitor Voltage");
axislegend(ax1)
ax2 = Axis(fig[2,1])
plot!(ax2, sol; idxs=VIndex(:vs, :i), label="Current produced by ideal v source")
axislegend(ax2)
fig # hide
```

