# ModelingToolkit Integration

NetworkDynamics.jl is compatible with [`ModelingToolkit.jl`](https://github.com/SciML/ModelingToolkit.jl) (MTK).
The general idea is to use MTK to define *component models* (i.e. edge and vertex dynamics)
which are then connected on network scale using NetworkDynamics.

The main entry point for this interop are the constructors
```julia
VertexModel(::System, inputs, outputs)
EdgeModel(::System, srcin, dstin, [srscout], dstout)
```
whose docstrings can be found in the [Component Models with MTK](@ref) section in the API.

These constructors will:
- transform the states marked as input to parameters and `mtkcompile`ing the system,
- generate the `f` and `g` functions,
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
The `VertexModel` can be generated from an `System` by providing names of the input and output states:

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

## Fully Implicit Outputs
When working with MTK systems in NetworkDynamics, you may encounter situations where
your desired output variables don't explicitly appear in the equations. This creates **fully
implicit outputs** - variables that are determined by the system's constraints but aren't
directly computed.

!!! tip "tl;dr"
    Introduce "fake" dependencies to your input-forcing equations `0 ~ in + implicit_output(y)`.
    Which is mathematically equivalent to `0 ~ in` but helps MTK to reason about dependencies.

Consider a system with a fully implicit output:
```
   u┌───────┐y
  ─→┤ 0 ~ u ├→─
    └───────┘
```
Here, $y$ does not appear in the equations at all. In general, that doesn't make too much sense.
During simplification, MTK will potentially get rid of the equation as it does not contribute to the system's state.

However, in NetworkDynamics, we're always dealing with **open loop models** on the equation level, which is not exactly what MTK was made for.
If you build a closed loop between a subsystem A which **has input forcing** and a subsystem
B which has **input feed forward**, the resulting system can be solved:
```
    (system with input forcing)
          ua┌─────────┐ya
        ┌──→┤  0 ~ ua ├→──┐
        │   └─────────┘   │
        │ yb┌─────────┐ub │
        └──←┤ yb ~ ub ├←──┘
            └─────────┘
(system with input feed forward)
```

Since MTK does not know about the closed loop (which is only introduced on the NetworkDynamics level once we leave the equation based domain) we need to help MTK to figure out those dependencies.
We can do so by introducing "fake" dependencies using [`implicit_output`](@ref).
This function is defined as
```julia
implicit_output(x) = 0
ModelingToolkit.@register_symbolic implicit_output(x)
```
which makes it numerically equivalent to zero (no effect on the simulation) but is
opaque to the Symbolic Simplification.

### Example

Consider a "Kirchhoff Node" between multiple resistors:
- the currents through the resistors directly depend on the voltage output of the node (input feed forward) and
- the Kirchhoff node requires the sum of all inflowing currents to be zero (input forcing).

We can model this type of node like this:
```@example mtk
@mtkmodel KirchhoffNode begin
    @variables begin
        v(t), [description="Node voltage", output=true]
        i_sum(t), [description="Sum of incoming currents", input=true]
    end
    @equations begin
        0 ~ i_sum + implicit_output(v)  # Kirchhoff's current law
    end
end
@named kirchhoff = KirchhoffNode()
VertexModel(kirchhoff, [:i_sum], [:v])
```
where we "trick" MTK into believing that the input forcing equation depends on the output too.

## Register Postprocessing Functions
MTK models are "composite" by design, i.e. your toplevel model might consist of multiple internal components.
Sometimes, you want to specify behavior of the **encapsulating model** on a subcomponent level. 
For that, NetworkDynamics provides the [`ComponentPostprocessing`](@ref) metadata mechanism.

Let's say you want to implement an integrator with anti windup
```asciiart
              __ outMax
             /
        ╭─────╮
     in │  K  │ out
    ╶───┤╶───╴├────╴
        │ s T │
        ╰─────╯
outMin __/
```

We could try to implement this using `ifelse`  functions, but the proper way is to do it with callbacks.
Since NetworkDynamics cannot understand events defined within MTK models, we need to attach the callbacks on the component level.
That's annoying, because the building block might appear in lots of different models.
To work around this shortcoming, ND.jl allows you to define postprocessing functions, which will be called on every model which makes use of the subsystem.

In our model, we add two "internal" parameters, to indicate if the system is in lower or upper saturation.
Depending on the saturation state, the integrator is "disabled".
The model definition looks like this. The important part is the `@metadata` block at the end.

```@example mtk
function attach_limint_callback! end # function needs to exist before model
@mtkmodel LimitedIntegrator begin
    @structural_parameters begin
        K # Gain
        T # Time constant
        outMin # Lower limit
        outMax # Upper limit
        guess=0
    end
    @parameters begin
        _callback_sat_max = 0
        _callback_sat_min = 0
    end
    @variables begin
        in(t), [description="Input signal", input=true]
        out(t), [guess=guess, description="limited integrator output state", output=true]
        min(t), [description="Lower limit"]
        max(t), [description="Upper limit"]
        forcing(t)
    end
    @equations begin
        min ~ outMin
        max ~ outMax
        forcing ~ K*in
        T*D(out) ~ (1 - _callback_sat_max - _callback_sat_min) * forcing
    end
    @metadata begin
        ComponentPostprocessing = attach_limint_callback!
    end
end
nothing #hide
```

The callback to generate is split in three separate conditions:
- if `out - max` registers an upcrossing, enter upper saturation
- if `-out + min` registers an upcrossing, enter lower saturation
- if currently in saturation, check for zero crossings in `forcing`, if in upper saturation and
forcing turns negative, unsaturate. If in lower saturation and forcing turns positive unsaturate.

The function below generates such a callback for the component at a given namespace:

```@example mtk
function attach_limint_callback!(cf, namespace)
    # generate the required symbols based on the namespace
    min = Symbol(namespace, "₊min")
    max = Symbol(namespace, "₊max")
    out = Symbol(namespace, "₊out")
    forcing = Symbol(namespace, "₊forcing")
    satmax = Symbol(namespace, "₊_callback_sat_max")
    satmin = Symbol(namespace, "₊_callback_sat_min")

    condition = ComponentCondition([min, max, out, forcing], [satmax, satmin]) do _out, u, p, _
        insatmax = !iszero(p[satmax])
        insatmin = !iszero(p[satmin])

        upcrossing_max =  u[out] - u[max]
        upcrossing_min = -u[out] + u[min]

        # enable upper saturation
        _out[1] = insatmax ? Inf : upcrossing_max
        # enable lower saturation
        _out[2] = insatmin ? Inf : upcrossing_min
        if insatmax || insatmin
            # when in saturation, check for zero crossing of forcing
            _out[3] = u[forcing]
        else
            _out[3] = Inf
        end
    end

    upcrossing_affect = ComponentAffect([], [satmax, satmin]) do u, p, eventidx, ctx
        if eventidx == 1
            println("$namespace: /⎺ reached upper saturation at $(round(ctx.t, digits=4))s")
            p[satmax] = 1.0
            p[satmin] = 0.0
        elseif eventidx == 2
            println("$namespace: \\_ reached lower saturation at $(round(ctx.t, digits=4))s")
            p[satmax] = 0.0
            p[satmin] = 1.0
        elseif eventidx == 3
            # upcrossing means, forcing went from negative to positive, i.e. we leave lower saturation
            insatmin = !iszero(p[satmin])
            if insatmin
                println("$namespace: _/ left lower saturation at $(round(ctx.t, digits=4))s")
                p[satmin] = 0.0
            end
        else
            error("Unknown event index $eventidx")
        end
    end

    downcrossing_affect = ComponentAffect([],[satmax]) do u, p, eventidx, ctx
        if eventidx == 3 # downcrossing means nothing for saturation affects
            # downcrossing means, forcing went from positive to negative, i.e. we leave upper saturation
            insatmax = !iszero(p[satmax])
            if insatmax
                println("$namespace: ⎺\\ left upper saturation at $(round(ctx.t, digits=4))s")
                p[satmax] = 0.0
            end
        else
            error("Unknown event index $eventidx")
        end
    end

    cb = VectorContinuousComponentCallback(condition, upcrossing_affect, 3; affect_neg! = downcrossing_affect)

    # finally add callback to component
    NetworkDynamics.add_callback!(cf, cb)
end
nothing #hide
```
!!! warning
    Even though fairly complex, the above function is a simplified example which may have performance problems 
    and does not add "safety" discrete callbacks. Under discrete jumps, the zero crossing might be skipped so 
    it is good practice to add additional discrete callbacks to detect those cases (i.e. check if `out >= max + 1e-10`).
    You'll find a "proper" version of the saturated integrator in the PowerDynamics.jl Library, feel free 
    to copy the code from there.

With the definition above, we can create a model which uses the limited integrator:
```@example mtk
@mtkmodel ComponentWithLimInt begin
    @components begin
        int = LimitedIntegrator(K=1, T=1, outMin=-1, outMax=1)
    end
    @variables begin
        out(t)
        in(t)
    end
    @equations begin
        int.in ~ in
        out ~ int.out
    end
end

@named testcomp = ComponentWithLimInt()
nothing #hide
```
When we build this (useless) vertex model
```@example mtk
vm = VertexModel(testcomp, [:in], [:out])
```
we see the callback was generated and attached automatically.
