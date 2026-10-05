# [Callbacks and Events](@id Callbacks)

Callback-functions are a way of handling discontinuities in differential equations.
In a nutshell, the solver checks for some "condition" (i.e. a zero crossing of some variable)
and calls some "affect" if the condition is fulfilled.
Within the affect function, it is safe to modify the integrator, e.g. changing some state or some parameter.

Since `NetworkDynamics.jl` provides nothing more than a RHS for DifferentialEquations.jl, please check
[their docs on event handling](https://docs.sciml.ai/DiffEqDocs/stable/features/callback_functions/)
as a general reference.
This page is introducing the general concepts, for a hands on example of a simulation with callbacks
refer to the [Cascading Failure](@ref) example.

!!! warning
    The `ODEProblem` contains a reference to exactly one copy of the *flat parameter array*.
    If you use callbacks to change those parameters (as we often do), it is advised to
    `copy` the parameter array before passing it to the ODEProblem!
    Also, this means you need to be careful when using the same `prob` for multiple subsequent
    `solve` calls, as the initial state of the `prob` object might have changed!

## Component-based Callback functions
In practice, events often act locally, meaning they only depend and act on a
specific component or type of component. `NetworkDynamics.jl` provides a way of
defining those callbacks on a component level and automatically combine them into performant
[`VectorContinuousCallback`](@extref SciMLBase.VectorContinuousCallback) and [`DiscreteCallback`](@extref SciMLBase.DiscreteCallback) for the whole network.

The main entry points are the types [`ContinuousComponentCallback`](@ref),
[`VectorContinuousComponentCallback`](@ref) and [`DiscreteComponentCallback`](@ref). All of those objects combine a [`ComponentCondition`](@ref) with an [`ComponentAffect`](@ref).

The "normal" [`ContinuousComponentCallback`](@ref) and [`DiscreteComponentCallback`](@ref) have a condition which returns a single value. The corresponding affect is triggered when the return value hits zero.
In contrast, the "vector" version has an in-place condition which writes `len` outputs. When any of those outputs hits zero, the affect is triggered once with an additional argument `event_signs`: a length-`len` vector where entry `i` is `0` (no crossing), `+1` (upcrossing) or `-1` (downcrossing) for the `i`-th output. The affect resolves the crossing direction (and any simultaneous crossings) itself. This mirrors the `affect!(integrator, event_signs)` interface of the underlying [`VectorContinuousCallback`](@extref SciMLBase.VectorContinuousCallback).

There is a special type [`PresetTimeComponentCallback`](@ref) which has no explicit condition and triggers the affect at given times.
This internally generates a [`PresetTimeCallback`](@extref DiffEqCallbacks.PresetTimeCallback) object from `DiffEqCallbacks.jl`.


### Defining the Callback
To construct a condition function, you need to tell NetworkDynamics which symbols of the component you'd like to "observe" within the condition. Any named symbol of the component model works: states, parameters, inputs, outputs and observed. Within the actual condition, those values are made available through `u`:
```julia
condition = ComponentCondition([:x, :y, :p1]) do u, t
    u[:x]  == u[1] # access a state or observable :x at current time
    u[:p1] == u[3] # parameters are listed like any other symbol
    return some_condition(u[:x], u[:y], u[:p1])
end
```
In case of a `VectorContinuousComponentCallback`, the function signature looks slightly different:
```julia
vectorcondition = ComponentCondition([:x, :y, :p1]) do out, u, t
    out[1] = some_condition(u[...])
    out[2] = some_condition(u[...])
    return nothing
end
```
The argument `u` will be passed as a [`SymbolicView`](@ref) object, which means
it is possible to use the getindex syntax to access the desired values by name.

The affect takes the same kind of symbol list:
```julia
affect = ComponentAffect([:u, :p, :obs]) do u, ctx
   t = ctx.t          # extract data from context
   u[:u] = u[:obs]    # states and parameters are writable, observed are readable
   u[:p] = 0
   println("Trigger affect at t=$t")
end
vectoraffect = ComponentAffect([:u, :p]) do u, event_signs, ctx
    for i in eachindex(event_signs)
        event_signs[i] == 0 && continue # skip outputs that did not cross
        if i == 1
            u[:u] = 0 # first output crossed: change state
        else
            u[:p] = 0 # second output crossed: change parameter
        end
        # event_signs[i] is +1 for an upcrossing and -1 for a downcrossing
    end
end
```
Entries of `u` which are states or parameters of the component can be written to, all other
entries (inputs, outputs, observed) are read only. The values in `u` are a snapshot taken when
the affect fires: writing a state or parameter updates the integrator immediately, but an
observed entry which depends on it is not recomputed within the same affect.

The affect gets passed a `ctx` "context" object, which is a named tuple which holds additional context like the integrator object, the component model, the index of the component model, the current time and so on. Please refer to the [`ComponentAffect`](@ref) docstring for a detailed list.

!!! note "Legacy form with separate parameter list"
    Earlier versions took two symbol lists, `ComponentCondition(f, sym, psym)` with `f(u, p, t)`
    and `ComponentAffect(f, sym, psym)` with `f(u, p, ctx)`, where `p` gave access to the
    parameters. This form is still accepted and behaves as before.

Lastly we need to define the actual callback object using [`ContinuousComponentCallback`](@ref)/[`VectorContinuousComponentCallback`](@ref):
```julia
ccb  = ContinuousComponentCallback(condition, affect; kwargs...)
vccb = VectorContinuousComponentCallback(condition, affect, len; kwargs...)
```
where the `kwargs` are passed to the underlying [`SciMLBase.VectorContinuousCallback`](@extref) to finetune the zerocrossing-detection.


### Registering the Callback
Once the callback is defined, we need to "attach" it to the component, for that you can use the methods [`add_callback!`](@ref) and [`set_callback!`](@ref):
```julia
vert = VertexModel(...)
add_callback!(vert, ccb)
add_callback!(vert, vccb)
```


### Extracting the Callback
Component-level callbacks are automatically extracted and combined when constructing an [`ODEProblem`](@ref SciMLBase.ODEProblem(::NetworkDynamics.Network, ::NetworkDynamics.NWState, ::Any)):
```julia
u0 = NWState(nw)
prob = ODEProblem(nw, u0, (0, 10))
sol = solve(prob, ...)
```

For more control over callback handling—such as adding network/system-level callbacks (e.g., `PeriodicCallback`),
temporary component callbacks, or overriding the default callbacks, the `ODEProblem` constructor provides
keyword arguments `add_comp_cb`, `add_nw_cb`, and `override_cb`. See the [`ODEProblem(nw::Network,...)`](@ref SciMLBase.ODEProblem(::NetworkDynamics.Network, ::NetworkDynamics.NWState, ::Any)) documentation for details.

When executing component callbacks, NetworkDynamics automatically checks whether states or parameters
changed during the affect and calls [`SciMLBase.auto_dt_reset!`](@extref) and [`save_parameters!`](@ref) if necessary.
An affect which only does bookkeeping, for example counting events in a parameter, can keep the
current step size by setting `ctx.dt_reset[] = false`. The parameter change is saved either way.


## Event Iteration
Discrete logic such as hysteresis switches, timers or sample-and-hold blocks keeps its memory in
parameters which callbacks set. Such blocks are meant to act synchronously. At one event instant
they should all see the same state, and one block's switch may cause another block to switch at
the same instant. Plain discrete callbacks don't work like that: separate callbacks run one after
another, so the result depends on their order, and a cascade within one instant is only noticed
in the next step.

For this, a discrete callback can join the event iteration of the network:
```julia
cb = DiscreteComponentCallback(condition, affect; iterative=true)
```
All iterative callbacks of the network run together as one set, after every other callback at an
event instant. They run in rounds:

1. All conditions are evaluated on the same state.
2. The affects of all fired callbacks run. They all read one snapshot, taken before the first of
   them writes.
3. If the network is a DAE, its algebraic states are reinitialized with the `initializealg` of
   the solve. Then the next round starts at step 1.

The loop ends as soon as no condition fires anymore. Step size reset and parameter saving happen
once for the whole instant.

Because the conditions are evaluated several times within one instant, they should be predicates
on the current state alone, for example "the stored switch state contradicts the input". The
affect which fixes the inconsistency then also makes its own condition false.
Iterative affects should read and write only through `u`. Writes through `ctx.integrator` are
invisible to the iteration, so they trigger no reinit and no further round.

The order of the callbacks at an event instant is:
1. the continuous callback which found the earliest root, if any,
2. the preset-time callbacks,
3. the other discrete callbacks,
4. the iterative set.

The callbacks passed as `add_nw_cb` to the `ODEProblem` are sorted into the same groups. Within a
group, the network's own callbacks come first.

So the event iteration reacts to the jumps of all other callbacks. This holds when the callbacks
come from the `ODEProblem` constructor. If you combine `get_callbacks(nw)` with your own callbacks
in a `CallbackSet` by hand, your discrete callbacks run after the iterative set.
If the set still fires after `event_maxiter` rounds (default 10), it keeps the state of the last
round and warns, naming the components which still fire. Pass `event_failure=:error` to throw
instead. Both are keyword arguments of [`get_callbacks`](@ref) and of the `ODEProblem`
constructor.


## Normal DiffEq Callbacks
Besides component based callbacks, it is also possible to use "normal" DiffEq
callbacks together with `NetworkDynamics.jl`.
It is far more powerful but also more cumbersome compared to the component based callback functions.
To access states and parameters of specific components, we heavily rely on the [Symbolic Indexing](@ref) features.

```julia
using SymbolicIndexingInterface as SII
nw = Network(#= some network =#)

condition = let getvalue = SII.getsym(nw, VIndex(1:5, :some_state))
    function(out, u, t, integrator)
        s = NWState(integrator, u, integrator.p, t)
        some_state = getvalue(s)
        out .= some_condition(some_state)
    end
end
```
Please note a few important things here:
 - Symbolic indexing can be costly, and the condition function gets called very
   often. By using [`SII.getsym`](@extref `SymbolicIndexingInterface.getsym`) we did
   some of the work *before* the callback by creating the accessor function.
   When handling with "normal states" and parameters consider using
   [`SII.variable_index`](@extref `SymbolicIndexingInterface.variable_index`) and
   [`SII.parameter_index`](@extref `SymbolicIndexingInterface.parameter_index`) for
   even better access patterns.
 - `t` refers to the current time of the zerocrossing-detection-algorithm. This is different from `integrator.t` which refers to the current timestep in which the zerocross-detectio takes place..

```julia
function affect!(integrator, vidx)
    p = NWParameter(integrator) # get symbolically indexable parameter object
    p.v[vidx, :some_vertex_parameter] = 0 # change some parameter
    auto_dt_reset!(integrator)
    save_parameters!(integrator)
end
```
The affect function is much more straight forward, as it (typically) is called far less frequent and thus less perfomance critical.

Once the `condition` and `affect!` is defined, you can use the [`SciMLBase.ContinuousCallback`](@extref) and [`SciMLBase.VectorContinuousCallback`](@extref) constructors to create the callback.

!!! note "Introducing discontinuities with adaptive timestepping"
    Since changes to `u` and `p` mostly introduce discontinuities in the
    solution, it is recommend to call [`auto_dt_reset!`](@extref
    `SciMLBase.auto_dt_reset!`) within the affect to restart integration with
    small steps afterwards.

!!! note "Changing Parameters and Observables"
    An "observable" is kind of a "virtual" state, which can be reconstructed for
    a given time `t`, a given state `u` and a given set of parameters `p`
    ```math
    o = f(u(t), p(t), t)
    ```
    To extract or plot timeseries of observed states under *time variant
    parameters* (i.e. parameters that are changed in a callback), those changes
    need to be recorded using the [`save_parameters!`](@ref) function whenever `p` is changed.
    When using [ComponentCallback](@ref NetworkDynamics.ComponentCallback), NetworkDynamics will automatically check for changes in `p` and save them if necessary.
