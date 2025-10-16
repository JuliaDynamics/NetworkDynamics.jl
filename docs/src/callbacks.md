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
In contrast, the "vector" version has an in-place condition which writes `len` outputs. When any of those outputs hits zero, the affect is triggered with an additional argument `event_idx` which tells the effect which dimension encountered the zerocrossing.

There is a special type [`PresetTimeComponentCallback`](@ref) which has no explicit condition and triggers the affect at given times.
This internally generates a [`PresetTimeCallback`](@extref DiffEqCallbacks.PresetTimeCallback) object from `DiffEqCallbacks.jl`.


### Defining the Callback
To construct a condition function, you need to tell network dynamics which states and parameters you'd like to "observe" within the condition. Within the actual condition, those states will be made available:
```julia
condition = ComponentCond([:x, :y], [:p1, :p2]) do u, p, t
    u[:x]  == u[1] # access a state or observable :x at current time
    p[:p2] == p[2] # access a parameter at current time
    return some_condition(u[:x], u[:y], ...)
end
```
In case of a `VectorContinuousComponentCallback`, the function signature looks slightly different:
```julia
vectorcondition = ComponentCond([:x, :y], [:p1, :p2]) do out, u, p, t
    out[1] = some_condition(u[...], p[...])
    out[2] = some_condition(u[...], p[...])
    return nothing
end
```
Note that the `syms` argument (here `[:x, :y]`) can be used to reference **any**
named state of the component model, this includes "ordinary" states, observed,
inputs and outputs.
The arguments `u` and `p` will be passed as [`SymbolicView`](@ref) objects, which mean
it is possible to use the getindex syntax to acces the desired states by name.

The affect takes a similar form:
```julia
affect = ComponentAffect([:u], [:p]) do u, p, ctx
   t = ctx.t # extract data from context
   obs = NWState(ctx.integrator)[VIndex(ctx.vidx, :obs)] # extract some observed state from context
   println("Trigger affect at t=$t")
end
vectoraffect = ComponentAffect([:u], [:p]) do u, p, event_idx, ctx
    if event_idx == 1
        u[:u] = 0 # change state
    else
        u[:p] = 0 # change parameter
    end
end
```
Notably, the `syms` (here `:u`) can *exclusively* refer to "ordinary" states, since they are now writable.
However the affect gets passed a `ctx` "context" object, which is a named tuple which holds additional context like the integrator object, the component model, the index of the component model, the current time and so on. Please refer to the [`ComponentAffect`](@ref) docstring for a detailed list.

Lastly we need to define the actual callback object using [`ContinuousComponentCallback`](@ref)/[`VectorContinuousComponentCallback`](@ref):
```julia
ccb  = ContinuousComponentCallback(condition, affect; kwargs...)
vccb = VectorContinuousComponentCallback(condition, affect; kwargs...)
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
