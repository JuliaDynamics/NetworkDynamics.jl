# Initialization
Initialization of the system describes the process of finding valid initial conditions, mostly a fixpoint of the system.
We distinguish between two types of initialization: full system initialziation and component initialization.

## Full-System Initialization
Full system initialization describs the process of finding a fixpoint/steady state of th entire system.

To do so, you can use [`find_fixpoint`](@ref), which creates a `SteadyStateProblem` of the whole network and tries do solve it. 

## Component-wise Initialization
In contrast to full-system initialization the goal of component-wise initialization is to find a valid initial condition for a single component first, given a network coupling.

This can be usefull in cases, where there are nontrivial internal dynamics and states within a single vertex or edge.
The idea of component-wise initialisation is to find internal states which match a given "network coupling" (fixed inputs and outputs).

Lets consider the following example of a Swing-equation generator model.
```@example compinit
using NetworkDynamics, ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as Dt

@mtkmodel Swing begin
    @variables begin
        u_r(t)=1, [description="bus d-voltage", output=true]
        u_i(t)=0.1, [description="bus q-voltage", output=true]
        i_r(t)=1, [description="bus d-current (flowing into bus)", input=true]
        i_i(t)=0.1, [description="bus d-current (flowing into bus)", input=true]
        ω(t), [guess=0.0, description="Rotor frequency"]
        θ(t), [guess=0.0, description="Rotor angle"]
        Pel(t), [guess=1, description="Electrical Power injected into the grid"]
    end
    @parameters begin
        M=0.005, [description="Inertia"]
        D=0.1, [description="Damping"]
        V=sqrt(u_r^2 + u_i^2), [description="Voltage magnitude"]
        ω_ref=0, [description="Reference frequency"]
        Pm, [guess=0.1,description="Mechanical Power"]
    end
    @equations begin
        Dt(θ) ~ ω - ω_ref
        Dt(ω) ~ 1/M * (Pm - D*ω - Pel)
        Pel ~ u_r*i_r + u_i*i_i
        u_r ~ V*cos(θ)
        u_i ~ V*sin(θ)
    end
end
sys = Swing(name=:swing)
vf = VertexModel(sys, [:i_r, :i_i], [:u_r, :u_i])
```
You can see in the provided [metadata](@ref), that we've set `default` values for the node outputs `u_r`, `u_i`, the node inputs `i_r`, `i_i` and most parameters.
For some states and parameters, we've onlye provided a `guess` rather than a default.
Variables which only have `guess`es are considered "tunable" for the initialization algorithm.

In order to initialize the remaining variables we use [`initialize_component!`](@ref), which is a mutating function which tries to solve the nonlinear initialization problem and store the found values for the "free" variables as `init` metadata.

```@example compinit
initialize_component!(vf; verbose=true)
nothing #hide
```
```@example compinit
vf #hide
```

Which lead to a successfull initialization of states `:θ` and `:ω` as well as parameter `:Pm`.
To retrieve the residual you can use [`init_residual`](@ref).

As a quick test we can ensure that the angle indeed matches the voltag angel:
```@example compinit
get_init(vf, :θ) ≈ atan(get_default(vf, :u_i), get_default(vf, :u_r))
```

It is possible to inspect initial states (also for observed symbols) using [`get_initial_state`](@ref). You can print out the whole state using [`dump_initial_state`](@ref).
```@example compinit
dump_initial_state(vf)
```
