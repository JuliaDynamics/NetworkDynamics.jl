# Initialization
Initialization of the system describes the process of finding valid initial conditions, primarily a fixed point of the system.
We distinguish between two types of initialization: full system initialization and component initialization.

## Full-System Initialization
Full system initialization describes the process of finding a fixed point/steady state of the entire system.

To do so, you can use [`find_fixpoint`](@ref), which creates a `SteadyStateProblem` of the whole network and attempts to solve it.

## Component-wise Initialization
In contrast to full-system initialization, the goal of component-wise initialization is to find a valid initial condition for a single component first, given a network coupling.

This can be useful in cases where there are nontrivial internal dynamics and states within a single vertex or edge.
The idea of component-wise initialization is to find internal states that match a given "network coupling" (fixed inputs and outputs).

### Mathematical Meaning
According to the [Mathematical Model](@ref) of NetworkDynamics.jl, a component forms an "input-output-system" of the form

```math
\begin{aligned}
M\,\frac{\mathrm{d}}{\mathrm{d}t}x &= f(x, i, p, t)\\
y &= g(x, i, p, t)
\end{aligned}
```
where $x$ are the internal states, $i$ are the inputs, $y$ are the outputs, and $p$ are the parameters.
To initialize at a fixed point, we require the RHS to be zero,
```math
\begin{aligned}
0 &= f(x, i, p, t)\\
0 &= g(x, i, p, t) - y
\end{aligned}
```
forming a nonlinear least squares problem for the residual.
Each variable in $x$, $i$, $y$ and $p$ is either considered **free** or **fixed** with respect to the nonlinear problem stated above.
Symbols that have a **default** value (see [Metadata](@ref)) are considered *fixed*.
All other symbols are considered *free* and must provide a **guess** value as an initial starting point for the nonlinear solver.

The **defaults** and **guesses** can be either obtained from the [Metadata](@ref) directly or provided as arguments.


### Typical Workflow
The following initialization workflow (as used in the [Tutorial on Initialization](@ref init-tutorial)) is quite common for complex dynamics:

  1. Define simple, quasi-static models with the same input-output structure as the dynamic models.
  2. Find a solution of the static model using [`find_fixpoint`](@ref).
  3. Define elaborate, dynamical models. Define dynamical network on the same graph topology.
  4. Use [`set_interface_defaults!`](@ref) to "copy" the inputs and outputs of all components from the static solution to the dynamical network.
  5. Use [`initialize_component!`](@ref) to determine "free" states and parameters within the dynamical models to reach a steady state while satisfying the constraints on inputs/outputs.


### Non-mutating vs Mutating Initialization

NetworkDynamics provides two approaches for component-wise initialization:

1. **Non-mutating approach** using [`initialize_component`](@ref): Returns a dictionary of values without modifying the component.
2. **Mutating approach** using [`initialize_component!`](@ref): Directly updates the component metadata with initialization results.

Both options take guesses and defaults from metadata by default; however, it is possible to specify otherwise (see method documentation).

The non-mutating version [`initialize_component`](@ref) is useful when you don't want to modify the metadata, for a more "stateless" approach:
```julia
# Get initialization results as a dictionary
init_state = initialize_component(vf; default_overrides=Dict(:x => 4))
```
It will return a `Dict{Symbol,Float64}` which contains values for **all** symbols in the model.


The mutating version [`initialize_component!`](@ref) directly updates the component's metadata with initialization results:

```julia
initialize_component!(vf; verbose=true) # set `init` metadata for free symbols
```

### Example
Let's consider the following example of a Swing-equation generator model.
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
        V, [guess=1.0, description="Voltage magnitude"]
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
You can see in the provided [metadata](@ref) that we've set `default` values for the node outputs `u_r`, `u_i`, the node inputs `i_r`, `i_i`, and most parameters.
For some states and parameters, we've only provided a `guess` rather than a default.
Variables that only have `guess`es are considered "tunable" for the initialization algorithm.

### Using the non-mutating initialization

We can use [`initialize_component`](@ref) to get the initialized values without modifying the component:
```@example compinit
init_values = initialize_component(vf; default_overrides=Dict(:u_i=>0), verbose=true)
```
The code returns a dictionary which pins *all* the variables of the component to some values which satisfy the initialization condition.

### Using the mutating initialization

Alternatively, we can make use of the mutating version to store the results of the initialization in the metadata:
```@example compinit
initialize_component!(vf; verbose=true)
```
Which stored the initialisation results as `:init` metadata of the component model `vf`:
```@example compinit
get_init(vf, :V) # get the value of :V at initialized state
```

It is possible to inspect initial states (works also for observed symbols) using [`get_initial_state`](@ref).
As a quick test we can ensure that the angle indeed matches the voltage angle:
```@example compinit
get_initial_state(vf, :θ) ≈ atan(get_initial_state(vf, :u_i), get_initial_state(vf, :u_r))
```
You can print out the whole state using [`dump_initial_state`](@ref).
```@example compinit
dump_initial_state(vf)
```

### Additional Initialization Constraints
Sometimes it is required to add additional initialization constraints for components. This is done by using [`set_initconstraint!`](@ref) in combination with an [`InitConstraint`](@ref) object.

In a nutshell, an additional initialization constraint is a function of states/inputs/... whose residual should be zero at the initial state:
```math
0 = f_\mathrm{additional}(x, i, p, t)
```
Such a function is constructed using [`InitConstraint`](@ref) constructor:
```@example compinit
additional_constraint = InitConstraint([:Pel, :u_r, :u_i], 2) do out, u
    out[1] = u[:Pel] - 1.0
    out[2] = sqrt(u[:u_r] ^ 2 + u[:u_i] ^ 2) - 1
end
nothing #hide
```
We need to pass a list of symbols we want to access, the number of additional equations (2 in this case), and a function which modifies the residual in-place.
You can access inputs, states, outputs, and observables within the constraint function.

For convenience, there is the [`@initconstraint`](@ref) macro to generate such constraints with less overhead. The definition below is equivalent to the definition above.
```@example compinit
additional_constraint = @initconstraint begin
    :Pel - 1.0
    sqrt(:u_r^2 + :u_i^2) - 1
end
nothing #hide
```
Once we attach the additional constraint to the component model, it is also indicated in the printout:
```@example compinit
set_initconstraint!(vf, additional_constraint)
vf
```
