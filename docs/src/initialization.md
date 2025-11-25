# [Initialization](@id initialization-guide)

Initialization is a critical step in simulation dynamical systems on networks, involving finding valid initial conditions that satisfy the system's constraints. NetworkDynamics provides several layers of initialization tools, from individual component initialization to full network initialization.

## Initialization Hierarchy

NetworkDynamics offers a tiered approach to initialization:

1. **Full-System Initialization**: Finding a steady state for the entire network at once
2. **Component-wise Network Initialization**: Initializing each component individually while respecting network coupling
3. **Single Component Initialization**: Finding valid internal states for a single component

## Full-System Initialization

Full system initialization aims to find a fixed point/steady state of the entire system simultaneously.

To do so, you can use [`find_fixpoint`](@ref), which creates a `SteadyStateProblem` of the whole network and attempts to solve it:

```julia
# Create a network
nw = Network(vertices, edges)

# Find a fixed point for the entire system
state = find_fixpoint(nw)
```

This approach works well for simpler systems but may face convergence challenges for complex networks with many interacting components.

## Component-wise Network Initialization

For more complex networks, a component-by-component approach is often more robust. NetworkDynamics provides [`initialize_componentwise`](@ref) and [`initialize_componentwise!`](@ref) functions that:

1. Initialize each component individually
2. Verify the combined solution works for the entire network

```julia
# Initialize each component in the network individually
state = initialize_componentwise(nw)

# Or using the mutating version that updates component metadata
state = initialize_componentwise!(nw)
```

### Two-Step Initialization Pattern

A common initialization pattern for complex networks involves:

1. Solving a simplified static model first
2. Using those results to initialize a more complex dynamic model

```julia
# 1. Solve a static/simplified model
static_model = create_static_network(...)
static_state = find_fixpoint(static_model)

# 2. Extract interface values
interface_vals = interface_values(static_state)

# 3. Use them to initialize a dynamic model
dynamic_model = create_dynamic_network(...)
dyn_state = initialize_componentwise(dynamic_model, default_overrides=interface_vals)
```

See the [Tutorial on Initialization](@ref init-tutorial) for a complete example of this approach.

## Single Component Initialization

At the lowest level, NetworkDynamics provides tools for initializing individual components based on their internal dynamics and interface constraints.

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

### Non-mutating vs Mutating Initialization

NetworkDynamics provides two approaches for component-wise initialization:

1. **Non-mutating approach** using [`initialize_component`](@ref): Returns a dictionary of values without modifying the component
2. **Mutating approach** using [`initialize_component!`](@ref): Directly updates the component metadata with initialization results

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

The same pattern applies at the network level with [`initialize_componentwise`](@ref) and [`initialize_componentwise!`](@ref).

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

### Advanced Component Initialization: Formulas and Constraints

NetworkDynamics provides three complementary mechanisms for customizing component initialization beyond the basic defaults and guesses: **initialization formulas** ,**initialization constraints** and **guess formulas**. These operate at different stages of the initialization pipeline and serve distinct purposes:
- init formulas **reduce the number of free variables** (by setting additional defaults),
- init constraints **increase the number of equations** for the init problem and
- guess formulas **refine starting values** for the initialization rootfinding problem.

The execution order is:
1. Collect defaults, guesses, and bounds from metadata
2. Apply init formulas → update defaults (fix more variables)
3. Apply guess formulas → update guesses (improve starting point)
4. Create and solve the nonlinear least squares problem with additonal constraints

#### Initialization Formulas (`InitFormula`)

Initialization formulas act early in the initialization pipeline to compute and set default values based on other known values. They are particularly useful for deriving dependent quantities or ensuring consistency between related variables.

Each formula can only reference symbols that are already available - it cannot use intermediate values computed within the same formula.

**Basic Usage**: Use the [`@initformula`](@ref) macro to define formulas with assignment syntax:

```@example compinit
# Example: Set voltage magnitude and electrical power based on voltage components
voltage_formula = @initformula begin
    :V = sqrt(:u_r^2 + :u_i^2)     # Voltage magnitude from components
    :Pel = :u_r * :i_r + :u_i * :i_i # Electrical power calculation
end
nothing #hide
```

**Applying Formulas**: Formulas can be either added to the metadata of components ([`set_initformula!`](@ref), [`add_initformula!`](@ref)) or passed as `additional_initformula` to the
[`initialize_component[!]`](@ref NetworkDynamics.initialize_component) functions.

**Dependency Resolution**: When applying multiple separate formulas, NetworkDynamics automatically sorts them topologically to ensure correct evaluation order.


#### Initialization Constraints (`InitConstraint`)

Initialization constraints add equations to the nonlinear system that must be satisfied during the initialization solve. Unlike formulas, they don't directly set values but impose mathematical relationships.

**Basic Usage**: Use the [`@initconstraint`](@ref) macro to define constraint equations:

```@example compinit
# Example: Constrain electrical power and voltage magnitude
power_constraint = @initconstraint begin
    :Pel - 1.0                        # Electrical power must equal 1.0
    sqrt(:u_r^2 + :u_i^2) - 1.0      # Voltage magnitude must equal 1.0
end
nothing #hide
```

**Applying Constraints**: Constraints can be either added to the metadata of components ([`set_initconstraint!`](@ref), [`add_initconstraint!`](@ref)) or passed as `additional_initconstraint` to the
[`initialize_component[!]`](@ref NetworkDynamics.initialize_component) functions.


#### Guess Formulas (`GuessFormula`)

Guess formulas operate later in the initialization pipeline than init formulas. While init formulas set default values (thereby reducing the number of free variables), guess formulas refine the initial guesses for free variables to improve solver convergence without changing the problem dimension.

While InitFormulas use defaults to update other defaults, GuessFormulas update guesses based on other defaults and guesses (if some variable has a default and guess defined, the default takes precedence).
They act **after** the InitFormulas, thus having access to the updated defaults.

Similar to InitFormulas, NetworkDyanmics makes sure that there are no circular dependencies between guess formulas and you can't have multiple guess formulas updating the same variables.
It performs some topological sorting on multiple guesses to update in order.

**Basic Usage**: Use the [`@guessformula`](@ref) macro with the same assignment syntax as init formulas:

```@example compinit
# Example: Improve angle and voltage guesses from interface values
setpoint_guesses = @guessformula begin
    :V_set = sqrt(:u_r^2 + :u_i^2)     # guess voltage mag setpoint close to actual voltage
    :P_set = :u_r * :i_r + :u_i * :i_i # guess power setpoint close to actual power
end
nothing #hide
```

**Applying GuessFormulas**: Like init formulas and constraints, guess formulas can be stored in metadata ([`set_guessformula!`](@ref)) or passed directly using the `additional_guessformula` keyword in [`initialize_componentwise`](@ref), [`initialize_component`](@ref) and friends.


## Analysing Fixpoints
In order to analyse fixpoints NetworkDynamis provides the functions [`isfixpoint`](@ref), [`is_linear_stable`](@ref) and [`jacobian_eigenvals`](@ref).
