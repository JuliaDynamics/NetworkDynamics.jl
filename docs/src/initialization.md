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

### Aliases: the name you pick doesn't matter

When a component is assembled from ModelingToolkit models, one physical quantity often ends up with several names. A bus voltage might be reachable both as `busbar₊u_r` and as `terminal₊u_r`, because an equation `busbar₊u_r ~ terminal₊u_r` ties the two together. NetworkDynamics recognizes such aliases and treats each group of them as a *single* variable.

For initialization this means you may attach a **default** or **guess** to *any* member of an alias group and it counts for the whole group — you never have to track down the "canonical" name. Giving the same value to two aliases of one quantity is fine; only genuinely *conflicting* values are reported. This happens automatically, before the initialization problem is assembled.

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
using NetworkDynamics
using SciCompDSL
using ModelingToolkitBase
using ModelingToolkitBase: t_nounits as t, D_nounits as Dt

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

A formula takes precedence over a plain default: if a symbol carries both a numeric `default` and an `initf`, the formula overwrites the default. That is by design — it lets a component ship a sensible standalone value which is superseded once the surrounding model pins the symbol down.

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

**From MTK models**: when a component is built from a ModelingToolkit `System`, init formulas are usually not written by hand but declared next to the equations using the `initf` variable option:

```julia
@component function avr(; name)
    @variables begin
        v_mag(t)
        v_meas(t), [guess=1, initf = v_mag]   # measurement equals actual value in steady state
    end
    @parameters begin
        v_ref, [guess=1, initf = v_meas]      # zero control error in steady state
    end
    # ...
end
```

`initf` reads as "at initialization, set this symbol to that expression". It works on unknowns **and** parameters — back-computing a free setpoint parameter from a known operating point is one of the main uses. The expression may reference any state, parameter or observable of the system, and formulas from different subcomponents connect to each other automatically (see [below](@ref backward-flow-init)), so a whole nested model can be initialized from a handful of interface values.

**Attaching a formula from outside (`set_initf`)**: the `initf` option annotates a variable *where it is declared*. Its postfix counterpart is [`set_initf`](@ref), which attaches the same kind of rule to an *already-built* system without touching the declaration:

```julia
sys = set_initf(sys, v_ref => v_meas)     # same effect as the `initf = v_meas` option above
```

Reach for `set_initf` when the rule comes from *outside* the block that owns the symbol — most often a parent component initializing one of its subsystems' symbols, which the `initf` option cannot express because it only annotates a system's *own* variables. `set_initf` returns a new system rather than mutating, so rebind it (`sys = set_initf(sys, …)`).

!!! warning "Symbolic expressions in metadata: be explicit about *when* they hold"
    An expression attached to a variable does not say **when** it is supposed to be true. Is `x = 2a` a statement about the initial state, or a relation that must hold for all time? Those mean very different things, so NetworkDynamics does not guess: a symbolic expression passed to `default` (on an unknown) or to `guess` is **rejected**.

    Say it explicitly instead:

    - `initf = <expr>` holds **once, at initialization**. It is evaluated while the init problem is assembled and is not enforced afterwards — during the simulation the state is free to drift away from it.
    - `guessf = <expr>` has the same timing but is only a *hint*: it seeds the solver's starting point and is never checked.

    **Exception: symbolic defaults on parameters.** `@parameters K = K_e` stays valid, because it means something categorically different — it holds **unconditionally, at all times**. MTK demotes `K` to an observable of `K_e`, so from then on there is only one parameter. This is the mechanism for wiring parameter dependencies between a model and its subcomponents: a parent's `voltage_kp` can be handed down to the generic `Kp` of a nested PI block (`@named pi = PIBlock(Kp = voltage_kp)`), which demotes `pi₊Kp` and leaves `voltage_kp` as the single tunable knob — rather than two separate parameters that must be kept in sync.

#### [Backward-Flow Initialization: Formulas That Find Each Other](@id backward-flow-init)

The real power of init formulas is that they **connect to one another**. Each formula is a small "set this from that" rule; NetworkDynamics gathers the rules of an entire nested model and wires them into a single dependency graph, automatically matching each value a formula *produces* to the formulas that *consume* it. You never connect them by hand — the outputs and inputs find each other by name (across aliases, too) and are evaluated in the right order.

This lets initialization flow **backwards** through a model, against the direction the equations normally compute. The motivating case is a control block that can initialize its own state *from its own output*. A PI controller's integrator state, in steady state, is:

```julia
x(t), [guess=0, initf = (y - K_p*err)/K_i]   # the PI law y = K_p⋅err + K_i⋅x, solved for x
```

On its own this rule is stuck: `y` is the block's *output*, normally computed *from* `x`, so the block cannot supply `y` by itself. The missing piece — what `y` has to be in steady state — is known only one layer up, where the surrounding component sets it with [`set_initf`](@ref):

```julia
@component function SimpleAVR(; name)
    @named pi = PIBlock()      # carries the x-from-y initf above
    # ... equations ...
    sys = System(eqs, t; name, systems=[pi])
    # in steady state the PI must hold the field circuit: set its output
    sys = set_initf(sys, pi.y => K_e*v_f)
end
```

Now the two rules link up: the parent's rule produces `y`, the child's rule consumes `y` and produces `x`. A dead end becomes a chain that determines every internal state from the operating point at the interface — the component summary lists the resulting formula clusters and `initialize_component` reports *"No free variables!"*.

**Setting an observable (pinning).** For the chain to close, a formula has to be allowed to set `y` — an *observable*, a quantity the model defines through an equation rather than storing it as a state. NetworkDynamics permits this and calls it *pinning*. A pinned value is only a stepping stone used while the formulas run; it is not stored on the component. Afterwards it is checked for consistency: the observable is recomputed from the solved state and must agree with what was pinned, otherwise you get a warning (and usually a failing residual). A plain default sitting on an observable is never used this way — pinning is a deliberate act of a formula, not ambient data.

**Guesses do the same, more softly.** A `GuessFormula` can drive the very same backward chain, but its values land among the *guesses*: they seed the solver near the right operating point without committing to it, and are never consistency-checked. That makes back-computing guess formulas a safe default to sprinkle throughout a model library — spell the whole backward chain as guesses and initialization simply starts close to the answer. Wherever you *do* know an exact value (a default, or an `InitFormula` pin), it takes precedence over the guessed one. The guess-side spellings mirror the init side exactly: the `guessf` variable option next to a declaration, and [`set_guessf`](@ref) to attach a rule from outside — see below.

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

**From MTK models**: guess formulas mirror the `initf` / [`set_initf`](@ref) spellings on the guess side. Next to a declaration use the `guessf` variable option; from outside the block use [`set_guessf`](@ref):

```julia
@component function avr(; name)
    @variables begin
        v_meas(t), [guess=1.0, guessf=v_mag]   # guess the measurement from the actual value
    end
    # ...
end

sys = set_guessf(sys, sub.x => sub.p)          # seed a subsystem's state from outside
```

`guessf` reads as "at initialization, *guess* this symbol from that expression". Because a guess is only a hint, it differs from `initf` in two ways: conflicting recipes for one target are a warning (not an error), and a formula whose inputs cannot be resolved is silently skipped. That skip is what makes `guessf` compose with a scalar `guess`: give a variable both `guess=0` (the fallback seed) and `guessf=<expr>` (the refined value), and the formula is used whenever it resolves while the scalar remains as a safe default. A symbolic value given to the plain `guess` option is rejected — spell it as `guessf` instead.


## Analysing Fixpoints
In order to analyse fixpoints NetworkDynamis provides the functions [`isfixpoint`](@ref), [`is_linear_stable`](@ref) and [`jacobian_eigenvals`](@ref).
