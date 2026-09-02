# # Tutorial on Stepwise Initialization of a Complex Model
#
# This example demonstrates how to initialize a complex network model with both static
# and dynamic components.
# The models are closely related to the ones used in the gas network example,
# but greatly simplified for the sake of this tutorial.
# We'll create a gas network model with three nodes
# and pipes connecting them, and show how to:
#
# 1. Create static models for initialization
# 2. Find a steady-state solution
# 3. Create corresponding dynamic models
# 4. Initialize the dynamic models with the steady-state solution
# 5. Simulate the system with dynamic behavior
#
#
# First, let's import the necessary packages:

using NetworkDynamics
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using OrdinaryDiffEqTsit5
using CairoMakie
nothing #hide

# ## Node Models
#
# We'll start by defining our node models using ModelingToolkit.
# First, let's create a template for common states and equations in all gas nodes:

@mtkmodel GasNode begin
    @variables begin
        p(t), [description="Pressure"] # node output
        q̃_nw(t), [description="aggregated flow from pipes into node"] # node input
        q̃_inj(t), [description="flow injected into the network"]
    end
    @equations begin
        q̃_inj ~ -q̃_nw
    end
end
nothing #hide

# Now we'll define three specific node types:
#
# **A) A constant pressure node that forces pressure to maintain a specific value**

@mtkmodel ConstantPressureNode begin
    @extend GasNode()
    @parameters begin
        p_set, [description="Constant pressure setpoint"]
    end
    @equations begin
        p ~ p_set
    end
end
nothing #hide

# **B) A static prosumer node which forces a certain flow**
#
# > **Fully Implicit Output**
# >
# > We need to use `implicit_output(p)` to handle the fully implicit pressure
# > output. See fully implicit outputs for details.

@mtkmodel StaticProsumerNode begin
    @extend GasNode()
    @parameters begin
        q̃_prosumer, [description="flow injected by prosumer"]
    end
    @equations begin
        -q̃_nw ~ q̃_prosumer + implicit_output(p)
    end
end
nothing #hide

# **C) A dynamic prosumer node with compliance, which introduces dynamic behavior to the pressure state**

@mtkmodel DynamicProsumerNode begin
    @extend GasNode()
    @parameters begin
        q̃_prosumer, [description="flow injected by prosumer"]
        C=0.1, [description="Compliance"]
    end
    @equations begin
        C*D(p) ~ q̃_prosumer + q̃_nw
    end
end
nothing #hide

# **D) A pressure control node that tries to maintain a set pressure by adjusting its injection**

@mtkmodel PressureControlNode begin
    @extend GasNode()
    @parameters begin
        p_set, [description="Pressure setpoint", guess=1]
        K_p=1, [description="Proportional gain"]
        K_i=1, [description="Integral gain"]
        C=0.1, [description="Compliance"]
    end
    @variables begin
        Δp(t), [description="Pressure error"]
        ξ(t), [description="Integral state", guess=0]
        q̃_prosumer(t), [description="flow injected by producer"]
    end
    @equations begin
        Δp ~ p_set - p
        D(ξ) ~ Δp
        q̃_prosumer ~ K_p*Δp + K_i*ξ
        C*D(p) ~ q̃_prosumer + q̃_nw
    end
end
nothing #hide

# ## Edge Models
#
# Now we'll define our edge models, starting with a template for the pipe:

@mtkmodel GasPipe begin
    @variables begin
        q̃(t), [description="flow through pipe"] #output
        p_src(t), [description="pressure at start of pipe"] #input
        p_dst(t), [description="pressure at end of pipe"] #input
    end
end
nothing #hide

# Next, we define a dynamic pipe with inertia (a simple delayed model):

@mtkmodel DynamicPipe begin
    @extend GasPipe()
    @parameters begin
        R=0.1, [description="Resistance"]
        M=0.1, [description="Inertia"]
    end
    @equations begin
        M*D(q̃) ~ (p_src - p_dst)/R - q̃ # some simple delayed model
    end
end
nothing #hide

# And finally, a quasistatic pipe model for initialization purposes. This model equals the
# dynamic model at steady state, making it ideal for finding initial conditions:

@mtkmodel QuasistaticPipe begin
    @extend GasPipe()
    @parameters begin
        R=0.1, [description="Resistance"]
    end
    @equations begin
        q̃ ~ (p_src - p_dst)/R
    end
end
nothing #hide

# ## Defining a Static Model for Initialization
#
# Our first step is to define a static model that we'll use to find the steady-state solution.
# This is a crucial step for initializing complex dynamic models.
#
# Step 1: Define all the components of our static model
# First, node 1 is our producer which will later be a controlled producer. For initialization, we use a static model:

@named v1_mod_static = ConstantPressureNode(p_set=1)
v1_static = VertexModel(v1_mod_static, [:q̃_nw], [:p], vidx=1)

# Nodes 2 and 3 are consumers. For them, we'll use static prosumer models:
@named v2_mod_static = StaticProsumerNode(q̃_prosumer=-0.6) # consumer
v2_static = VertexModel(v2_mod_static, [:q̃_nw], [:p], vidx=2)

@named v3_mod_static = StaticProsumerNode(q̃_prosumer=-0.4) # consumer
v3_static = VertexModel(v3_mod_static, [:q̃_nw], [:p], vidx=3)
nothing #hide

# Now we define the static pipe models connecting our nodes:

@named p_mod_static = QuasistaticPipe()
p12_static = EdgeModel(p_mod_static, [:p_src], [:p_dst], AntiSymmetric([:q̃]), src=1, dst=2)
p13_static = EdgeModel(p_mod_static, [:p_src], [:p_dst], AntiSymmetric([:q̃]), src=1, dst=3)
p23_static = EdgeModel(p_mod_static, [:p_src], [:p_dst], AntiSymmetric([:q̃]), src=2, dst=3)
nothing #hide

# Assemble all components into a static network:

nw_static = Network([v1_static, v2_static, v3_static], [p12_static, p13_static, p23_static])

# Create an initial guess for the steady state and modify it with reasonable values:

u_static_guess = NWState(nw_static)
u_static_guess.v[2, :p] = 1.0
u_static_guess.v[3, :p] = 1.0
nothing #hide

# Find the steady-state solution using our initial guess:

u_static = find_fixpoint(nw_static, u_static_guess)

# ## Defining a Dynamic Model
#
# Now we'll define our dynamic model using more complex components:

@named v1_mod_dyn = PressureControlNode()
v1_dyn = VertexModel(v1_mod_dyn, [:q̃_nw], [:p], vidx=1)

@named v2_mod_dyn = DynamicProsumerNode(q̃_prosumer=-0.6)
v2_dyn = VertexModel(v2_mod_dyn, [:q̃_nw], [:p], vidx=2)

@named v3_mod_dyn = DynamicProsumerNode(q̃_prosumer=-0.4)
v3_dyn = VertexModel(v3_mod_dyn, [:q̃_nw], [:p], vidx=3)
nothing #hide

# Create dynamic pipe models with inertia:

@named p_mod_dyn = DynamicPipe()
p12_dyn = EdgeModel(p_mod_dyn, [:p_src], [:p_dst], AntiSymmetric([:q̃]), src=1, dst=2)
p13_dyn = EdgeModel(p_mod_dyn, [:p_src], [:p_dst], AntiSymmetric([:q̃]), src=1, dst=3)
p23_dyn = EdgeModel(p_mod_dyn, [:p_src], [:p_dst], AntiSymmetric([:q̃]), src=2, dst=3)
nothing #hide

# Assemble the dynamic network:

nw_dyn = Network([v1_dyn, v2_dyn, v3_dyn], [p12_dyn, p13_dyn, p23_dyn])

# ## Initializing the Dynamic Model with the Static Solution
#
# Now comes the important part: we need to initialize the dynamic model with the results from
# the static model. To do so, we need to make use of two functions:
#
# 1. `interface_values`: Extracts all interface values (inputs and outputs) from a network state
# 2. `initialize_componentwise!`: Initializes all components in a network one by one
#
# First, we extract all interface values from our static solution:

interface_vals = interface_values(u_static)

# Next, we initialize the dynamic model using these interface values as defaults:

u0_dyn = initialize_componentwise!(nw_dyn, default_overrides=interface_vals, verbose=true)

# Internally, this function uses `initialize_component!` on every single component.
# For each component, it overwrites the `default`s to be consistent with the interface values
# of the static model. Therefore, we ensure that the dynamic model is initialized near the
# steady state of the static model.
#
# We can inspect individual components if needed:

dump_initial_state(nw_dyn[VIndex(1)])
nothing #hide

# Let's verify that our initialization is correct by checking that the derivatives are close to zero:

du = ones(dim(nw_dyn))
nw_dyn(du, uflat(u0_dyn), pflat(u0_dyn), 0.0)
extrema(du .- zeros(dim(nw_dyn))) # very close to zero, confirming we have a steady state!

# Alternatively, we can used the `isfixpoint` function to check if the state is a fixpoint:

@assert isfixpoint(nw_dyn, u0_dyn)
nothing #hide

# ## Simulating the Dynamic Model
#
# With our properly initialized model, we can now simulate the system to observe its behavior.
# To make the simulation more interesting, we'll introduce a disturbance to see how the system
# responds from its steady state.
#
# We'll use a callback to increase consumer demand at a specific time. For more information on
# callbacks, see the documentation on Callbacks.

affect = ComponentAffect([], [:q̃_prosumer]) do u, p, ctx
    @info "Increase consumer demand at t=$(ctx.t)"
    p[:q̃_prosumer] -= 0.1
end
cb = PresetTimeComponentCallback([1.0], affect)
set_callback!(nw_dyn[VIndex(2)], cb) # attach disturbance to second node
nothing #hide

# Create and solve the ODE problem with the callback. Note that we're using the flat
# representation of our initialized state (via `uflat` and `pflat`) as input to the ODE solver:

prob = ODEProblem(nw_dyn, copy(uflat(u0_dyn)), (0, 7), copy(pflat(u0_dyn));
    callback=get_callbacks(nw_dyn))
sol = solve(prob, Tsit5())
nothing #hide

# ## Visualizing the Results
#
# Finally, let's visualize the results of our simulation.
# The plots show how our gas network responds to the increased consumer demand at t=1:
#
# 1. **Pressure at nodes**: We see a pressure drop at all nodes after the disturbance before the pressure is stabilized by the controller.
#
# 2. **Injection by producer**: Node 1 increases its injection to compensate for the higher demand.
#
# 3. **Consumption by consumers**: The solid lines show the actual flows at nodes 2 and 3, while the dashed lines show the set consumer demands. At t=1, we see the step change in consumer demand at node 2.
#
# 4. **Flows through pipes**: Shows how the flows in all pipes adjust to the new demand pattern.

let
    fig = Figure(size=(1000,1000))
    ax = Axis(fig[1, 1]; title="Pressure at nodes")
    for i in 1:3
        lines!(ax, sol; idxs=VIndex(i, :p), label="Node $i", color=Cycled(i))
    end

    ax = Axis(fig[2, 1]; title="Injection by producer")
    lines!(ax, sol; idxs=VIndex(1, :q̃_inj), label="Node 1", color=Cycled(1))

    ax = Axis(fig[3, 1]; title="Consumption by consumers")
    for i in 2:3
        lines!(ax, sol; idxs=@obsex(-1*VIndex(i, :q̃_inj)), label="Node $i", color=Cycled(i))
        lines!(ax, sol; idxs=@obsex(-1*VIndex(i, :q̃_prosumer)), label="Node $i", linestyle=:dash, color=Cycled(i))
    end

    ax = Axis(fig[4, 1]; title="Flows through pipes")
    for i in 1:3
        lines!(ax, sol; idxs=@obsex(abs(EIndex(i, :q̃))), label="Pipe $i", color=Cycled(i))
    end

    fig
end

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
