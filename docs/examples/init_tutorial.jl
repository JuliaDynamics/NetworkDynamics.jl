#=
# [Tutorial on Stepwise Initialization of a Complex Model](@id init-tutorial)

This example demonstrates how to initialize a complex network model with both static
and dynamic components.
The models are closely related to the ones used in the [gas network example](@ref gas-example),
but greatly simplified for the sake of this tutorial.
We'll create a gas network model with three nodes
and pipes connecting them, and show how to:

1. Create static models for initialization
2. Find a steady-state solution
3. Create corresponding dynamic models
4. Initialize the dynamic models with the steady-state solution
5. Simulate the system with dynamic behavior

This script can be downloaded as a normal Julia script [here](@__NAME__.jl). #md

First, let's import the necessary packages:
=#

using NetworkDynamics
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using OrdinaryDiffEqTsit5
using CairoMakie
nothing #hide

#=
## Node Models

We'll start by defining our node models using ModelingToolkit.
First, let's create a template for common states and equations in all gas nodes:
=#
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

#=
Now we'll define three specific node types:

**A) A constant pressure node that forces pressure to maintain a specific value**
=#
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

#=
**B) A static prosumer node which forces a certain flow (pressure is fully implicit)**
=#
@mtkmodel StaticProsumerNode begin
    @extend GasNode()
    @parameters begin
        q̃_prosumer, [description="flow injected by prosumer"]
    end
    @equations begin
        -q̃_nw ~ q̃_prosumer
    end
end
nothing #hide

#=
**C) A dynamic prosumer node with compliance, which adds dynamics to the pressure state**
=#
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

#=
**D) A pressure control node that tries to maintain a set pressure by adjusting its injection**
=#
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

#=
## Edge Models

Now we'll define our edge models, starting with a template for the pipe:
=#
@mtkmodel GasPipe begin
    @variables begin
        q̃(t), [description="flow through pipe"] #output
        p_src(t), [description="pressure at start of pipe"] #input
        p_dst(t), [description="pressure at end of pipe"] #input
    end
end
nothing #hide

#=
Next, we define a dynamic pipe with inertia (a simple delayed model):
=#
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

#=
And finally a quasistatic pipe model for initialization purposes. This equals the 
dynamic model in steady state, making it ideal for finding initial conditions:
=#
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

#=
## Defining a Static Model for Initialization

Our first step is to define a static model that we'll use to find the steady-state solution.
This is a crucial step for initializing complex dynamic models.

Step 1: Define all the components of our static model
First, node 1 is our producer which will later be a controlled producer. For initialization, we use a static model:
=#
@named v1_mod_static = ConstantPressureNode(p_set=1)
v1_static = VertexModel(v1_mod_static, [:q̃_nw], [:p], vidx=1)

## Nodes 2 and 3 are consumers. For them, we'll use static prosumer models:
@named v2_mod_static = StaticProsumerNode(q̃_prosumer=-0.6) # consumer
v2_static = VertexModel(v2_mod_static, [:q̃_nw], [:p], vidx=2)

@named v3_mod_static = StaticProsumerNode(q̃_prosumer=-0.4) # consumer
v3_static = VertexModel(v3_mod_static, [:q̃_nw], [:p], vidx=3)
nothing #hide

#=
Now we define the static pipe models connecting our nodes:
=#
@named p_mod_static = QuasistaticPipe()
p12_static = EdgeModel(p_mod_static, [:p_src], [:p_dst], AntiSymmetric([:q̃]), src=1, dst=2)
p13_static = EdgeModel(p_mod_static, [:p_src], [:p_dst], AntiSymmetric([:q̃]), src=1, dst=3)
p23_static = EdgeModel(p_mod_static, [:p_src], [:p_dst], AntiSymmetric([:q̃]), src=2, dst=3)
nothing #hide
#=
Assemble all components into a static network:
=#
nw_static = Network([v1_static, v2_static, v3_static], [p12_static, p13_static, p23_static])

#=
Create an initial guess for the steady state and modify it with reasonable values:
=#
u_static_guess = NWState(nw_static)
u_static_guess.v[2, :p] = 1.0
u_static_guess.v[3, :p] = 1.0
nothing #hide

#=
Find the steady-state solution using our initial guess:
=#
u_static = find_fixpoint(nw_static, u_static_guess)

#=
## Defining a Dynamic Model

Now we'll define our dynamic model using more complex components:
=#
@named v1_mod_dyn = PressureControlNode()
v1_dyn = VertexModel(v1_mod_dyn, [:q̃_nw], [:p], vidx=1)

@named v2_mod_dyn = DynamicProsumerNode(q̃_prosumer=-0.6)
v2_dyn = VertexModel(v2_mod_dyn, [:q̃_nw], [:p], vidx=2)

@named v3_mod_dyn = DynamicProsumerNode(q̃_prosumer=-0.4)
v3_dyn = VertexModel(v3_mod_dyn, [:q̃_nw], [:p], vidx=3)
nothing #hide

#=
Create dynamic pipe models with inertia:
=#
@named p_mod_dyn = DynamicPipe()
p12_dyn = EdgeModel(p_mod_dyn, [:p_src], [:p_dst], AntiSymmetric([:q̃]), src=1, dst=2)
p13_dyn = EdgeModel(p_mod_dyn, [:p_src], [:p_dst], AntiSymmetric([:q̃]), src=1, dst=3)
p23_dyn = EdgeModel(p_mod_dyn, [:p_src], [:p_dst], AntiSymmetric([:q̃]), src=2, dst=3)
nothing #hide

#=
Assemble the dynamic network:
=#
nw_dyn = Network([v1_dyn, v2_dyn, v3_dyn], [p12_dyn, p13_dyn, p23_dyn])

#=
## Initializing the Dynamic Model with the Static Solution

Now comes the important part: we need to initialize the interface values (pressures and flows)
of the dynamic model with the results from the static model.

We can do this manually:
=#
## Vertex 1: output
set_default!(nw_dyn[VIndex(1)], :p, u_static.v[1, :p])
## Vertex 1: input
set_default!(nw_dyn[VIndex(1)], :q̃_nw, u_static.v[1, :q̃_nw])
nothing #hide

#=
But there is also a built-in method [`set_interface_defaults!`](@ref) which we can use
automatically:
=#
set_interface_defaults!(nw_dyn, u_static; verbose=true)
nothing #hide

#=
With the interfaces all set, we can "initialize" the internal states of the dynamic models.

For example, let's inspect the state of our first vertex:
=#
nw_dyn[VIndex(1)]

#=
We observe that both the initial state `ξ` as well as the pressure setpoint `p_set`
are left "free". Using [`initialize_component!`](@ref), we can try to find values for the
"free" states and parameters such that the interface constraints are fulfilled.
=#
initialize_component!(nw_dyn[VIndex(1)])
#=
We may also use [`dump_initial_state`](@ref) to get a more detailed view of the state:
=#
dump_initial_state(nw_dyn[VIndex(1)])
nothing #hide

#=
We can also initialize the other two vertices, however it is unnecessary
since their state is already completely determined by the fixed input/output:
=#
initialize_component!(nw_dyn[VIndex(2)])
initialize_component!(nw_dyn[VIndex(3)])
nothing #hide

#=
Similarly, we can initialize the dynamic pipe models. However, since their dynamic state
equals the output, once again there is nothing to initialize.
=#
initialize_component!(nw_dyn[EIndex(1)])
initialize_component!(nw_dyn[EIndex(2)])
initialize_component!(nw_dyn[EIndex(3)])
nothing #hide

#=
Now, everything is initialized, which means every input, output, state and parameter
either has a `default` metadata or an `init` metadata. When constructing the `NWState`
for this network, it will be filled with all those values which should now correspond
to a steady state of the system:
=#
u0_dyn = NWState(nw_dyn)

#=
Let's verify that our initialization is correct by checking that the derivatives are close to zero:
=#
du = ones(dim(nw_dyn))
nw_dyn(du, uflat(u0_dyn), pflat(u0_dyn), 0.0)
extrema(du .- zeros(dim(nw_dyn))) # very close to zero, confirming we have a steady state!

#=
## Simulating the Dynamic Model

Now we can solve the dynamic model and add a disturbance to see how the system responds:
=#
affect = ComponentAffect([], [:q̃_prosumer]) do u, p, ctx
    @info "Increase consumer demand at t=$(ctx.t)"
    p[:q̃_prosumer] -= 0.1
end
cb = PresetTimeComponentCallback([1.0], affect)
set_callback!(nw_dyn[VIndex(2)], cb) # attach disturbance to second node
nothing #hide

#=
Create and solve the ODE problem with the callback:
=#
prob = ODEProblem(nw_dyn, copy(uflat(u0_dyn)), (0, 7), copy(pflat(u0_dyn));
    callback=get_callbacks(nw_dyn))
sol = solve(prob, Tsit5())
nothing #hide

#=
## Visualizing the Results

Finally, let's visualize the results of our simulation.
The plots show how our gas network responds to the increased consumer demand at t=1:

1. **Pressure at nodes**: We see a pressure drop at all nodes after the disturbance before the pressure is stabilized by the controller.

2. **Injection by producer**: Node 1 increases its injection to compensate for the higher demand.

3. **Draw by consumers**: The solid lines show the actual flows at nodes 2 and 3, while the dashed lines show the set consumer demands. At t=1, we see the step change in consumer demand at node 2.

4. **Flows through pipes**: Shows how the flows in all pipes adjust to the new demand pattern.
=#

let
    fig = Figure(size=(1000,1000))
    ax = Axis(fig[1, 1]; title="Pressure at nodes")
    for i in 1:3
        lines!(ax, sol; idxs=VIndex(i, :p), label="Node $i", color=Cycled(i))
    end

    ax = Axis(fig[2, 1]; title="Injection by producer")
    lines!(ax, sol; idxs=VIndex(1, :q̃_inj), label="Node 1", color=Cycled(1))

    ax = Axis(fig[3, 1]; title="Draw by consumers")
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
