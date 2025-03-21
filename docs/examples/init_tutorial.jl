#=
# Tutorial on Stepwise Initialization of a Complex Model

This example demonstrates how to initialize a complex network model with both static
and dynamic components. We'll create a simple gas network model with three nodes
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

1. A constant pressure node that forces pressure to maintain a specific value
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
2. A static prosumer node which forces a certain flow (pressure is fully implicit)
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
3. A dynamic prosumer node with compliance, which adds dynamics to the pressure state
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
4. A pressure control node that tries to maintain a set pressure by adjusting its injection
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
@named v1_mod_dyn = PressureControlNode(;p_set=1)
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

First, let's handle node 1 (the pressure control node):
=#

set_default!(nw_dyn[VIndex(1)], :p, u_static.v[1, :p])
set_default!(nw_dyn[VIndex(1)], :q̃_nw, u_static.v[1, :q̃_nw])
v1_dyn # hide
#=
In the output, you can see that state ξ is "approx" 1 (guess value) while the rest is fixed.
Now we can initialize the dynamic model's internal states:
=#
initialize_component!(v1_dyn)
v1_dyn #hide

#=
For the other two vertices (which are simpler), we just need to set the default values:
=#
set_default!(nw_dyn[VIndex(2)], :p, u_static.v[2, :p])
set_default!(nw_dyn[VIndex(3)], :p, u_static.v[3, :p])
nothing #hide

#=
For the pipe models, we manually initialize the flow state of the dynamic line model:
=#
set_default!(nw_dyn[EIndex(1)], :q̃, u_static.e[1, :q̃])
set_default!(nw_dyn[EIndex(2)], :q̃, u_static.e[2, :q̃])
set_default!(nw_dyn[EIndex(3)], :q̃, u_static.e[3, :q̃])
nothing #hide

#=
Now that we've set all the "default" values for all the states, we can call `NWState` on the
network to get a fully initialized state vector:
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

#=
## Interactive Visualization

You can also visualize the results interactively using NetworkDynamicsInspector:

```julia
using NetworkDynamicsInspector
inspect(sol)
```
=#
