# Tutorial on stepwise initialization of a complex model

using NetworkDynamics
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using OrdinaryDiffEqTsit5
using CairoMakie

## Node Models
# template for commont states and equations in all gas nodes
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

# node that forces pressure to be constant
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

# node which forces a certain flow (pressure fully implicit)
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

# dynamic prosumer node with compliance
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

# node with a pressure controller
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
# template for pipe
@mtkmodel GasPipe begin
    @variables begin
        q̃(t), [description="flow through pipe"] #output
        p_src(t), [description="pressure at start of pipe"] #input
        p_dst(t), [description="pressure at end of pipe"] #input
    end
end
nothing #hide

# pipe with a simple delayed model
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

# quasistatic pipe model for initialization (equals delayed model in steady state!)
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
## Define a static model
=#

# step 1: definition of a static model
# node 1 is our producer which will later be a controlled producer. For initialization we use a static model
@named v1_mod_static = ConstantPressureNode(p_set=1)
v1_static = VertexModel(v1_mod_static, [:q̃_nw], [:p], vidx=1)

# node 2 and 3 are consumers, there is no difference between static and dynamic model for them
@named v2_mod_static = StaticProsumerNode(q̃_prosumer=-0.6) # consumer, initialize pressure
v2_static = VertexModel(v2_mod, [:q̃_nw], [:p], vidx=2)

@named v3_mod_static = StaticProsumerNode(q̃_prosumer=-0.4) # consumer
v3_static = VertexModel(v3_mod, [:q̃_nw], [:p], vidx=3)

# static pipes
@named p_mod_static = QuasistaticPipe()
p12_static = EdgeModel(p_mod_static, [:p_src], [:p_dst], AntiSymmetric([:q̃]), src=1, dst=2)
p13_static = EdgeModel(p_mod_static, [:p_src], [:p_dst], AntiSymmetric([:q̃]), src=1, dst=3)
p23_static = EdgeModel(p_mod_static, [:p_src], [:p_dst], AntiSymmetric([:q̃]), src=2, dst=3)

nw_static = Network([v1_static, v2_static, v3_static], [p12_static, p13_static, p23_static])

u_static_guess = NWState(nw_static)
u_static_guess.v[2, :p] = 1.0
u_static_guess.v[3, :p] = 1.0

u_static = find_fixpoint(nw_static, u_static_guess)

# ## Define a dynamic model
@named v1_mod_dyn = PressureControlNode(;p_set=1)
v1_dyn = VertexModel(v1_mod_dyn, [:q̃_nw], [:p], vidx=1)

@named v2_mod_dyn = DynamicProsumerNode(q̃_prosumer=-0.6)
v2_dyn = VertexModel(v2_mod_dyn, [:q̃_nw], [:p], vidx=2)

@named v3_mod_dyn = DynamicProsumerNode(q̃_prosumer=-0.4)
v3_dyn = VertexModel(v3_mod_dyn, [:q̃_nw], [:p], vidx=3)

@named p_mod_dyn = DynamicPipe()
p12_dyn = EdgeModel(p_mod_dyn, [:p_src], [:p_dst], AntiSymmetric([:q̃]), src=1, dst=2)
p13_dyn = EdgeModel(p_mod_dyn, [:p_src], [:p_dst], AntiSymmetric([:q̃]), src=1, dst=3)
p23_dyn = EdgeModel(p_mod_dyn, [:p_src], [:p_dst], AntiSymmetric([:q̃]), src=2, dst=3)

nw_dyn = Network([v1_dyn, v2_dyn, v3_dyn], [p12_dyn, p13_dyn, p23_dyn])

# no we need to do some magic. we want to initialize the interface values (pressures and flow)
# of the dynamic model with the results of the static model

set_default!(nw_dyn[VIndex(1)], :p, u_static.v[1, :p])
set_default!(nw_dyn[VIndex(1)], :q̃_nw, u_static.v[1, :q̃_nw])
v1_dyn
# in the output you can see now that state ξ is "approx" 1 (guess value) while the rest is fixed
# so we can initialize the ydnamic model
initialize_component!(v1_dyn)
v1_dyn

# the ther two vertices are simplere, the we need t oset the default values
set_default!(nw_dyn[VIndex(2)], :p, u_static.v[2, :p])
set_default!(nw_dyn[VIndex(3)], :p, u_static.v[3, :p])

# we can manually "initialize" the only open state of the dynamic line model
set_default!(nw_dyn[EIndex(1)], :q̃, u_static.e[1, :q̃])
set_default!(nw_dyn[EIndex(2)], :q̃, u_static.e[2, :q̃])
set_default!(nw_dyn[EIndex(3)], :q̃, u_static.e[3, :q̃])

# we have set all the "default" values for all teh states now. so we call `NWState` on the
# network we should get a fully initialized state
u0_dyn = NWState(nw_dyn)

du = ones(dim(nw_dyn))
nw_dyn(du, uflat(u0_dyn), pflat(u0_dyn), 0.0)
extrema(du .- zeros(dim(nw_dyn))) # very close to zero, is steady state!

# now we can solve the dynamic model
affect = ComponentAffect([], [:q̃_prosumer]) do u, p, ctx
    @info "Increas consumer at t=$(ctx.t)"
    p[:q̃_prosumer] -= 0.1
end
cb = PresetTimeComponentCallback([1.0], affect)
set_callback!(nw_dyn[VIndex(2)], cb) # attach disturabnce to second node

prob = ODEProblem(nw_dyn, copy(uflat(u0_dyn)), (0, 7), copy(pflat(u0_dyn));
    callback=get_callbacks(nw_dyn))
sol = solve(prob, Tsit5())

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


# btw have you tried the  GUi yet? :)
# using NetworkDynamicsInspector
# inspect(sol)
