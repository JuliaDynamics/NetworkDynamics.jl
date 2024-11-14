#=
# Dynamic Flow in simple Gas Network

This Example is based on the paper

    Albertus J. Malan, Lukas Rausche, Felix Strehle, Sören Hohmann,
    Port-Hamiltonian Modelling for Analysis and Control of Gas Networks,
    IFAC-PapersOnLine, Volume 56, Issue 2, 2023, https://doi.org/10.1016/j.ifacol.2023.10.193.

and tries replicate a simple simulation of flow in a 3-node gas network.

This example can be dowloaded as a normal Julia script [here](@__NAME__.jl). #md

We start by importing the necessary packages:
=#
using NetworkDynamics
using ModelingToolkit
using DynamicQuantities
using ModelingToolkit: D as Dt, t as t
using Test
using StaticArrays
using LinearInterpolations
using OrdinaryDiffEqTsit5
using CairoMakie
CairoMakie.activate!(type="svg") #hide

#=
## Node Models

There are 2 node models used in the paper. The first node type has a constant
pressure. We model this as an first order ODE with zero rhs.

In this example, we use the equation based modeling using `ModelingToolkit.jl`.
To verify the equations on a basic level we also provide units to eveything to
perform dimensionality checks.
=#
@mtkmodel ConstantPressureNode begin
    @parameters begin
        p_set, [description="Constant pressure setpoint", unit=u"Pa"]
    end
    @variables begin
        p(t) = p_set, [description="Pressure", unit=u"Pa", output=true]
        q̃_nw(t), [description="aggregated flow from pipes into node", unit=u"m^3/s", input=true]
        q̃_inj(t), [description="internal help variable", unit=u"m^3/s"]
    end
    @equations begin
        Dt(p) ~ 0
        q̃_inj ~ -q̃_nw
    end
end
nothing #hide

#=
The second node model is a variable pressure node. It has one output state, the pressure and one input state,
the aggregated flows from the connected pipes.
As an internal state we have the injected flow from our source/load.
The source/load behaviour itself is provided via a time dependent function.
=#
@mtkmodel VariablePressureNode begin
    @structural_parameters begin
        load_profile # time dependent load profile
    end
    @parameters begin
        C, [description="Lumped capacitance of connected pipes", unit=u"m^4 * s^2 / kg"]
    end
    @variables begin
        p(t)=5e6, [description="Pressure", unit=u"Pa", output=true]
        q̃_inj(t), [description="external injection into node", unit=u"m^3/s"]
        q̃_nw(t), [description="aggregated flow from pipes into node", unit=u"m^3/s", input=true]
    end
    @equations begin
        q̃_inj ~ load_profile(t)
        C * Dt(p) ~ q̃_inj + q̃_nw # (30)
    end
end
nothing #hide

#=
## Pipe Model

The pipe is modeld as a first order ODE for the volumetric flow at the `dst` end. It has two inputs:
the pressure at the source and and the pressure at the destination end.
Later on, we'll specify the model to be antisymmetric, thus the flow is calculated explicitly for the
destination end, but the source end will just recive just that times (-1).
=#
@mtkmodel Pipe begin
    @parameters begin
        L, [description="Length of pipe", unit=u"m"]
        A, [description="Cross-sectional area of pipe", unit=u"m^2"]
        ρ̃, [description="constant standard density", unit=u"kg/m^3"]
        D, [description="Diameter of pipe", unit=u"m"]
        g, [description="Gravitational acceleration", unit=u"m/s^2"]
        sinθ, [description="Angle of inclination" ]
        c, [description="Speed of sound in fluid", unit=u"m/s"]
        γ, [description="Friction efficiency factor"]
        η, [description="Dynamic viscosity", unit=u"kg/(m*s)"]
        r, [description="Pipe roughness", unit=u"m"]
    end
    @variables begin
        p_src(t), [description="Pressure at source end", unit=us"Pa", input=true]
        p_dst(t), [description="Pressure at destination end", unit=us"Pa", input=true]
        q̃(t)=1, [description="Flow through pipe", unit=u"m^3/s", output=true]
        Re(t), [description="Reynolds number"]
        λ(t), [description="Friction factor"]
        λe(t), [description="Effective friction factor"]
        pM(t), [description="mean pressure", unit=us"Pa"]
    end
    @equations begin
        Re ~ (ρ̃ * abs(q̃) * D) / (η * A) # (6)
        λ ~ ifelse(Re < 2300,
            64/Re, # laminar (7)
            (2*log10(4.518/Re * log10(Re/7) + r/(3.71*D)))^(-2) # turbulent (8)
        )
        λe ~ λ/γ^2 # (10)
        pM ~ 2/3*(p_src + p_dst - (p_src*p_dst)/(p_src + p_dst)) # (20)

        ρ̃ * L / A * Dt(q̃) ~ -(λe * ρ̃^2 * c^2 * L * abs(q̃))/(2 * D * A^2 * pM) * q̃ - (g * L * sinθ)/(c^2) * pM + (p_src - p_dst) # (31)
    end
end
nothing #hide

#=
## Parametrization

The parameterization turned out to be a bit tricky. There might be errors in there.

Some of them are quite cleare and explicitly given.
=#
g = 9.81u"m/s^2"       # that we just know
Rs = 518.28u"J/(kg*K)" # Specific gas constant for hydrogen
η  = 1e-5u"kg/(m*s)"   # Dynamic viscosity
pc = 46.5u"bar"        # Critical pressure
p̃  = 1.01325u"bar"     # standard pressure
Tc = 190.55u"K"        # critical temperature
T̃  = 273.15u"K"        # standard temperature
γ  = 0.98              # friction efficiency factor
r  = 0.012u"mm"        # pipe roughness
D  = 0.6u"m"           # pipe diameter

L₁₂ = 80u"km"
L₁₃ = 90u"km"
L₂₃ = 100u"km"
Δh₁ = 0u"km"           # this value is different for different sims in the paper
p₁_set = 50u"bar"
nothing # hide

#=
Some neede parameters for the equations are easily derived.
=#
A = π/4 * D^2
sinθ₁₂ = ustrip(Δh₁ / L₁₂)
sinθ₁₃ = ustrip(Δh₁ / L₁₃)
sinθ₂₃ = 0.0
nothing # hide

#=
Other parameters are more complicated. The compressibility factor Z according to (5)
depends on temperature and pressure. While the temperature is not considered in the model
at all, I am not sure whether we sould include the pressure dependency.
Here, i just calculated it using the "standard" pressure and "standard" temperature.
=#
Z̃ = 1 - 3.52 * p₁_set/pc * exp(-2.26*(T̃/Tc)) + 0.274 * (p₁_set/pc)^2 * exp(-1.878*(T̃/Tc)) # (5)
nothing # hide

#=
Similarily, we can use (4) to derive the speed of sound. This again is
temperature and compressibility dependent. Here I use the "standard" compressibility again
to calculate a constant "standard" sound speed.
Equivalently, we use the "standard" compressibility to calculate the standard density.
=#
c = sqrt(T̃ * Rs * Z̃) # (4) but in theory this is also time dependend
ρ̃ = p̃ / (Rs * T̃ * Z̃) # (4)
nothing #hide

#=
The equivalent "pressure capacity" at the nodes is calculated as a sum of the connected
pipe parameters according ot (28).
=#
C₂ = L₁₂*A/(2*ρ̃*c^2) + L₂₃*A/(2*ρ̃*c^2) # (28)
C₃ = L₁₃*A/(2*ρ̃*c^2) + L₂₃*A/(2*ρ̃*c^2) # (28)
nothing #hide

#=
## Load Profile

The paper specifies the load profile at two nodes. We use the package `LinearInterpolations`
to get a callable object which represents this picewise linear interpolation.

However, this function is not Symbolics.jl compatible, so we need to stop Symbolics.jl/ModelingToolkit.jl
from tracing it. To do so, we use `@register_symbolic` to declare it as a symbolic function which is treated
as a blackbox.

Additionally, we need to tell ModelingToolkit about the units of this object. This is just used
for the static unit check during construction of the model. Later one, when we generate the
Julia code from the symbolic reepresentation all units will be stripped.

!!! note "Discontinuities in RHS"
    The picewise linear interpolated function creates discontinuities in the RHS of the system.
    However since we know the times exactly, we can handle this by simply giving a list of explicit
    tstops to the solve command, to make sure those are hit exactly.
=#
load2(t) = -Interpolate(SA[0, 4, 12, 20, 24]*3600, SA[20, 30, 10, 30, 20], extrapolate=LinearInterpolations.Constant(20))(t)
load3(t) = -Interpolate(SA[0, 4, 12, 20, 24]*3600, SA[40, 50, 30, 50, 40], extrapolate=LinearInterpolations.Constant(40))(t)
@register_symbolic load2(t)
@register_symbolic load3(t)
ModelingToolkit.get_unit(op::typeof(load2), _) = u"m^3/s"
ModelingToolkit.get_unit(op::typeof(load3), _) = u"m^3/s"
nothing #hide

#=
## Building the Network

To bild the Network we first need to define the components. This is a two step process:

- first create the symbolic `ODESystem` using ModelingToolkit
- secondly build a NetworkDynamics component function (`ODEVertex`/`ODEEdge`) based on the symbolic system.

In the first step we can use the keyword arguments to pass "default" values for our parameters and states.
Those values will be automaticially transfered to the metadata of the component function the second step.

The second step requires to define the interface variables, i.e. what are the "input" states of your
component function and what are the "output" states.
For `ODEVertex` the input state is the aggregated flow of all connected pipes. The output state is the pressure of the node.
=#
@named v1_mtk = ConstantPressureNode(p_set=p₁_set)
v1 = ODEVertex(v1_mtk, [:q̃_nw], [:p]; name=:v1, vidx=1)
#

@named v2_mtk = VariablePressureNode(C=C₂, load_profile=load2)
v2 = ODEVertex(v2_mtk, [:q̃_nw], [:p]; name=:v2, vidx=2)
#

@named v3_mtk = VariablePressureNode(C=C₃, load_profile=load3)
v3 = ODEVertex(v3_mtk, [:q̃_nw], [:p]; name=:v3, vidx=3)

#=
For `ODEEdge` we have two inputs: the pressure on both source and destination end.
There is a single output state: the volumetric flow. However we also need to tell
NetworkDynamics about the coupling type. In this case we use `AntiSymmetric`, which
meas that the source end will recieve the same flow, just inverted sign.
=#
@named e12_mtk = Pipe(; L=L₁₂, sinθ=sinθ₁₂, A, ρ̃, D, g, c, γ, η, r)
@named e13_mtk = Pipe(; L=L₁₃, sinθ=sinθ₁₃, A, ρ̃, D, g, c, γ, η, r)
@named e23_mtk = Pipe(; L=L₂₃, sinθ=sinθ₂₃, A, ρ̃, D, g, c, γ, η, r)

e12 = ODEEdge(e12_mtk, [:p_src], [:p_dst], [:q̃], AntiSymmetric(); name=:e12, src=1, dst=2)
e13 = ODEEdge(e13_mtk, [:p_src], [:p_dst], [:q̃], AntiSymmetric(); name=:e13, src=1, dst=3)
e23 = ODEEdge(e23_mtk, [:p_src], [:p_dst], [:q̃], AntiSymmetric(); name=:e23, src=2, dst=3)

#=
To build the network object we just need to pass the vertices and edges to the constructor.

Note that we've used the `vidx` and `src`/`dst` keywords in the constructors to define for
each component to which "part" of the network it belongs.

This means, the constructor can automaticially construct a graph based on those informations and
we don't need to pass it explicitly.
=#

nw = Network([v1, v2, v3], [e12, e13, e23])

#=
As a result, we recive a network with 3 unique types (v2 and v3 are similar but
structurally different, because both functions capure a unique loadprofile function).

## Finding a Steady State

To simulate the systme, we first need to find a steadystate. As a "guess" for that
we create a `NWState` object from the network.
This will allocate flat arrays for states `u` and parameters `p` and fill them with the
default values.
=#
uguess = NWState(nw)

#=
This is not a steadystate of the system however. To find a true steadystate we want to ensure
that the lhs of the system is zero.
We can solve for a steady state numerically by defining a Nonlinear Rootfind problem.

To do so, we need to wrap the Network object in a closure. This closure is important to
overwrite all changes of the solver to `u[1]`, to ensure that it is not altered while
solving.
=#

nwwrap = (du, u, p) -> begin
    u[1] = uflat(uguess)[1]
    nw(du, u, p, 0)
    nothing
end
initprob = NonlinearProblem(nwwrap, uflat(uguess), pflat(uguess))
initsol = solve(initprob)

#=
The solver complains a a bit about rank deficient matrices but is able to find a solution.

We can create a new `NWState` object by wrapping the solution from the nonlinear problem and the
original prameters in a new `NWState` object.
=#
u0 = NWState(nw, initsol.u, uguess.p)

#=
## Solving the ODE

Using this as our initial state we can create the actual `ODEProblem`.
Since the ode allways operates on flat state and aprameter arrays we use `uflat` and `pflat` to extract them.
=#
prob = ODEProblem(nw, uflat(u0), (0.0,24*3600), copy(pflat(u0)))
sol = solve(prob, Tsit5(), tstops=[0,4,12,20,24]*3600)
nothing #hide

#=
## Inspect the Solution

Inspecting the solution is all which is left to do.
=#
xticks = ((0:4:24)*3600, string.(0:4:24)) # its nice to display hours
fig = begin
    _fig = Figure()
    row = 1
    ax = Axis(_fig[row, 1]; xlabel="time [h]", ylabel="pressure [Pa]", title="Pressure at nodes", xticks)
    xlims!(ax, sol.t[begin], sol.t[end])
    for i in 1:3
        lines!(ax, sol, idxs=vidxs(nw, i, :p); label="v$i", color=Cycled(i))
    end
    axislegend(ax)
    row += 1

    ax = Axis(_fig[row, 1]; xlabel="time [h]", ylabel="flow [m³/s]", title="Flow through pipes", xticks)
    xlims!(ax, sol.t[begin], sol.t[end])
    for i in 1:2
        lines!(ax, sol, idxs=eidxs(nw, i, :q̃); label="e$i flow", color=Cycled(i))
    end
    axislegend(ax)
    row += 1
    _fig
end


#=
Notably, the "internal" states defined in the symbolic models are not "states" in the sense of the ODE.
For example, we captured the load profile in the `q̃_inj` state of the `VariablePressureNode`.
The only dynamic state of the model however is `p`.
Using the "observables" mechanism from SciML, which is implemented by NetworkDynamics, we can reconstruct
those "optimized" states which have been removed symbolicially.
Here we plot the reconstructed load profile of nodes 2 and 3. Also, we know that node 1 is infinetly stiff,
acting as an infinite source of volumetric flow. We can reconstruct this flow too.
=#
fig = begin
    _fig = Figure()
    row = 1
    ax = Axis(_fig[row, 1]; xlabel="time [h]", ylabel="flow [m³/s]", title="Flow at nodes", xticks)
    xlims!(ax, sol.t[begin], sol.t[end])
    lines!(ax, sol, idxs=vidxs(nw, 1, :q̃_inj); label="v1 compensation", color=Cycled(1))
    for i in 2:3
        lines!(ax, sol, idxs=vidxs(nw, i, :q̃_inj); label="v$i load profile", color=Cycled(i))
    end
    axislegend(ax, position=:rc)
    _fig
end

#=
Lastly we want to observe two internal states of the pipes: the Reynolds number and the mean pressure.
We see, that we're purely in the turbulent flow regime.
=#
fig = begin
    _fig = Figure()
    row = 1
    ax = Axis(_fig[row, 1]; xlabel="time [h]", ylabel="Reynolds number", title="Reynolds number", xticks)
    xlims!(ax, sol.t[begin], sol.t[end])
    for i in 1:3
        lines!(ax, sol, idxs=eidxs(nw, i, :Re); label="e $i", color=Cycled(i))
    end
    hlines!(ax, 2300, color=:black, linestyle=:dash, label="L/T transition")
    axislegend(ax, position=:rb)
    row += 1

    ax = Axis(_fig[row, 1]; xlabel="time [h]", ylabel="Mean pressure [Pa]", title="Mean pressure in pipes", xticks)
    xlims!(ax, sol.t[begin], sol.t[end])
    for i in 1:3
        lines!(ax, sol, idxs=eidxs(nw, i, :pM); label="e $i", color=Cycled(i))
    end
    axislegend(ax, position=:rb)
    _fig
end