# External Inputs
External inputs for components are way to pass signals between components outside of the network structure.
The most common usecase for that are control systems: make your vertex `i` depend on some state of vertex `j`.

If used, this will essentially wides the number of received inputs of a component function. I.e. the baseline mathematical models for vertex models

```math
\begin{aligned}
M^{\mathrm v}\,\frac{\mathrm{d}}{\mathrm{d}t}x^{\mathrm v} &= f^{\mathrm v}(u^{\mathrm v}, i^{\mathrm v}, i^{\mathrm{ext}}, p^{\mathrm v}, t)\\
y^{\mathrm v} &= g^{\mathrm v}(u^{\mathrm v}, i^{\mathrm v}, i^{\mathrm{ext}}, p^{\mathrm v}, t)
\end{aligned}
```
```julia
fᵥ(dxᵥ, xᵥ, e_aggr, ext, pᵥ, t)
gᵥ(yᵥ, xᵥ, [e_aggr, ext,] pᵥ, t)
```
and edge models
```math
\begin{aligned}
M^{\mathrm e}\,\frac{\mathrm{d}}{\mathrm{d}t}x^{\mathrm e} &= f^{\mathrm e}(u^{\mathrm e}, y^{\mathrm v}_{\mathrm{src}}, y^{\mathrm v}_{\mathrm{dst}}, i^{\mathrm{ext}}, p^{\mathrm e}, t)\\
y^{\mathrm e}_{\mathrm{dst}} &= g_\mathrm{dst}^{\mathrm e}(u^{\mathrm e}, y^{\mathrm v}_{\mathrm{src}}, y^{\mathrm v}_{\mathrm{dst}}, i^{\mathrm{ext}}, p^{\mathrm e}, t)\\
y^{\mathrm e}_{\mathrm{src}} &= g_\mathrm{src}^{\mathrm e}(u^{\mathrm e}, y^{\mathrm v}_{\mathrm{src}}, y^{\mathrm v}_{\mathrm{dst}}, i^{\mathrm{ext}}, p^{\mathrm e}, t)
\end{aligned}
```
```julia
fₑ(dxₑ, xₑ, v_src, v_dst, ext, pₑ, t)
gₑ(y_src, y_dst, xᵥ, [v_src, v_dst, ext,] pₑ, t)
```
change. You may still oomit the input section from `g` according to the different [Feed Forward Behavior](@ref)s. However you either have to use *all* inputs (including `ext`) or none.

## Usage
Vertex and Edge models with external inputs can be created using the `extin` keyword of the [`EdgeModel`](@ref) and [`VertexModel`](@ref) constructors.

You need to pass a vector of `SymbolicIndices` ([`VIndex`](@ref) and [`EIndex`](@ref)), e.g. 
```julia
VertexModel(... , extin=[VIndex(12,:x), VIndex(9, :x)], ...)
``` 
means your vertex receives a 2 dimensional external input vector with the states `x` of vertices 12 and 9.
See below for hands on example.

## Limitations
There are some limitations in place. You can only reference **states** (i.e. things that appear in `xᵥ` or `xₑ` of some component model) or **outputs of non-feed-forward components**, i.e. states `yᵥ` or `yₑ` of some component model which does not have feed forward behavior in their `g` function.

## Example
As an example system, we'll consider two capacitors with a resistor between them.
Vertex 1 `v1` has a controllable current source. 
Using a PI controller for the current source, it tries to keep the voltage at the 
second vertex stable under the disturbance of some time periodic current sind at `v2`.

```

                  v1    Resistor   v2
PI controlled   ─→─o─←────MMM────→─o─→─ time dependent 
current source     ┴               ┴    current sink
                   ┬               ┬
                   ⏚               ⏚
```

The example will be implemented 2 times: in plain NetworkDynamcics and using MTK.

### Plain NetworkDynamics

First we define the resistor:
```@example extin
using NetworkDynamics
using OrdinaryDiffEqTsit5
using CairoMakie

function resistor_g(idst, vsrc, vdst, (R,), t)
    idst[1] = (vsrc[1] - vdst[1])/R
end
resistor = EdgeModel(g=AntiSymmetric(resistor_g), 
                     outsym=:i, insym=:v, psym=:R=>0.1, 
                     src=:source, dst=:load, name=:resistor)
```

Then we define the "load" vertex with sinusoidial load profile:
```@example extin
function f_load(dv, v, iin, (C,), t)
    dv[1] = 1/C*(iin[1] - (1 + 0.1*sin(t)))
end
load = VertexModel(f=f_load, g=1, 
                   sym=:V=>0.5, insym=:i, psym=:C=>1,
                   vidx=2, name=:load)
```
Lastly, we define the "source" vertex
```@example extin
function f_source(dv, v, iin, extin, (C,Vref,Ki,Kp), t)
    Δ = Vref - extin[1]      # tracking error
    dv[2] = Δ                # integrator state of PI
    i_inj = Kp*Δ + Ki*v[2]   # controller output
    dv[1] = 1/C*(iin[1] + i_inj)
end
source = VertexModel(f=f_source, g=1, 
                     sym=[:V=>0.5,:ξ=>0], insym=:i, 
                     psym=[:C=>1, :Vref=>1, :Ki=>0.5, :Kp=>10],
                     extin=[VIndex(:load, :V)], 
                     vidx=1, name=:source)
```
Then we can create the network and simulate:
```@example extin
nw = Network([source, load], [resistor])
u0 = NWState(nw) # everything has default values
prob = ODEProblem(nw, uflat(u0), (0,100), pflat(u0))
sol = solve(prob, Tsit5())

fig, ax, p = lines(sol, idxs=VIndex(:load, :V), label="Voltage @ load");
lines!(ax, sol, idxs=VPIndex(:source, :Vref), label="Reference", color=Cycled(2));
axislegend(ax; position=:rb);
fig # hide
```

### MTK Models
First we define the resistor:
```@example extin_mtk
using NetworkDynamics
using OrdinaryDiffEqTsit5
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as Dt
using CairoMakie

@mtkmodel Resistor begin
    @variables begin
        i(t), [description="Current at dst end"]
        V_src(t), [description="Voltage at src end"]
        V_dst(t), [description="Voltage at dst end"]
    end
    @parameters begin
        R=0.1, [description="Resistance"]
    end
    @equations begin
        i ~ (V_src - V_dst)/R
    end
end
@named resistor = Resistor()
resistor_edge = EdgeModel(resistor, [:V_src], [:V_dst], AntiSymmetric([:i]); src=:load, dst=:source)
```

Then we define the "load" vertex with sinusoidial load profile:
```@example extin_mtk
@mtkmodel Load begin
    @variables begin
        V(t)=0.5, [description="Voltage"]
        i_load(t), [description="Load current"]
        i_grid(t), [description="Current from grid"]
    end
    @parameters begin
        C=1, [description="Capacitance"]
    end
    @equations begin
        i_load ~ 1 + 0.1*sin(t)
        Dt(V) ~ 1/C*(i_grid - i_load)
    end
end
@named load = Load()
load_vertex = VertexModel(load, [:i_grid], [:V]; vidx=2)
```
Lastly, we define the "source" vertex
```@example extin_mtk
@mtkmodel Source begin
    @variables begin
        V(t)=0.5, [description="Voltage"]
        ξ(t)=0, [description="Integrator state"]
        i_grid(t), [description="Current from grid"]
        i_source(t), [description="Current from source"]
        Δ(t), [description="Tracking Error"]
        V_load(t), [description="Voltage at load"]
    end
    @parameters begin
        C=1, [description="Capacitance"]
        Vref=1, [description="Reference voltage"]
        Ki=0.5, [description="Integral gain"]
        Kp=10, [description="Proportional gain"]
    end
    @equations begin
        Δ ~ Vref - V_load
        Dt(ξ) ~ Δ
        i_source ~ Kp*Δ + Ki*ξ
        Dt(V) ~ 1/C*(i_grid + i_source)
    end
end
@named source = Source()
source_vertex = VertexModel(source, [:i_grid], [:V]; 
                            extin=[:V_load => VIndex(:load, :V)], vidx=1)
```
Then we can create the network and simulate:
```@example extin_mtk
nw = Network([source_vertex, load_vertex], [resistor_edge])
u0 = NWState(nw) # everything has default values
prob = ODEProblem(nw, uflat(u0), (0,100), pflat(u0))
sol = solve(prob, Tsit5())

let
    fig = Figure();
    ax1 = Axis(fig[1,1]);
    lines!(ax1, sol, idxs=VIndex(:load, :V), label="Voltage @ load");
    lines!(ax1, sol, idxs=VPIndex(:source, :Vref), label="Reference", color=Cycled(2));
    axislegend(ax1; position=:rb);
    ax2 = Axis(fig[2,1]);
    lines!(ax2, sol, idxs=VIndex(:load, :i_load), label="load current");
    lines!(ax2, sol, idxs=VIndex(:source, :i_source), label="source current", color=Cycled(2));
    axislegend(ax2);
    fig
end
```
Using MTK for modeling, we can also inspect the currents `i_load` and `i_source` the MTK interface preserves the [Observables](@ref).




