is_precompiling() = ccall(:jl_generating_output, Cint, ()) == 1

using NetworkDynamics
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as Dt

@connector Terminal begin
    u_r(t), [description="d-voltage"]
    u_i(t), [description="q-voltage"]
    i_r(t), [description="d-current", connect=Flow]
    i_i(t), [description="q-current", connect=Flow]
end
@mtkmodel Swing begin
    @components begin
        terminal = Terminal()
    end
    @variables begin
        ω(t)=0.0, [description="Rotor frequency"]
        θ(t)=0.0, [description="Rotor angle"]
        Pel(t), [description="Electrical Power injected into the grid"]
    end
    @parameters begin
        M=1, [description="Inertia"]
        D=0.1, [description="Damping"]
        V=1.0, [description="Voltage magnitude"]
        ω_ref=0, [description="Reference frequency"]
        Pm, [description="Mechanical Power"]
    end
    @equations begin
        Dt(θ) ~ ω - ω_ref
        Dt(ω) ~ 1/M * (Pm - D*ω - Pel)
        Pel ~ terminal.u_r*terminal.i_r + terminal.u_i*terminal.i_i
        terminal.u_r ~ V*cos(θ)
        terminal.u_i ~ V*sin(θ)
    end
end
@mtkmodel BusBase begin
    @variables begin
        u_r(t)=1, [description="bus d-voltage", output=true]
        u_i(t)=0, [description="bus q-voltage", output=true]
        i_r(t), [description="bus d-current (flowing into bus)", input=true]
        i_i(t), [description="bus d-current (flowing into bus)", input=true]
        P(t), [description="bus active power (flowing into network)"]
        Q(t), [description="bus reactive power (flowing into network)"]
        u_mag(t), [description="bus voltage magnitude"]
        u_arg(t), [description="bus voltage argument"]
        i_mag(t), [description="bus current magnitude"]
        i_arg(t), [description="bus current argument"]
        ω(t), [description="bus angular frequency"]
    end
    @equations begin
        P ~ u_r * (-i_r) + u_i * (-i_i)
        Q ~ u_i * (-i_r) - u_r * (-i_i)
        u_mag ~ sqrt(u_r^2 + u_i^2)
        u_arg ~ atan(u_i, u_r)
        i_mag ~ sqrt(i_r^2 + i_i^2)
        i_arg ~ atan(i_i, i_r)
        ω ~ Dt(u_arg)
    end
end
@mtkmodel BusBar begin
    @extend BusBase()
    @components begin
        terminal = Terminal()
    end
    @equations begin
        u_r ~ terminal.u_r
        u_i ~ terminal.u_i
        i_r ~ terminal.i_r
        i_i ~ terminal.i_i
    end
end
@mtkmodel Bus begin
    @components begin
        swing = Swing()
        busbar = BusBar()
    end
    @equations begin
        connect(swing.terminal, busbar.terminal)
    end
end
@named swingbus = Bus()

if is_precompiling()
    VertexFunction(swingbus, [:busbar₊i_r, :busbar₊i_i], [:busbar₊u_r, :busbar₊u_i]; verbose=false)
else
    @info "ODEVertex"
    @time @eval VertexFunction(swingbus, [:busbar₊i_r, :busbar₊i_i], [:busbar₊u_r, :busbar₊u_i]; verbose=false)
end
