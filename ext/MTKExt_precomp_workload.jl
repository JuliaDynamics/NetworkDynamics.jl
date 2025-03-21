using NetworkDynamics
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as Dt

@connector Terminal begin
    u_r(t)
    u_i(t)
    i_r(t), [connect=Flow]
    i_i(t), [connect=Flow]
end
@mtkmodel Swing begin
    @components begin
        terminal = Terminal()
    end
    @variables begin
        ω(t), [description="Rotor frequency", guess=0]
        θ(t), [description="Rotor angle", guess=0]
        Pel(t)
    end
    @parameters begin
        M=0.005
        D=0.1
        V, [guess=1]
        ω_ref=0
        Pm, [description="Mechanical Power", guess=-1]
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
        u_r(t)=1.0
        u_i(t)=0.0
        i_r(t)=1.0
        i_i(t)=0.0
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
vm = VertexModel(swingbus, [:busbar₊i_r, :busbar₊i_i], [:busbar₊u_r, :busbar₊u_i]; verbose=false)
NetworkDynamics.initialize_component!(vm; verbose=false)
