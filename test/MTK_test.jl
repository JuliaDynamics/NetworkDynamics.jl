import ModelingToolkit as MTK
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as Dt
using NetworkDynamics
using OrdinaryDiffEqTsit5
using LinearAlgebra
using Graphs
using Chairmarks: @b
using Test

@mtkmodel Bus begin
    @variables begin
        θ(t), [description = "voltage angle", output=true]
        P(t), [description = "Electical Powerflow into Network", input=true]
    end
end;

@mtkmodel SwingNode begin
    @extend Bus()
    @variables begin
        ω(t) = 0.0, [description = "Rotor frequency"]
    end
    @parameters begin
        M = 1, [description = "Inertia"]
        D = 0.1, [description = "Damping"]
        Pmech, [description = "Mechanical Power"]
    end
    @equations begin
        Dt(θ) ~ ω
        Dt(ω) ~ 1/M * (Pmech - D*ω + P)
    end
end;

@named swing = SwingNode()
v = VertexModel(swing, [:P], [:θ])

data = NetworkDynamics.rand_inputs_fg(v)
b = @b $(NetworkDynamics.compfg(v))($data...)
@test b.allocs == 0

@mtkmodel Line begin
    @variables begin
        srcθ(t), [description = "voltage angle at src end", input=true]
        srcP(t), [description = "flow towards not at src end", output=true]
        dstθ(t), [description = "voltage angle at dst end", input=true]
        dstP(t), [description = "flow towards not at dst end", output=true]
    end
end
@mtkmodel StaticPowerLine begin
    @extend Line()
    @parameters begin
        K = 100.0, [description = "Line conductance"]
    end
    @variables begin
       Δθ(t)
    end
    @equations begin
        Δθ ~ dstθ - srcθ
        srcP ~ -K*sin(Δθ)
        dstP ~ -srcP
    end
end
@named line = StaticPowerLine()
e = EdgeModel(line, [:srcθ], [:dstθ], [:srcP], [:dstP])
@test NetworkDynamics.insym(e).src == [:srcθ]
@test NetworkDynamics.insym(e).dst == [:dstθ]

data = NetworkDynamics.rand_inputs_fg(e)
b = @b $(NetworkDynamics.compfg(e))($data...)
@test b.allocs == 0

g = complete_graph(4)
nw = Network(g, v, e)
u0 = NWState(nw)
u0.v[:,:θ] .= 0
u0.p.v[1, :Pmech] = 1.0
u0.p.v[2, :Pmech] = 1.0
u0.p.v[3, :Pmech] = -1.0
u0.p.v[4, :Pmech] = -1.0

prob = ODEProblem(nw, uflat(u0), (0, 10), pflat(u0))
sol = solve(prob, Tsit5())

# plot(sol; idxs=vidxs(nw, :, :θ))
# plot(sol; idxs=vidxs(nw, :, :ω))

####
#### dq model
####
@mtkmodel DQBus begin
    @variables begin
        u_r(t), [description = "d-voltage", output=true]
        u_i(t), [description = "q-voltage", output=true]
        i_r(t), [description = "d-current", input=true]
        i_i(t), [description = "d-current", input=true]
    end
end

rotm(θ) = [cos(θ) -sin(θ); sin(θ) cos(θ)]


@mtkmodel DQSwing begin
    @extend DQBus()
    @variables begin
        ω(t) = 0.0, [description = "Rotor frequency"]
        θ(t) = 0.0, [description = "Rotor angle"]
        Pel(t), [description = "Electrical Power"]
    end
    @parameters begin
        M = 1, [description = "Inertia"]
        D = 0.1, [description = "Damping"]
        Pmech, [description = "Mechanical Power"]
        V = 1.0, [description = "Voltage magnitude"]
    end
    @equations begin
        Dt(θ) ~ ω
        Pel ~ real((u_r + im*u_i) * (i_r - im*i_i))
        Dt(ω) ~ 1/M * (Pmech - D*ω + Pel)
        u_r ~ V*cos(θ)
        u_i ~ V*sin(θ)
    end
end

@named dqswing = DQSwing()
v = VertexModel(dqswing, [:i_r, :i_i], [:u_r, :u_i])
@test v.mass_matrix == Diagonal([1,1])

@mtkmodel DQLine begin
    @variables begin
        src_u_r(t), [description = "src d-voltage", input=false]
        src_u_i(t), [description = "src q-voltage", input=false]
        dst_u_r(t), [description = "dst d-voltage", input=false]
        dst_u_i(t), [description = "dst q-voltage", input=false]
        src_i_r(t), [description = "src d-current", output=true]
        src_i_i(t), [description = "src d-current", output=true]
        dst_i_r(t), [description = "dst d-current", output=true]
        dst_i_i(t), [description = "dst d-current", output=true]
    end
end

@mtkmodel DQYLine begin
    @extend DQLine()
    @parameters begin
        Y_r, [description = "Real part of admittance"]
        Y_i, [description = "Imaginary part of admittance"]
    end
    @equations begin
        icomplex ~ (Y_r + im*Y_i) * (dst_u_r + im*dst_u_i - src_u_r - im*src_u_i)
        dst_i_r ~ real(icomplex)
        dst_i_i ~ imag(icomplex)
    end
end

@mtkmodel DQPiLine begin
    @extend DQLine()
    @parameters begin
        R = 1.0, [description = "Resistance"]
        X = 1.0, [description = "Reactance"]
        src_B = 0.0, [description = "Shunt susceptance at src end"]
        dst_B = 0.0, [description = "Shunt susceptance at dst end"]
    end
    @equations begin
        src_i_r ~ -real(-(src_u_r + im*src_u_i - dst_u_r - im*dst_u_i)/(R + im*X) - (im*src_B)*(src_u_r + im*src_u_i))
        src_i_i ~ -imag(-(src_u_r + im*src_u_i - dst_u_r - im*dst_u_i)/(R + im*X) - (im*src_B)*(src_u_r + im*src_u_i))
        dst_i_r ~ -real(+(src_u_r + im*src_u_i - dst_u_r - im*dst_u_i)/(R + im*X) - (im*dst_B)*(dst_u_r + im*dst_u_i))
        dst_i_i ~ -imag(+(src_u_r + im*src_u_i - dst_u_r - im*dst_u_i)/(R + im*X) - (im*dst_B)*(dst_u_r + im*dst_u_i))
    end
end

@named piline = DQPiLine()
l = EdgeModel(piline, [:src_u_r, :src_u_i], [:dst_u_r, :dst_u_i], [:src_i_r, :src_i_i], [:dst_i_r, :dst_i_i])

## test nested model
@mtkmodel Shaft begin
    @variables begin
        θ(t), [description = "Shaft angle"]
        ω(t), [description = "Shaft speed"]
        Pmech(t), [description = "Mechanical Power"]
        Pel(t), [description = "Electrical Power"]
    end
    @parameters begin
        M = 1, [description = "Inertia"]
        D = 0.1, [description = "Friction Damping"]
    end
    @equations begin
        Dt(θ) ~ ω
        Dt(ω) ~ 1/M * (Pmech - D*ω + Pel)
    end
end

@mtkmodel Gov begin
    @variables begin
        Pmech(t), [description = "Mechanical Power"]
        ω_meas(t),[description = "Measured Rotor Frequency"]
    end
    @parameters begin
        D = 0.1, [description = "Governor Droop"]
        Pref = 1.0, [description = "Reference Power"]
    end
    @equations begin
        Pmech ~ Pref - D*ω_meas
    end
end

@mtkmodel NestedSwing begin
    @components begin
        shaft = Shaft()
        gov = Gov()
    end
    @variables begin
        u_r(t), [description = "d-voltage", output=true]
        u_i(t), [description = "q-voltage", output=true]
        i_r(t), [description = "d-current", input=true]
        i_i(t), [description = "d-current", input=true]
    end
    @parameters begin
        V = 1.0, [description = "Voltage magnitude"]
    end
    @equations begin
        shaft.Pel ~ real((u_r + im*u_i) * (i_r - im*i_i))
        gov.ω_meas ~ shaft.ω
        gov.Pmech ~ shaft.Pmech
        u_r ~ V*cos(shaft.θ)
        u_i ~ V*sin(shaft.θ)
    end
end

@named nestedswing = NestedSwing()
v = VertexModel(nestedswing, [:i_r, :i_i], [:u_r, :u_i])
@test v.mass_matrix == Diagonal([1,1])
data = NetworkDynamics.rand_inputs_fg(v)
b = @b $(NetworkDynamics.compfg(v))($data...)
@test b.allocs == 0


@testset "Test constants in MTK models" begin
    # from MTK@10 onwars, constants are just non-tunable parameters!
    @mtkmodel DQSwing_Constants begin
        @extend DQBus()
        @variables begin
            ω(t) = 0.0, [description = "Rotor frequency"]
            θ(t) = 0.0, [description = "Rotor angle"]
            Pel(t), [description = "Electrical Power"]
        end
        @constants begin
            V = 1, [description = "Voltage magnitude"]
            useless = 0
        end
        @parameters begin
            M = 1, [description = "Inertia"]
            D = 0.1, [description = "Damping"]
            Pmech, [description = "Mechanical Power"]
        end
        @equations begin
            Dt(θ) ~ ω
            Dt(ω) ~ 1/M * (Pmech - D*ω + Pel)
            Pel ~ real((u_r + im*u_i) * (i_r - im*i_i)) + useless
            u_r ~ V*cos(θ)
            u_i ~ V*sin(θ)
        end
    end
    @named withconst = DQSwing_Constants()
    v = VertexModel(withconst, [:i_r, :i_i], [:u_r, :u_i])

    data = NetworkDynamics.rand_inputs_fg(v)
    NetworkDynamics.compfg(v)(data...) # no error

    data = NetworkDynamics.rand_inputs_obsf(v)
    v.obsf(data...) # no error
end

# Recreate OpPoDyn's Terminal connector
@connector Terminal begin
    u_r(t), [description="d-voltage"]
    u_i(t), [description="q-voltage"]
    i_r(t), [guess=0, description="d-current", connect=Flow]
    i_i(t), [guess=0, description="q-current", connect=Flow]
end

# Recreate OpPoDyn's PVConstraint model
@mtkmodel PVConstraint begin
    @components begin
        terminal = Terminal()
    end
    @parameters begin
        P
        V
    end
    @equations begin
        P ~ terminal.u_r*terminal.i_r + terminal.u_i*terminal.i_i
        V^2 ~ terminal.u_r^2 + terminal.u_i^2
    end
end

# Recreate OpPoDyn's BusBase model
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
    end
    @equations begin
        # Observed equations - flipped sign in P and Q, flow direction opposite to i
        P ~ u_r * (-i_r) + u_i * (-i_i)
        Q ~ u_i * (-i_r) - u_r * (-i_i)
        u_mag ~ sqrt(u_r^2 + u_i^2)
        u_arg ~ atan(u_i, u_r)
        i_mag ~ sqrt(i_r^2 + i_i^2)
        i_arg ~ atan(i_i, i_r)
    end
end

@testset "Test transitive output->internal->state alias" begin
    # Recreate OpPoDyn's BusBar model
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

    @named pv = PVConstraint(; P=2, V=1)
    @named busbar = BusBar()
    mtkbus = System(connect(busbar.terminal, pv.terminal), t; systems=[busbar, pv], name=:pvbus)
    vf = VertexModel(mtkbus, [:busbar₊i_r, :busbar₊i_i], [:busbar₊u_r, :busbar₊u_i], verbose=false)

    @test vf.sym == [:busbar₊u_r,:busbar₊u_i]
end

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
@mtkmodel StaticProsumerNode begin
    @extend GasNode()
    @parameters begin
        q̃_prosumer, [description="flow injected by prosumer"]
    end
    @equations begin
        -q̃_nw ~ q̃_prosumer
    end
end
@mtkmodel Wrapper begin
    @components begin
        prosumer = StaticProsumerNode()
    end
    @variables begin
        p(t), [description="Pressure at prosumer"]
        q̃_nw(t), [description="aggregated flow from pipes into node"] # node input
    end
    @equations begin
        p ~ prosumer.p # connect the pressure output of the prosumer to the wrapper
        q̃_nw ~ prosumer.q̃_nw # connect the flow input of the prosumer to the wrapper
    end
end
@mtkmodel WrapperFixed begin
    @components begin
        prosumer = StaticProsumerNode()
    end
    @variables begin
        p(t), [description="Pressure at prosumer"]
        q̃_nw(t), [description="aggregated flow from pipes into node"] # node input
    end
    @equations begin
        p ~ prosumer.p # connect the pressure output of the prosumer to the wrapper
        q̃_nw + implicit_output(p) ~ prosumer.q̃_nw # connect the flow input of the prosumer to the wrapper
    end
end
@testset "Test transformation of implicit outputs" begin
    # first, lets test if the underlying MTK problem still exists
    # see https://github.com/SciML/ModelingToolkit.jl/pull/3686
    @mtkmodel ImplicitForcing begin
        @variables begin
            u(t), [description = "Input Variable", input=true]
            y(t), [description = "fully implicit output", output=true]
        end
        @equations begin
            0 ~ sqrt(u) # implicitly forces output y because u=f(y) in  closed loop
        end
    end
    @named implicit = ImplicitForcing()
    simp = mtkcompile(implicit; inputs = ModelingToolkit.unbound_inputs(implicit))
    @test isempty(equations(simp)) # the equation was dropped!

    # test fully implicit outputs
    @mtkmodel FullyImplicit begin
        @variables begin
            u(t), [description = "Input Variable", input=true]
            x(t), [description = "Explicit Variable"]
            y(t), [description = "Implicit Variable, present in equations"]
            z(t), [description = "fully implicit variable, not present but output", output=true]
        end
        @equations begin
            Dt(x) ~ -x
            0 ~ sqrt(y+x)
            0 ~ u # implicitly forces z becaus u(z)
        end
    end
    @named fullyimplicit = FullyImplicit()
    @test_throws ArgumentError VertexModel(fullyimplicit, [:u], [:z])

    # from init_tutorial
    # dependent MTKModels need to be defined at top level, so they are in front of the testset
    @named prosumer = StaticProsumerNode() # consumer
    @test_throws ArgumentError VertexModel(prosumer, [:q̃_nw], [:p])
    @named prosumer_wrapped = Wrapper()
    @test_throws ArgumentError VertexModel(prosumer_wrapped, [:q̃_nw], [:p])
    @named prosumer_fixed = WrapperFixed()
    VertexModel(prosumer_fixed, [:q̃_nw], [:p]) # no throw
end
