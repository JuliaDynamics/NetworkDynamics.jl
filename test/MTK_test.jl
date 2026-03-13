using ModelingToolkitBase
using ModelingToolkitBase: t_nounits as t, D_nounits as Dt
using NetworkDynamics
using OrdinaryDiffEqTsit5
using LinearAlgebra
using Graphs
using Chairmarks: @b
using Test
using SciCompDSL
mtkext = Base.get_extension(NetworkDynamics, :NetworkDynamicsMTKExt)

@testset "get_variables_deriv test" begin
    @variables x(t) y(t)
    # if this test gets fixed, we can get rid of workaround get_variables_deriv
    @test get_variables(Dt(x) + y) == Set([Dt(x), y])
    @test mtkext.get_variables_deriv(Dt(x) + y) == Set([x, y])
    @test mtkext.get_variables_deriv(Dt(x) + y) == Set([x, y])
end

@testset "eqtype test" begin
    @variables begin
        x(t)
        y(t)
        z(t)
        a
        b
        c
    end
    eq = Dt(x) ~ y + a
    @test mtkext.eq_type(eq) == (:explicit_diffeq, x.val)

    eq = y ~ x + b
    @test mtkext.eq_type(eq) == (:explicit_algebraic, y.val)

    eq = y ~ x + b + y
    @test mtkext.eq_type(eq) == (:implicit_algebraic, y.val)

    eq = 0 ~ x+y+b
    @test mtkext.eq_type(eq) == (:implicit_algebraic, nothing)

    # non zero on the lhs is not expected
    eq = 1 ~ x+y+b
    @test_throws "non-zero constant on the lhs" mtkext.eq_type(eq)

    eq = Dt(x) ~ Dt(x)
    @test_throws "differentials on the rhs" mtkext.eq_type(eq)

    eq = 0 ~ Dt(x)
    @test_throws "differentials on the rhs" mtkext.eq_type(eq)

    # paraemter on the lhs is not expected
    eq = a ~ x+y
    @test_throws "Can't determine eq type" mtkext.eq_type(eq)
end

@testset "get_scaled_diff test" begin
    @variables x(t) y(t) z

    # scaled by symbolic variable
    r = mtkext.get_scaled_diff((z*Dt(x)).val)
    @test !isnothing(r)
    @test isequal(r[1], Dt(x).val)
    @test isequal(r[2], z.val)

    # scaled by literal number
    r = mtkext.get_scaled_diff((2*Dt(x)).val)
    @test !isnothing(r)
    @test isequal(r[1], Dt(x).val)
    @test isequal(r[2], 2)

    # scaled by symbolic power
    r = mtkext.get_scaled_diff((z^2*Dt(x)).val)
    @test !isnothing(r)
    @test isequal(r[1], Dt(x).val)
    @test isequal(r[2], z.val^2)

    # negative coefficient
    r = mtkext.get_scaled_diff((-z*Dt(x)).val)
    @test !isnothing(r)
    @test isequal(r[1], Dt(x).val)
    @test isequal(r[2], -z.val)

    # plain differential — not a scaled diff
    @test isnothing(mtkext.get_scaled_diff(Dt(x).val))

    # squared differential — not a scaled diff
    @test isnothing(mtkext.get_scaled_diff((z*Dt(x)^2).val))

    # non-differential product — not a scaled diff
    @test isnothing(mtkext.get_scaled_diff((z*x).val))
end

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
        M = 1, [guess=0.1, description = "Inertia"]
        D = 0.1, [guess=M, description = "Damping"]
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
    vf = VertexModel(mtkbus, [:busbar₊i_r, :busbar₊i_i], [:busbar₊u_r, :busbar₊u_i], verbose=true)

    @test Set(vf.sym) == Set([:busbar₊u_r,:busbar₊u_i])
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
    @test_throws r"outputs .* do not appear in the equations" VertexModel(fullyimplicit, [:u], [:z])
    # works when we assume io coupling
    VertexModel(fullyimplicit, [:u], [:z]; assume_io_coupling=true)
    # from init_tutorial
    # dependent MTKModels need to be defined at top level, so they are in front of the testset
    @named prosumer = StaticProsumerNode() # consumer
    @test_throws r"outputs .* do not appear in the equations" VertexModel(prosumer, [:q̃_nw], [:p])
    @named prosumer_wrapped = Wrapper()
    # the below command used to fail, but works now with custom matching/tearing in NetworkDynamics
    VertexModel(prosumer_wrapped, [:q̃_nw], [:p]; verbose=false)
    # the fixed version works too
    @named prosumer_fixed = WrapperFixed()
    VertexModel(prosumer_fixed, [:q̃_nw], [:p])
end

@testset "Error on vector variables" begin
    @mtkmodel VectorModel begin
        @variables begin
            # out(t)[1:2]
            out(t)[1:2]
        end
        @parameters begin
            A[1:2], [description = "vector parameter"]
        end
        @equations begin
            out[1] ~ A[1] + A[2]
            out[2] ~ A[1] - A[2]
        end
    end
    @named vectormodel = VectorModel()
    @test_logs (:warn, r"does not support vector-variables .* A\[1\], A\[2\], \(out\(t\)\)\[1\], \(out\(t\)\)\[2\]") begin
        try
            VertexModel(vectormodel, [\], [:out])
        catch
        end
    end
end

####
#### ComponentPostprocesssing tests
#### needs to be on top level because of modules and mtkmacro
####

module ToLate
    struct CustomMetadata end
    using ModelingToolkitBase
    using ModelingToolkitBase: t_nounits as t, D_nounits as Dt
    using NetworkDynamics: ComponentPostprocessing
    using SciCompDSL
    @mtkmodel LateModel begin
        @variables begin
            in(t)
            out(t)
        end
        @equations begin
            Dt(out) ~ in
        end
        @metadata begin
            ComponentPostprocessing = to_late_defined
            # CoponentPostprocessing = :foo
            # CustomMetadata = :foo
        end
    end
    function to_late_defined(cf, namespace)
        # this function is created after the component
    end
end
module InTime
    struct CustomMetadata end
    using ModelingToolkitBase
    using ModelingToolkitBase: t_nounits as t, D_nounits as Dt
    using SciCompDSL
    using NetworkDynamics: ComponentPostprocessing
    cfref = Ref{Any}(nothing)
    nsref = Ref{Any}(nothing)
    function in_time end
    @mtkmodel Model begin
        @variables begin
            in(t)
            out(t)
        end
        @equations begin
            Dt(out) ~ in
        end
        @metadata begin
            ComponentPostprocessing = in_time
        end
    end
    function in_time(cf, namespace)
        cfref[] = cf
        nsref[] = namespace
    end
end
@testset "postrpocessing callback defined after model" begin
    testmodule = @__MODULE__
    mod = testmodule.ToLate.LateModel(; name=:mod)
    try
        VertexModel(mod, [:in], [:out])
        @test false
    catch e
        @test contains(e.msg, "postprocessing function `to_late_defined` included in")
    end

    @info "Befor in time model"
    mod = testmodule.InTime.Model(; name=:mod)
    vm = VertexModel(mod, [:in], [:out])

    @test InTime.cfref[] === vm
    @test InTime.nsref[] == ""
end

postprocessing_called = Any[]
function sub1_postproc(cf, namespace)
    push!(postprocessing_called, (:sub1, cf, namespace))
end
@mtkmodel SubModel1 begin
    @metadata begin
        ComponentPostprocessing = sub1_postproc
    end
end
function sub2_postproc(cf, namespace)
    push!(postprocessing_called, (:sub2, cf, namespace))
end
@mtkmodel SubModel2 begin
    @components begin
        sub1_in_sub2 = SubModel1()
    end
    @metadata begin
        ComponentPostprocessing = sub2_postproc
    end
end
function toplevel_postproc(cf, namespace)
    push!(postprocessing_called, (:toplevel, cf, namespace))
end
@mtkmodel Toplevel begin
    @components begin
        sub1 = SubModel1()
        sub2 = SubModel2()
    end
    @metadata begin
        ComponentPostprocessing = toplevel_postproc
    end
    @variables begin
        in(t)
        out(t)
    end
    @equations begin
        Dt(out) ~ in
    end
end

@testset "test in submodel" begin
    empty!(postprocessing_called)
    @named mod = Toplevel()
    vm = VertexModel(mod, [:in], [:out])
    @test postprocessing_called[1] == (:sub1, vm, "sub1")
    @test postprocessing_called[2] == (:sub1, vm, "sub2₊sub1_in_sub2")
    @test postprocessing_called[3] == (:sub2, vm, "sub2")
    @test postprocessing_called[4] == (:toplevel, vm, "")
end

@testset "test assume_io_coupling=true" begin
    @mtkmodel GridFollowingInverter begin
        @parameters begin
            w_cu
            omega_ini
            Kp_omega
            Ki_omega
            K_omega
            K_v
            V_ref
            omega_ref
            P_ref
            Q_ref
            KpP
            KiP
            KpQ
            KiQ
        end
        @variables begin
            θ(t)
            z_omega(t)
            u_fil_d(t)
            u_fil_q(t)
            z_P(t)
            z_Q(t)
            V_d(t)
            V_q(t)
            ω_meas(t)
            V_fil_mag(t)
            P_rd(t)
            Q_rd(t)
            P(t)
            Q(t)
            i_d(t)
            i_q(t)
            u_r(t)
            u_i(t)
            i_r(t)
            i_i(t)
        end
        begin
            T_to_loc(α)  = [ cos(α)   sin(α);
                            -sin(α)   cos(α) ]
            T_to_glob(α) = T_to_loc(-α)
        end
        @equations begin
            # voltage to dq
            [V_d, V_q] .~ T_to_loc(θ) * [u_r, u_i]

            # filter
            Dt(u_fil_d) ~ w_cu * (V_d - u_fil_d)
            Dt(u_fil_q) ~ w_cu * (V_q - u_fil_q)

            # PLL
            Dt(z_omega) ~ V_q
            ω_meas ~ omega_ini + Kp_omega * V_q + Ki_omega * z_omega

            # voltage magnitude
            V_fil_mag ~ sqrt(u_fil_d^2 + u_fil_q^2)

            # droop
            P_rd ~ -K_omega * (ω_meas - omega_ref) + P_ref
            Q_rd ~ -K_v     * (V_fil_mag - V_ref)   + Q_ref

            # PI controllers for active and reactive power
            Dt(z_P) ~ (P_rd - P)
            Dt(z_Q) ~ (Q_rd - Q)

            P ~ KpP * (P_rd - P) + KiP * z_P
            Q ~ KpQ * (Q_rd - Q) + KiQ * z_Q

            # dq currents
            i_d ~ (P * u_fil_d + Q * u_fil_q) / (u_fil_d^2 + u_fil_q^2)
            i_q ~ (P * u_fil_q - Q * u_fil_d) / (u_fil_d^2 + u_fil_q^2)

            # to global frame and connect to terminal
            [i_r, i_i] .~ T_to_glob(θ) * [i_d, i_q]

            # angle integration
            Dt(θ) ~ ω_meas
        end
    end

    @named inverter = GridFollowingInverter()

    # Without assume_io_coupling, this should throw RHSDifferentialsError
    # because i_r/i_i outputs depend on θ which has a derivative
    # NOT neede anymore sice custom linear alg reduction does not attempt to remove this
    # @test_throws NetworkDynamics.RHSDifferentialsError begin
    #     vm = VertexModel(inverter, [:i_r, :i_i], [:u_r, :u_i]; verbose=true)
    #     vm.metadata[:equations]
    # end
    VertexModel(inverter, [:i_r, :i_i], [:u_r, :u_i]; verbose=false, assume_io_coupling=false)

    # With assume_io_coupling=true, the construction should succeed
    # This forces MTK to recognize the input->output coupling even with derivatives
    v = VertexModel(inverter, [:i_r, :i_i], [:u_r, :u_i];
                    verbose=false, assume_io_coupling=true)
    data = NetworkDynamics.rand_inputs_fg(v)
    NetworkDynamics.compfg(v)(data...)
end

function topologicical_sorted(eqs)
    lhs = [eq.lhs for eq in eqs]
    for i in eachindex(eqs)
        not_allowed = Set(lhs[i:end])
        rhs_syms = ModelingToolkitBase.get_variables(eqs[i].rhs)
        if !isempty(intersect(not_allowed, rhs_syms))
            return false
        end
    end
    return true
end

@testset "Algebraic system reduction" begin
    @variables a1(t) a2(t) a3(t) s1(t) s2(t) s3(t)
    @parameters p1 p2 p3

    algebraic_states = [a1, a2, a3]
    eqs = [
        0 ~ p1*a1 + p2*a2
        0 ~ s1*s2*a2 + cos(a3)
        0 ~ sin(a3)
    ]

    red_eqs, red_obs, red_states = mtkext.reduce_linear_algebraic(eqs, Equation[], algebraic_states; verbose=false)
    # eq1 can solve for a1 (linear in a1 with const coeff p1), eq2 can solve for a2 (linear in a2 with symbolic coeff s1*s2)
    # eq3 is nonlinear in a3 → stays as constraint
    @test length(red_states) == 1
    @test isequal(only(red_states), a3.val)
    @test length(red_obs) == 2

    @variables A(t) B(t) C(t) D(t) E(t) F(t) G(t) H(t) I(t) J(t) K(t) L(t) M(t) N(t)
    @parameters in1 in2 p1 p2
    eqs = [
        0 ~ sqrt(G^2 + I^2) - L     # eq1: linear in L (coeff -1)
        0 ~ atan(I, G) - J           # eq2: linear in J (coeff -1)
        0 ~ -N + in1*I + G*in2       # eq3: linear in N (coeff -1)
        0 ~ -M - in1*G + I*in2       # eq4: linear in M (coeff -1)
        0 ~ B - G                    # eq5: linear in B,G (coeff 1,-1)
        0 ~ A - I                    # eq6: linear in A,I (coeff 1,-1)
        0 ~ K - in2                  # eq7: linear in K (coeff 1)
        0 ~ -in1 + C                 # eq8: linear in C (coeff 1)
        0 ~ p1^2 - (D^2) - (E^2)    # eq9: nonlinear in D,E
        0 ~ p2 - D*H - F*E          # eq10: nonlinear in D,E,F,H (products of states)
        0 ~ -B + E                   # eq11: linear in B,E (coeff -1,1)
        0 ~ D - A                    # eq12: linear in D,A (coeff 1,-1)
        0 ~ F + K                    # eq13: linear in F,K (coeff 1,1)
        0 ~ C + H                    # eq14: linear in C,H (coeff 1,1)
    ]
    all_states = [A, B, C, D, E, F, G, H, I, J, K, L, M, N]
    red_eqs, red_obs, red_states = mtkext.reduce_linear_algebraic(eqs, Equation[], all_states; verbose=false)
    # 13 of 14 equations matched via SCC decomposition (each individually linear).
    # Only eq9 (p1^2 - D^2 - E^2) is fully nonlinear and can't be matched.
    # One state remains as a constraint variable.
    @test length(red_states) == 1
    @test length(red_eqs) == 1
    @test topologicical_sorted(red_eqs)
    @test topologicical_sorted(red_obs)

    @variables p(t) prosumer_p(t) prosumer_q_nw(t) prosumer_q_inj(t)
    @parameters q_nw prosumer_q_prosumer
    eqs = [
        0 ~ -p + prosumer_p
        0 ~ -q_nw + prosumer_q_nw
        0 ~ -prosumer_q_inj - prosumer_q_nw
        0 ~ prosumer_q_prosumer + prosumer_q_nw
    ]
    red_eqs, red_obs, red_states = mtkext.reduce_linear_algebraic(eqs, Equation[], [p, prosumer_p, prosumer_q_nw, prosumer_q_inj]; verbose=false)
    @test length(red_eqs) == 1
    @test length(red_obs) == 3
    @test topologicical_sorted(red_obs)

    # Minimal regression: chain A→B→C where C is already in obseqs.
    # The batch insertion in _insert_sorted! places A (position 0) before B (position k>0)
    # even though A depends on B. Fix requires sequential insertion + correct SCC order.
    @variables sA(t) sB(t) sC(t)
    @parameters sp
    existing_obs = Equation[sC ~ sp]
    eqs = [
        0 ~ -sA + sB,  # sA depends on sB
        0 ~ -sB + sC,  # sB depends on sC (already in existing_obs)
    ]
    red_eqs, red_obs, red_states = mtkext.reduce_linear_algebraic(eqs, existing_obs, [sA, sB]; verbose=false)
    @test isempty(red_states)
    @test topologicical_sorted(red_obs)

    # Longer 4-state dependency chain A→B→C→D (D=parameter).
    # Tests that correct SCC order is maintained across a longer chain.
    @variables chA(t) chB(t) chC(t) chD(t)
    @parameters chp
    eqs = [
        0 ~ -chA + chB,   # chA depends on chB
        0 ~ -chB + chC,   # chB depends on chC
        0 ~ -chC + chD,   # chC depends on chD
        0 ~ -chD + chp,   # chD = chp (pure parameter)
    ]
    red_eqs, red_obs, red_states = mtkext.reduce_linear_algebraic(eqs, Equation[], [chA, chB, chC, chD]; verbose=false)
    @test isempty(red_states)
    @test isempty(red_eqs)
    @test topologicical_sorted(red_obs)
    # chD must come before chC, chC before chB, chB before chA
    obs_lhs = [eq.lhs for eq in red_obs]
    @test findfirst(isequal(chD.val), obs_lhs) < findfirst(isequal(chC.val), obs_lhs)
    @test findfirst(isequal(chC.val), obs_lhs) < findfirst(isequal(chB.val), obs_lhs)
    @test findfirst(isequal(chB.val), obs_lhs) < findfirst(isequal(chA.val), obs_lhs)

    # Nonlinear dependency ordering regression:
    # nlx ~ atan(nla, nlb) has a *nonlinear* dep on nla and nlb.
    # nla and nlb in turn depend on nlc which is already in obseqs, so
    # _insert_sorted! places them at a high position. Without tracking
    # nonlinear deps in the SCC graph, nlx lands at position 1 (before nla/nlb).
    @variables nlx(t) nla(t) nlb(t) nlc(t)
    @parameters nlp
    eqs = [
        0 ~ -nlx + atan(nla, nlb),  # nlx nonlinearly depends on nla, nlb
        0 ~ -nla + nlc,             # nla linearly depends on nlc (already in obseqs)
        0 ~ -nlb + nlc,             # nlb linearly depends on nlc
    ]
    red_eqs, red_obs, red_states = mtkext.reduce_linear_algebraic(eqs, Equation[nlc ~ nlp], [nlx, nla, nlb]; verbose=false)
    @test isempty(red_states)
    @test topologicical_sorted(red_obs)
    obs_lhs = [eq.lhs for eq in red_obs]
    @test findfirst(isequal(nla.val), obs_lhs) < findfirst(isequal(nlx.val), obs_lhs)
    @test findfirst(isequal(nlb.val), obs_lhs) < findfirst(isequal(nlx.val), obs_lhs)
end

@testset "FF-blocking in reduce_linear_algebraic" begin
    # Setup: u_r, u_i are outputs (e.g. bus voltages), i_r, i_i are inputs (currents)
    # n is a non-output internal state, p1/p2 are pure parameters
    @variables u_r(t) u_i(t) n(t) foo(t)
    @parameters P Q i_r i_i p1 p2

    pq_eqs = [
        0 ~ P + u_r*i_r + i_i*u_i   # P-injection: depends on inputs i_r, i_i
        0 ~ Q + u_r*i_i - i_r*u_i   # Q-injection: depends on inputs i_r, i_i
    ]

    ff_set = Set([i_r, i_i])

    # A: Output states in input-dependent equations → both stay as constraints
    red_eqs, red_obs, red_states = mtkext.reduce_linear_algebraic(
        pq_eqs, Equation[], [u_r, u_i];
        outputs=[u_r, u_i], ff_inputs=ff_set, verbose=false)
    @test length(red_states) == 2
    @test isempty(red_obs)
    @test length(red_eqs) == 2
    # Equations should be the originals (no divisions introduced)
    @test red_eqs == pq_eqs

    # B: Same equations, ff_inputs=Set() → normal solving (blocking disabled)
    red_eqs, red_obs, red_states = mtkext.reduce_linear_algebraic(
        pq_eqs, Equation[], [u_r, u_i];
        outputs=[u_r, u_i], ff_inputs=Set(), verbose=false)
    @test isempty(red_states)
    @test length(red_obs) == 2

    # C: Non-output state in input-dependent equation → still solved (blocking only for outputs)
    eqs_n = [0 ~ n - i_r*u_r - i_i*u_i]  # n is not an output
    red_eqs, red_obs, red_states = mtkext.reduce_linear_algebraic(
        eqs_n, Equation[], [n];
        outputs=[], ff_inputs=ff_set, verbose=false)
    @test isempty(red_states)
    @test length(red_obs) == 1
    @test isequal(only(red_obs).lhs, n.val)

    # D: Output state in input-free equation → still solved (no input dependency)
    eqs_free = [0 ~ u_r - p1*p2]
    red_eqs, red_obs, red_states = mtkext.reduce_linear_algebraic(
        eqs_free, Equation[], [u_r];
        outputs=[u_r], ff_inputs=ff_set, verbose=false)
    @test isempty(red_states)
    @test length(red_obs) == 1
    @test isequal(only(red_obs).lhs, u_r.val)

    # E: Cramer's rule for non-output 2×2 coupled states: verify solution correctness
    # and absence of spurious i_r or i_i denominators
    @variables x(t) y(t)
    eqs_2x2 = [0 ~ P + x*i_r + i_i*y, 0 ~ Q + x*i_i - i_r*y]
    red_eqs, red_obs, red_states = mtkext.reduce_linear_algebraic(
        eqs_2x2, Equation[], [x, y]; verbose=false)
    @test isempty(red_states)
    @test length(red_obs) == 2
    x_sol = only(filter(eq -> isequal(eq.lhs, x.val), red_obs)).rhs
    y_sol = only(filter(eq -> isequal(eq.lhs, y.val), red_obs)).rhs
    # Expected: x = -(P*i_r + Q*i_i)/(i_r^2+i_i^2), y = (Q*i_r - P*i_i)/(i_r^2+i_i^2)
    # Verify numerically: P=3, Q=4, i_r=1, i_i=2 → denom=5, x=-11/5, y=-2/5
    sub = Dict(P => 3.0, Q => 4.0, i_r => 1.0, i_i => 2.0)
    @test unwrap_const(Symbolics.substitute(x_sol, sub)) ≈ -(3*1 + 4*2)/5
    @test unwrap_const(Symbolics.substitute(y_sol, sub)) ≈ (4*1 - 3*2)/5
    # At i_i=0, i_r=1: old symbolic_linear_solve divided by i_i → NaN; Cramer uses i_r^2+i_i^2=1
    sub2 = Dict(P => 3.0, Q => 4.0, i_r => 1.0, i_i => 0.0)
    @test unwrap_const(Symbolics.substitute(x_sol, sub2)) ≈ -3.0
    @test unwrap_const(Symbolics.substitute(y_sol, sub2)) ≈  4.0

    # F: Mixed: output blocked from input-dep eq, non-output still solved
    # n (non-output) and u_r (output) both appear in same equation set
    eqs_mixed = [
        0 ~ P + u_r*i_r + i_i*u_i   # u_r, u_i are outputs → blocked
        0 ~ Q + u_r*i_i - i_r*u_i   # same
        0 ~ n - u_r - u_i            # n is non-output, depends only on states (not inputs)
    ]
    red_eqs, red_obs, red_states = mtkext.reduce_linear_algebraic(
        eqs_mixed, Equation[], [u_r, u_i, n];
        outputs=[u_r, u_i], ff_inputs=ff_set, verbose=true)
    @test length(red_states) == 2   # u_r, u_i stay as constraints
    @test length(red_obs) == 1      # n solved from its equation
    @test isequal(only(red_obs).lhs, n.val)

    # G: PV-like case: output blocked from linear P-eq, nonlinear V-eq stays
    @parameters Vset Pset
    eqs_pv = [
        0 ~ Pset + u_r*i_r + i_i*u_i   # linear in u_r, depends on inputs → blocked
        0 ~ -(Vset^2) + u_r^2 + u_i^2  # nonlinear → can't be solved for either state
    ]
    red_eqs, red_obs, red_states = mtkext.reduce_linear_algebraic(
        eqs_pv, Equation[], [u_r, u_i];
        outputs=[u_r, u_i], ff_inputs=ff_set, verbose=false)
    @test length(red_states) == 2   # both stay
    @test isempty(red_obs)
    @test length(red_eqs) == 2      # both original equations preserved intact

    # H: 4-SCC FF chain: input→LC→LS→LS→LC→output
    # cn1, cn2 are dynamic states whose values appear as coefficients,
    # making SCCs 2 and 3 :linear_state respectively.
    # Walking from the output SCC (LC) backward, SCC3 is the first :linear_state hit
    # → SCC3 must be forbidden; SCCs 1, 2 and 4 must still be solved.
    @variables c1(t) c2(t) c3(t) c4(t) cn1(t) cn2(t)
    @parameters c_inp
    ff_chain_eqs = [
        Dt(cn1) ~ 0,          # cn1 is a dynamic state (coefficient role → SCC2 is :LS)
        Dt(cn2) ~ 0,          # cn2 is a dynamic state (coefficient role → SCC3 is :LS)
        0 ~ c1 - c_inp,       # SCC1: c1 = c_inp,  coeff of c1 = 1    (LC), has_ff
        0 ~ cn1*c2 - c1,      # SCC2: cn1*c2 = c1, coeff of c2 = cn1  (LS)
        0 ~ cn2*c3 - c2,      # SCC3: cn2*c3 = c2, coeff of c3 = cn2  (LS)
        0 ~ c4 - c3,          # SCC4: c4 = c3,     coeff of c4 = 1    (LC), output
    ]
    all_chain_states = [cn1, cn2, c1, c2, c3, c4]
    red_eqs, red_obs, red_states = mtkext.reduce_linear_algebraic(
        ff_chain_eqs, Equation[], all_chain_states;
        outputs=[c4], ff_inputs=Set([c_inp]), verbose=false)
    obs_lhs = Set(eq.lhs for eq in red_obs)
    # SCCs 1 and 2 are not forbidden
    @test c1.val ∈ obs_lhs
    @test c2.val ∈ obs_lhs
    # SCC3 is forbidden (first :LS from output end): c3 stays as a constraint state
    @test c3.val ∉ obs_lhs
    @test c3.val ∈ Set(red_states)
    # SCC4 (LC) is not forbidden: c4 solved as obs (depends on state c3, no direct FF)
    @test c4.val ∈ obs_lhs
    @test topologicical_sorted(red_obs)
end

@testset "CPL bus: matching prefers linear_const over linear_state" begin
    # Reproduces the constant-power-load bus structure with ₊-separated names.
    # 16 implicit algebraic equations, 16 states. Outputs: busbar₊u_r, busbar₊u_i.
    # FF inputs: busbar₊i_r, busbar₊i_i.
    #
    # The bug: the column sort key `-count('₊', ...)` causes deeply-nested states
    # like load₊terminal₊u_r (2 ₊) to be tried before load₊P (1 ₊) in the bipartite
    # matching. The P-balance equation gets matched to load₊terminal₊u_r (linear_state)
    # instead of load₊P (linear_const). This pulls output-aliased states into an SCC
    # that is both ff_source and output_sink, causing the group to be skipped.
    #
    # After reduction only 2 states should remain (the output voltages as algebraic
    # constraints defined by the CPL equations).
    @variables begin
        busbar₊u_r(t); busbar₊u_i(t); busbar₊P(t); busbar₊Q(t); busbar₊u_mag(t); busbar₊u_arg(t)
        busbar₊terminal₊u_r(t); busbar₊terminal₊u_i(t); busbar₊terminal₊i_r(t); busbar₊terminal₊i_i(t)
        load₊P(t); load₊Q(t)
        load₊terminal₊u_r(t); load₊terminal₊u_i(t); load₊terminal₊i_r(t); load₊terminal₊i_i(t)
    end
    @parameters busbar₊i_r busbar₊i_i load₊Pset load₊Qset

    eqs = [
        0 ~ -busbar₊P - busbar₊i_i*busbar₊u_i - busbar₊u_r*busbar₊i_r
        0 ~ -busbar₊Q + busbar₊i_i*busbar₊u_r - busbar₊u_i*busbar₊i_r
        0 ~ sqrt(busbar₊u_r^2 + busbar₊u_i^2) - busbar₊u_mag
        0 ~ atan(busbar₊u_i, busbar₊u_r) - busbar₊u_arg
        0 ~ busbar₊terminal₊u_r - busbar₊u_r
        0 ~ busbar₊terminal₊u_i - busbar₊u_i
        0 ~ busbar₊terminal₊i_r - busbar₊i_r
        0 ~ -busbar₊i_i + busbar₊terminal₊i_i
        0 ~ -load₊P + load₊terminal₊u_i*load₊terminal₊i_i + load₊terminal₊i_r*load₊terminal₊u_r
        0 ~ -load₊Q + load₊terminal₊u_i*load₊terminal₊i_r - load₊terminal₊i_i*load₊terminal₊u_r
        0 ~ -load₊terminal₊i_r + (load₊Pset*load₊terminal₊u_r + load₊Qset*load₊terminal₊u_i) / (load₊terminal₊u_i^2 + load₊terminal₊u_r^2)
        0 ~ (load₊Pset*load₊terminal₊u_i - load₊Qset*load₊terminal₊u_r) / (load₊terminal₊u_i^2 + load₊terminal₊u_r^2) - load₊terminal₊i_i
        0 ~ -busbar₊terminal₊u_r + load₊terminal₊u_r
        0 ~ load₊terminal₊u_i - busbar₊terminal₊u_i
        0 ~ load₊terminal₊i_r + busbar₊terminal₊i_r
        0 ~ busbar₊terminal₊i_i + load₊terminal₊i_i
    ]
    all_states = [busbar₊u_r, busbar₊u_i, busbar₊P, busbar₊Q, busbar₊u_mag, busbar₊u_arg,
                  busbar₊terminal₊u_r, busbar₊terminal₊u_i, busbar₊terminal₊i_r, busbar₊terminal₊i_i,
                  load₊P, load₊Q, load₊terminal₊u_r, load₊terminal₊u_i, load₊terminal₊i_r, load₊terminal₊i_i]
    outputs = [busbar₊u_r, busbar₊u_i]
    ff_set = Set([busbar₊i_r, busbar₊i_i])

    red_eqs, red_obs, red_states = mtkext.reduce_linear_algebraic(
        eqs, Equation[], all_states;
        outputs=outputs, ff_inputs=ff_set, verbose=false)

    # After reduction only 2 states should remain as algebraic constraints.
    # The remaining states may be the output variables or their aliases (before
    # pick_best_alias_names resolves them), but there must be exactly 2.
    @test length(red_states) == 2
    @test length(red_eqs) == 2
    @test length(red_obs) == 14
    @test topologicical_sorted(red_obs)
end

@testset "Test get_alias function" begin
    @variables a(t) b(t) c(t)
    @parameters p

    @test mtkext.get_alias(a ~ b) == (a, b)
    @test mtkext.get_alias(0 ~ b - b) == nothing
    @test mtkext.get_alias(p ~ 1) == nothing
    # parameters must not be treated as aliases — equation like `x ~ V` defines x,
    # it doesn't alias two unknowns (regression: VoltageSource p₊v ~ V was falsely aliased)
    @test mtkext.get_alias(p ~ a) == nothing
    @test mtkext.get_alias(a ~ p) == nothing
    @test mtkext.get_alias(0 ~ p - a) == nothing
    @test mtkext.get_alias(p ~ -a) == nothing
    @test mtkext.get_alias(0 ~ p-sin(a)) == nothing
end

@testset "Output defined by parameter (NWTerminal regression)" begin
    @mtkmodel NWTerminal begin
        @variables begin
            v(t), [description="Voltage at node"]
            i(t), [description="Current flowing into node"]
        end
    end
    @mtkmodel VoltageSource begin
        @components begin
            p = NWTerminal()
        end
        @parameters begin
            V = 1.0
        end
        @variables begin
            i(t), [description="produced current"]
        end
        @equations begin
            i ~ -p.i
            p.v ~ V
        end
    end
    @named vs = VoltageSource()
    # p₊v ~ V (output defined by parameter) must not be treated as an alias
    vs_vertex = VertexModel(vs, [:p₊i], [:p₊v])
    @test Set(NetworkDynamics.outsym(vs_vertex)) == Set([:p₊v])
    @test :V ∈ Set(NetworkDynamics.psym(vs_vertex))
end

@__MODULE__()==Main ? includet(joinpath(pkgdir(NetworkDynamics), "test", "ComponentLibrary.jl")) : (const Lib = Main.Lib)
@testset "Test Component Lib" begin
    m = Lib.swing_mtk()
    @test Set(sym(m)) == Set([:θ, :ω])
    m = Lib.line_mtk()
    @test isempty(sym(m))
    m = Lib.dqbus_swing()
    @test Set(sym(m)) == Set([:θ, :ω])
    m = Lib.dqbus_pq()
    @test Set(sym(m)) == Set([:u_r, :u_i])
    m = Lib.dqbus_timedeppq(Pfun = _->1.0)
    @test Set(sym(m)) == Set([:u_r, :u_i])
    m = Lib.dqbus_pv()
    @test Set(sym(m)) == Set([:u_r, :u_i])
    m = Lib.dqbus_pv(injector=true)
    @test Set(sym(m)) == Set([:i_r, :i_i])
    m = Lib.dqbus_slack()
    @test isempty(sym(m))
    m = Lib.dqbus_slack(injector=true)
    @test Set(sym(m)) == Set([:i_r, :i_i])
    m = Lib.dqline(R=0.1, X=0.1)
    @test isempty(sym(m))
    m = Lib.dqbus_swing_and_load()
    @test Set(sym(m)) == Set([:swing₊θ, :swing₊ω])
    m = Lib.dqbus_swing_injector()
    @test Set(sym(m)) == Set([:θ, :ω, :i_r, :i_i])
    m = Lib.dqbus_pq_injector()
    @test isempty(sym(m))
    m = Lib.dqbus_shunt_hub()
    @test Set(sym(m)) == Set([:u_r, :u_i])

    Lib.powergridlike_network()
    Lib.powergridlike_injector_network()
end

@testset "Test warn on missing features" begin
    @component function testcomp()
        @variables x(t) y(t)
        @discretes s(t)
        System([x ~ y + s, s ~ 1], t, [x, y, s], []; name=:test)
    end
    sys = testcomp()
    @test_logs (:warn, r"Model contains .* @discretes") VertexModel(testcomp(), [:x], [:y])
end

@testset "promotion of bindings to init_formulas" begin
    @testset "no nesting" begin
        @component function slack_differential(; name)
            @parameters u_init_r=1 u_init_i=0
            @named busbar = BusBase(; u_r = u_init_r + 0.1, u_i = u_init_i)
            eqs = [Dt(busbar.u_r) ~ 0, Dt(busbar.u_i) ~ 0]
            System(eqs, t, [], [u_init_r, u_init_i]; name, systems=[busbar])
        end
        @named sys = slack_differential()
        vm = VertexModel(sys, [:busbar₊i_r, :busbar₊i_i], [:busbar₊u_r, :busbar₊u_i]; verbose=false)

        # formulas are automatically attached to the vertex model
        @test has_initformula(vm)
        formulas = collect(get_initformulas(vm))
        @test length(formulas) == 2

        f_ur = only(filter(f -> f.outsym == [:busbar₊u_r], formulas))
        f_ui = only(filter(f -> f.outsym == [:busbar₊u_i], formulas))

        out = NetworkDynamics.SymbolicView(zeros(1), f_ur.outsym)
        f_ur(out, NetworkDynamics.SymbolicView([2.0], f_ur.sym))
        @test out[:busbar₊u_r] ≈ 2.1

        out = NetworkDynamics.SymbolicView(zeros(1), f_ur.outsym)
        f_ur(out, NetworkDynamics.SymbolicView([0.0], f_ur.sym))
        @test out[:busbar₊u_r] ≈ 0.1

        out = NetworkDynamics.SymbolicView(zeros(1), f_ui.outsym)
        f_ui(out, NetworkDynamics.SymbolicView([3.5], f_ui.sym))
        @test out[:busbar₊u_i] ≈ 3.5
    end

    @testset "nested model" begin
        @component function slack_diff_inner(; name)
            @parameters u_init_r=1
            @named busbar = BusBase(; u_r = u_init_r)
            eqs = [Dt(busbar.u_r) ~ 0, Dt(busbar.u_i) ~ 0]
            System(eqs, t, [], [u_init_r]; name, systems=[busbar])
        end
        @component function slack_diff_outer(; name)
            @variables u_r(t) u_i(t) i_r(t) i_i(t)
            @named inner = slack_diff_inner()
            eqs = [
                u_r ~ inner.busbar.u_r
                u_i ~ inner.busbar.u_i
                i_r ~ inner.busbar.i_r
                i_i ~ inner.busbar.i_i
            ]
            System(eqs, t; name, systems=[inner])
        end
        # Outer model has no bindings of its own; inner model has u_r = u_init_r binding.
        # bindings_to_initformulas should propagate the inner binding with the namespaced symbol.
        @named outer = slack_diff_outer()
        vm = VertexModel(outer, [:i_r, :i_i], [:u_r, :u_i]; verbose=false)

        @test has_initformula(vm)
        formulas = collect(get_initformulas(vm))
        @test length(formulas) == 1
        f = only(formulas)
        # inner₊busbar₊u_r is aliased to u_r after pick_best_alias_names
        @test f.outsym == [:u_r]
        @test f.sym    == [:inner₊u_init_r]

        out = NetworkDynamics.SymbolicView(zeros(1), f.outsym)
        f(out, NetworkDynamics.SymbolicView([4.2], f.sym))
        @test out[:u_r] ≈ 4.2

        @testset "warn on duplicate formula targets" begin
            @component function slack_diff_outer_dup(; name)
                @parameters u_init_outer=2.0
                @variables u_r(t) u_i(t) i_r(t) i_i(t)
                @named inner = slack_diff_inner()
                eqs = [
                    u_r ~ inner.busbar.u_r
                    u_i ~ inner.busbar.u_i
                    i_r ~ inner.busbar.i_r
                    i_i ~ inner.busbar.i_i
                ]
                # outer also binds u_r: after obs_subs both this and inner₊busbar₊u_r alias to u_r
                System(eqs, t; name, systems=[inner], bindings=[u_r => u_init_outer])
            end
            @named outer_dup = slack_diff_outer_dup()
            @test_logs (:warn, r"Multiple InitFormula") match_mode=:any begin
                VertexModel(outer_dup, [:i_r, :i_i], [:u_r, :u_i]; verbose=false)
            end
        end
    end
end

@testset "promotion of guesses to guess_formulas" begin
    @component function vertex_with_symbolic_guess(; name)
        @parameters u_init=1.0
        @variables x(t) y(t)
        eqs = [Dt(x) ~ -x; Dt(y) ~ -y]
        System(eqs, t, [x, y], [u_init]; name, guesses=[y => 2*u_init])
    end
    @named sys = vertex_with_symbolic_guess()
    vm = VertexModel(sys, [], [:x, :y]; verbose=false)

    @test has_guessformula(vm)
    formulas = collect(get_guessformulas(vm))
    @test length(formulas) == 1
    f = only(formulas)
    @test f.outsym == [:y]
    @test f.sym    == [:u_init]

    out = NetworkDynamics.SymbolicView(zeros(1), f.outsym)
    f(out, NetworkDynamics.SymbolicView([3.0], f.sym))
    @test out[:y] ≈ 6.0

    out = NetworkDynamics.SymbolicView(zeros(1), f.outsym)
    f(out, NetworkDynamics.SymbolicView([0.5], f.sym))
    @test out[:y] ≈ 1.0

    # constant guess is NOT promoted to a formula
    @test all(gf -> :x ∉ gf.outsym, get_guessformulas(vm))
end

@testset "bound parameters survive mtkcompile" begin
    # Bound parameters (params whose default references another param, like Sn = S_b)
    # are removed from parameters(sys) by mtkcompile/complete via
    # remove_bound_parameters_from_ps, but still appear symbolically in equations.
    # ND must include them in allparams so build_function sees them and
    # InitFormula targeting them can be validated.
    @mtkmodel BoundParamNode begin
        @variables begin
            u(t), [description = "input"]
            x(t) = 0.0, [description = "state"]
        end
        @parameters begin
            S_b = 100.0, [description = "system base power"]
            Sn  = S_b,   [description = "machine nominal power (bound to S_b)"]
        end
        @equations begin
            Dt(x) ~ u / Sn - x   # Sn must survive as a real parameter
        end
    end
    @named node = BoundParamNode()
    vm = VertexModel(node, [:u], [:x])
    # Sn is independently settable (not just an alias for S_b)
    @test :Sn  ∈ Set(NetworkDynamics.psym(vm))
    @test :S_b ∈ Set(NetworkDynamics.psym(vm))
    # The InitFormula Sn = S_b must be attached
    @test has_initformula(vm)
    f = only(filter(f -> f.outsym == [:Sn], collect(get_initformulas(vm))))
    out = NetworkDynamics.SymbolicView(zeros(1), f.outsym)
    f(out, NetworkDynamics.SymbolicView([42.0], f.sym))
    @test out[:Sn] ≈ 42.0   # Sn gets value from S_b
end
