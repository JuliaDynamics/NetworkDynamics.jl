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

@testset "Woraround tests" begin
    @variables x(t)
    # if this test gets fixed, we can get rid of workaround get_variables_fix
    @test get_variables(Dt(x)) == Set(Dt(x))
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
    mtkext.eq_type(eq) == (:explicit_algebraic, y.val)

    eq = y ~ x + b + y
    mtkext.eq_type(eq) == (:implicit_algebraic, y.val)

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
    vf = VertexModel(mtkbus, [:busbar₊i_r, :busbar₊i_i], [:busbar₊u_r, :busbar₊u_i], verbose=false)

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

    # E: Input dependency via obs chain → output should be blocked
    # foo is observed as foo ~ i_r + 1 (depends on input i_r)
    # equation 0 ~ u_r + foo looks input-free without expanding obs
    existing_obs = Equation[foo ~ i_r + 1]
    eqs_chain = [0 ~ u_r + foo]
    red_eqs, red_obs, red_states = mtkext.reduce_linear_algebraic(
        eqs_chain, existing_obs, [u_r];
        outputs=[u_r], ff_inputs=ff_set, verbose=false)
    @test length(red_states) == 1  # u_r stays as state (obs chain detected)
    @test length(red_eqs) == 1

    # F: Cramer's rule for non-output 2×2 coupled states: verify solution correctness
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

    # G: Mixed: output blocked from input-dep eq, non-output still solved
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

    # H: PV-like case: output blocked from linear P-eq, nonlinear V-eq stays
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

####
#### remove_aliases tests
####
@testset "remove_aliases" begin
    using Symbolics: unwrap
    @variables a(t) b(t) c(t) d(t) x(t) y(t) z(t)
    @parameters p q

    # shorthand: call remove_aliases on copies, check length invariant, return results
    function RA(eqs, obseqs, states, outputs; kw...)
        eqs_c, obs_c, st_c = copy(eqs), copy(obseqs), copy(states)
        ret = mtkext.remove_aliases(eqs_c, obs_c, st_c, outputs; verbose=false, kw...)
        @assert length(ret[1]) == length(ret[3]) "length mismatch: $(length(ret[3])) states vs $(length(ret[1])) eqs"
        ret
    end
    # Order-independent alias check: true if obs contains an alias equation between v1 and v2
    function is_alias_between(obs, v1, v2)
        any(o -> (isequal(o.lhs, unwrap(v1)) && isequal(o.rhs, unwrap(v2))) ||
                 (isequal(o.lhs, unwrap(v2)) && isequal(o.rhs, unwrap(v1))), obs)
    end

    @testset "No aliases — system unchanged" begin
        eqs = [Dt(x) ~ -x, 0 ~ y + x]
        states = unwrap.([x, y])
        ret_eqs, ret_obs, ret_st = RA(eqs, Equation[], states, [])

        @test isequal(ret_eqs, eqs)
        @test isequal(ret_st, states)
        @test isempty(ret_obs)
    end

    @testset "Explicit alias (a ~ b), no outputs" begin
        eqs    = [a ~ b]
        states = [unwrap(b)]
        ret_eqs, ret_obs, ret_st = RA(eqs, Equation[], states, [])
        @test isempty(ret_eqs)
        @test isempty(ret_st)
        @test length(ret_obs) == 1
        @test is_alias_between(ret_obs, a, b)
    end

    @testset "Explicit alias (a ~ b), b is output → b is main, obs: a ~ b" begin
        eqs    = [a ~ b]
        states = unwrap.([b])
        ret_eqs, ret_obs, ret_st = RA(eqs, Equation[], states, [unwrap(b)])
        @test isempty(ret_eqs)
        @test isempty(ret_st)
        # b is output → b is main → a becomes alias of b
        @test ret_obs == [a ~ b]
    end

    @testset "Explicit alias (a ~ b), a is output → a is main, obs: b ~ a" begin
        eqs    = [a ~ b]
        states = [unwrap(b)]
        ret_eqs, ret_obs, ret_st = RA(eqs, Equation[], states, [unwrap(a)])
        @test isempty(ret_eqs)
        @test isempty(ret_st)
        @test length(ret_obs) == 1
        # a is output → a is main → b becomes alias of a
        @test isequal(only(ret_obs).lhs, unwrap(b))
        @test isequal(only(ret_obs).rhs, unwrap(a))
    end

    @testset "Namespace hierarchy: fewer ₊ wins as main" begin
        # Simulate connector aliasing: sub₊a ~ b (sub₊a has 1 '₊', b has 0)
        # → b is main, sub₊a is rest → obs: sub₊a ~ b
        @variables sub₊a(t)
        eqs    = [sub₊a ~ b]
        states = [unwrap(b)]
        ret_eqs, ret_obs, ret_st = RA(eqs, Equation[], states, [])
        @test isempty(ret_eqs)
        @test isempty(ret_st)
        @test length(ret_obs) == 1
        @test isequal(only(ret_obs).lhs, unwrap(sub₊a))
        @test isequal(only(ret_obs).rhs, unwrap(b))
    end

    @testset "Explicit alias coexists with diff eq — diff eq unchanged" begin
        # Dt(x) ~ -x  and  a ~ b
        eqs    = [Dt(x) ~ -x, a ~ b]
        states = [unwrap(x), unwrap(b)]
        ret_eqs, ret_obs, ret_st = RA(eqs, Equation[], states, [])
        @test length(ret_eqs) == 1
        @test length(ret_st)  == 1
        @test isequal(only(ret_eqs), Dt(x) ~ -x)
        @test isequal(only(ret_st), unwrap(x))
        @test length(ret_obs) == 1
    end

    @testset "Multiple explicit aliases, no implicit_algebraic eqs (the bug case)" begin
        # 3 alias eqs + 1 diff eq: before fix the @assert fired because no implicit_algebraic
        # eqs existed to pair with the alias state removal search
        eqs    = [a ~ b, c ~ d, Dt(x) ~ -x, y ~ z]
        states = [unwrap(b), unwrap(d), unwrap(x), unwrap(z)]
        ret_eqs, ret_obs, ret_st = RA(eqs, Equation[], states, [])
        @test isequal(ret_eqs, [Dt(x) ~ -x])
        @test isequal(ret_st,  [unwrap(x)])
        @test length(ret_obs) == 3
        @test is_alias_between(ret_obs, a, b)
        @test is_alias_between(ret_obs, c, d)
        @test is_alias_between(ret_obs, y, z)
    end

    @testset "Multiple explicit aliases, one output pinned as main" begin
        # a ~ b, c ~ d: d is output → d is deterministically the main → obs: c ~ d
        eqs    = [a ~ b, c ~ d]
        states = [unwrap(b), unwrap(d)]
        ret_eqs, ret_obs, ret_st = RA(eqs, Equation[], states, [unwrap(d)])
        @test isempty(ret_eqs)
        @test isempty(ret_st)
        @test length(ret_obs) == 2
        @test any(o -> isequal(o, c ~ d), ret_obs)   # d pinned as output → d is main
        @test is_alias_between(ret_obs, a, b)
    end

    @testset "Implicit alias (0 ~ a - b), state from pool" begin
        # eq_type(0 ~ a - b) = (:implicit_algebraic, nothing)
        # simulate: pool gives `a` as the state for this equation
        eqs    = [0 ~ a - b]
        states = [unwrap(a)]   # a popped from implicit_states pool
        ret_eqs1, ret_obs1, ret_st1 = RA(eqs, Equation[], states, [])
        @test isempty(ret_eqs1)
        @test isempty(ret_st1)
        @test length(ret_obs1) == 1
        @test is_alias_between(ret_obs1, a, b)

        # same canonical alias regardless of which variable was the pool state
        eqs    = [0 ~ a - b]
        states = [unwrap(b)]
        ret_eqs2, ret_obs2, ret_st2 = RA(eqs, Equation[], states, [])
        @test isempty(ret_eqs2)
        @test isempty(ret_st2)
        @test isequal(ret_obs1, ret_obs2)
    end

    @testset "Implicit alias + explicit alias together" begin
        # 0 ~ a - b (implicit, pool gives a), c ~ d (explicit, state=d), Dt(x) ~ -x
        eqs    = [0 ~ a - b, c ~ d, Dt(x) ~ -x]
        states = [unwrap(a), unwrap(d), unwrap(x)]
        ret_eqs, ret_obs, ret_st = RA(eqs, Equation[], states, [])
        @test isequal(ret_eqs, [Dt(x) ~ -x])
        @test isequal(ret_st,  [unwrap(x)])
        @test length(ret_obs) == 2
        @test is_alias_between(ret_obs, a, b)
        @test is_alias_between(ret_obs, c, d)
    end

    @testset "Transitive chain a ~ b, b ~ c → single alias group" begin
        # a ~ b  (state=b),  b ~ c  (state=c): all three must end up in one alias group
        eqs    = [a ~ b, b ~ c]
        states = [unwrap(b), unwrap(c)]
        ret_eqs, ret_obs, ret_st = RA(eqs, Equation[], states, [])
        @test isempty(ret_eqs)
        @test isempty(ret_st)
        # two non-mains aliased to a single main (all three vars covered)
        @test length(ret_obs) == 2
        @test isequal(ret_obs[1].rhs, ret_obs[2].rhs)  # same main
        covered = [ret_obs[1].lhs, ret_obs[2].lhs, ret_obs[1].rhs]
        @test all(v -> any(isequal(v, c2) for c2 in covered), unwrap.([a, b, c]))
    end

    @testset "Alias in obseqs picked up for aliasgroup" begin
        # no aliases in main eqs; obseqs contains a ~ b — should be picked up and canonicalized
        eqs    = [Dt(x) ~ -x]
        states = [unwrap(x)]
        obseqs = [a ~ b]
        ret_eqs, ret_obs, ret_st = RA(eqs, obseqs, states, [])
        @test isequal(ret_eqs, [Dt(x) ~ -x])
        @test isequal(ret_st,  [unwrap(x)])
        @test length(ret_obs) == 1
        @test is_alias_between(ret_obs, a, b)
    end

    @testset "Alias variable substituted into remaining equations" begin
        # a is output → a is deterministically the main; b→a substitution in Dt(x) ~ b
        eqs    = [a ~ b, Dt(x) ~ b]
        states = [unwrap(b), unwrap(x)]
        ret_eqs, ret_obs, ret_st = RA(eqs, Equation[], states, [unwrap(a)])
        @test isequal(ret_eqs, [Dt(x) ~ a])
        @test isequal(ret_st,  [unwrap(x)])
        @test isequal(ret_obs, [b ~ a])
    end

    @testset "Alias substituted into obseqs" begin
        # a is output → a is deterministically the main (b→a); obs y ~ a*x unchanged (no b)
        # alias obs b ~ a inserted before y ~ a*x in topological order
        eqs    = [a ~ b, Dt(x) ~ -x]
        states = [unwrap(b), unwrap(x)]
        obseqs = [y ~ a*x]
        ret_eqs, ret_obs, ret_st = RA(eqs, obseqs, states, [unwrap(a)])
        @test isequal(ret_eqs, [Dt(x) ~ -x])
        @test isequal(ret_st,  [unwrap(x)])
        @test isequal(ret_obs, [b ~ a, y ~ a*x])
    end

    @testset "SwingBus-like: multiple connector aliases, diff states" begin
        # Mirrors the failing precomp_workload case:
        # explicit alias eqs for connector variables alongside diff eqs for θ and ω.
        # ur, ui are outputs → mains; sub₊ur, sub₊ui become obs
        @variables θ(t) ω(t) Pel(t) ur(t) ui(t) sub₊ur(t) sub₊ui(t)
        diff_eqs   = [Dt(θ) ~ ω, Dt(ω) ~ -(0.1*ω + Pel)]
        alias_eqs  = [sub₊ur ~ ur, sub₊ui ~ ui]
        eqs    = vcat(diff_eqs, alias_eqs)
        states = [unwrap(θ), unwrap(ω), unwrap(ur), unwrap(ui)]
        ret_eqs, ret_obs, ret_st = RA(eqs, Equation[], states, [unwrap(ur), unwrap(ui)])
        @test isequal(ret_eqs, diff_eqs)
        @test isequal(ret_st,  [unwrap(θ), unwrap(ω)])
        @test length(ret_obs) == 2
        @test any(o -> isequal(o, sub₊ur ~ ur), ret_obs)
        @test any(o -> isequal(o, sub₊ui ~ ui), ret_obs)
    end
end;

