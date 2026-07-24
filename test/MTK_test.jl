using ModelingToolkitBase
using ModelingToolkitBase: t_nounits as t, D_nounits as Dt
using NetworkDynamics
using NetworkDynamics: AliasMap, get_aliasmap
using OrdinaryDiffEqTsit5
using LinearAlgebra
using Graphs
using Chairmarks: @b
using Test
using SciCompDSL
mtkext = Base.get_extension(NetworkDynamics, :NetworkDynamicsMTKExt)

@testset "get_variables_deriv test" begin
    @variables x(t) y(t)
    # get_variables treats D(x) as atomic (intended behavior); get_variables_deriv unwraps it to x
    @test get_variables(Dt(x) + y) == Set([Dt(x), y])
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
    @test mtkext.eq_type(eq) == (:implicit_algebraic, nothing)

    eq = 0 ~ x+y+b
    @test mtkext.eq_type(eq) == (:implicit_algebraic, nothing)

    eq = y^2 ~ x
    @test mtkext.eq_type(eq) == (:implicit_algebraic, nothing)

    # non zero on the lhs is not expected
    eq = 1 ~ x+y+b
    @test mtkext.eq_type(eq) == (:implicit_algebraic, nothing)

    eq = Dt(x) ~ Dt(x)
    @test mtkext.eq_type(eq) == (:implicit_diffeq, nothing)

    eq = 0 ~ Dt(x)
    @test mtkext.eq_type(eq) == (:implicit_diffeq, nothing)

    eq = a ~ x+y
    @test mtkext.eq_type(eq) == (:explicit_algebraic, a.val)
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
        D = 0.1, [guess=0.1, description = "Damping"]
        Pmech, [description = "Mechanical Power"]
    end
    @equations begin
        Dt(θ) ~ ω
        Dt(ω) ~ 1/M * (Pmech - D*ω + P)
    end
end;

@named swing = SwingNode()
v = VertexModel(swing, [:P], [:θ])
@test v.mass_matrix == Diagonal([1,1])

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
@test dim(e) == 0

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
    # @test_throws r"outputs .* do not appear in the equations" VertexModel(fullyimplicit, [:u], [:z])
    VertexModel(fullyimplicit, [:u], [:z])
    # works when we assume io coupling
    VertexModel(fullyimplicit, [:u], [:z]; assume_io_coupling=true)
    # from init_tutorial
    # dependent MTKModels need to be defined at top level, so they are in front of the testset
    @named prosumer = StaticProsumerNode() # consumer
    VertexModel(prosumer, [:q̃_nw], [:p])
    @named prosumer_wrapped = Wrapper()
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

    v = VertexModel(inverter, [:i_r, :i_i], [:u_r, :u_i]; verbose=false, assume_io_coupling=false)
    @test sum(v.mass_matrix) == dim(v) - 2

    # With assume_io_coupling=true, the construction should succeed
    # This forces MTK to recognize the input->output coupling even with derivatives
    v = VertexModel(inverter, [:i_r, :i_i], [:u_r, :u_i];
                    verbose=false, assume_io_coupling=true)
    @test sum(v.mass_matrix) == dim(v) - 2
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

ST = mtkext.ST
function _reduce_equations(eqs, obs, states; outset=[], ff_inputs=[], verbose=false)
    mtkext.reduce_equations(
        Vector{Equation}(eqs),
        Vector{Equation}(obs),
        Vector{ST}(states);
        outset=Set{ST}(outset),
        ff_inputs=Set{ST}(ff_inputs),
        verbose
    )
end

@testset "Algebraic system reduction" begin
    using Symbolics: unwrap
    @variables a1(t) a2(t) a3(t) s1(t) s2(t) s3(t)
    @parameters p1 p2 p3

    algebraic_states = [a1, a2, a3]
    eqs = [
        0 ~ p1*a1 + p2*a2
        0 ~ s1*s2*a2 + cos(a3)
        0 ~ sin(a3)
    ]

    red_eqs, red_obs, red_states = _reduce_equations(eqs, Equation[], algebraic_states; verbose=false)
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
    red_eqs, red_obs, red_states = _reduce_equations(eqs, Equation[], all_states; verbose=false)
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
    red_eqs, red_obs, red_states = _reduce_equations(eqs, Equation[], [p, prosumer_p, prosumer_q_nw, prosumer_q_inj]; verbose=false)
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
    red_eqs, red_obs, red_states = _reduce_equations(eqs, existing_obs, [sA, sB]; verbose=false)
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
    red_eqs, red_obs, red_states = _reduce_equations(eqs, Equation[], [chA, chB, chC, chD]; verbose=false)
    @test isempty(red_states)
    @test isempty(red_eqs)
    @test topologicical_sorted(red_obs)

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
    red_eqs, red_obs, red_states = _reduce_equations(eqs, Equation[nlc ~ nlp], [nlx, nla, nlb]; verbose=false)
    @test isempty(red_states)
    @test topologicical_sorted(red_obs)
    obs_lhs = [eq.lhs for eq in red_obs]
    @test findfirst(isequal(nla.val), obs_lhs) < findfirst(isequal(nlx.val), obs_lhs)
    @test findfirst(isequal(nlb.val), obs_lhs) < findfirst(isequal(nlx.val), obs_lhs)

    # Handle solved extra states:
    # x_hs and y_hs are both differential states. The algebraic constraint eq3 contains
    # only x_hs and y_hs — both "extended states" (inner vars of D(.)) — so
    # needs_extension=true in _build_coeff_mat.
    # The matching eliminates one of {x_hs, y_hs} (the choice is symmetric; which one
    # is picked depends on internal tie-breaking). The "handle solved extra states" path
    # detects the conflict between the diff obs and the alg obs for the eliminated state,
    # differentiates the alg obs, and substitutes known differentials to derive 0 ~ ±2*z_hs
    # → z_hs = 0.  Final: one of {x_hs, y_hs} is the sole diff state, the other mirrors
    # it, and z_hs ~ 0.
    @variables x_hs(t) y_hs(t) z_hs(t)
    eqs = [
        Dt(x_hs) ~ z_hs,    # diff eq for x_hs, D(x_hs) → x_hs extended
        Dt(y_hs) ~ -z_hs,   # diff eq for y_hs, D(y_hs) → y_hs extended
        0 ~ x_hs - y_hs,    # algebraic: allsyms={x_hs,y_hs} ⊆ extended_states_set → needs_extension
    ]
    red_eqs, red_obs, red_states = _reduce_equations(eqs, Equation[], [x_hs, y_hs, z_hs], verbose=true)
    @test length(red_states) == 1
    # x_hs and y_hs are symmetric; the matching may eliminate either one
    remaining_hs = only(red_states)
    other_hs = isequal(remaining_hs, x_hs.val) ? y_hs.val : x_hs.val
    @test isequal(remaining_hs, x_hs.val) || isequal(remaining_hs, y_hs.val)
    @test length(red_eqs) == 1 && mtkext.isdifferential(only(red_eqs).lhs)
    @test length(red_obs) == 2
    obs_dict_hs = Dict(eq.lhs => eq.rhs for eq in red_obs)
    @test isequal(obs_dict_hs[other_hs], remaining_hs)     # eliminated state mirrors remaining
    @test isequal(obs_dict_hs[z_hs.val], Num(0).val)       # z_hs forced to 0 by consistency

    # Substitute known differentials:
    # sx and sy_sd are two states; sy_sd's equation contains D(sx) with coefficient 2
    # (:linear_const, match cost 1), so eq1 (explicit, cost 0) wins the match for D(sx).
    # sy_sd cannot be solved in round 1 because the presence of D(sx) makes eq2
    # fall into the "has differential" branch (all other states marked :unsolvable).
    # After D(sx) is solved from eq1, the "substitute known differentials" path
    # replaces D(sx) with -sx in the remaining eq2: 0 ~ 2*(-sx) + sy_sd.
    # sy_sd is then solved in round 2 as sy_sd ~ 2*sx.
    @variables sx(t) sy_sd(t)
    eqs = [
        Dt(sx) ~ -sx,           # D(sx) solved first (explicit coeff → cost 0)
        0 ~ 2*Dt(sx) + sy_sd,   # D(sx) with coeff 2 (:linear_const) → blocked until substitution
    ]
    red_eqs, red_obs, red_states = _reduce_equations(eqs, Equation[], [sx, sy_sd])
    @test length(red_states) == 1 && isequal(only(red_states), sx.val)
    @test length(red_eqs) == 1 && mtkext.isdifferential(only(red_eqs).lhs)
    @test length(red_obs) == 1 && isequal(only(red_obs).lhs, sy_sd.val)
    # sy_sd = -2*D(sx) evaluated at D(sx)=-sx → sy_sd = 2*sx
    @test isequal(only(red_obs).rhs, (2*sx).val)
end

@testset "FF-blocking in reduce_equations" begin
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
    red_eqs, red_obs, red_states = _reduce_equations(
        pq_eqs, Equation[], [u_r, u_i];
        outset=[u_r, u_i], ff_inputs=ff_set, verbose=false)
    @test length(red_states) == 2
    @test isempty(red_obs)
    @test length(red_eqs) == 2
    # Equations should be the originals (no divisions introduced)
    @test red_eqs == pq_eqs

    # B: Same equations, ff_inputs=Set() → normal solving (blocking disabled)
    red_eqs, red_obs, red_states = _reduce_equations(
        pq_eqs, Equation[], [u_r, u_i];
        outset=[u_r, u_i], ff_inputs=Set(), verbose=false)
    @test isempty(red_states)
    @test length(red_obs) == 2

    # C: Non-output state in input-dependent equation → still solved (blocking only for outputs)
    eqs_n = [0 ~ n - i_r*u_r - i_i*u_i]  # n is not an output
    red_eqs, red_obs, red_states = _reduce_equations(
        eqs_n, Equation[], [n];
        outset=[], ff_inputs=ff_set, verbose=false)
    @test isempty(red_states)
    @test length(red_obs) == 1
    @test isequal(only(red_obs).lhs, n.val)

    # D: Output state in input-free equation → still solved (no input dependency)
    eqs_free = [0 ~ u_r - p1*p2]
    red_eqs, red_obs, red_states = _reduce_equations(
        eqs_free, Equation[], [u_r];
        outset=[u_r], ff_inputs=ff_set, verbose=false)
    @test isempty(red_states)
    @test length(red_obs) == 1
    @test isequal(only(red_obs).lhs, u_r.val)

    # E: Cramer's rule for non-output 2×2 coupled states: verify solution correctness
    # and absence of spurious i_r or i_i denominators
    @variables x(t) y(t)
    eqs_2x2 = [0 ~ P + x*i_r + i_i*y, 0 ~ Q + x*i_i - i_r*y]
    red_eqs, red_obs, red_states = _reduce_equations(
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
       -P ~ u_r*i_r + i_i*u_i   # u_r, u_i are outputs → blocked
        0 ~ Q + u_r*i_i - i_r*u_i   # same
        0 ~ n - u_r - u_i            # n is non-output, depends only on states (not inputs)
    ]
    red_eqs, red_obs, red_states = _reduce_equations(
        eqs_mixed, Equation[], [u_r, u_i, n];
        outset=[u_r, u_i], ff_inputs=ff_set, verbose=false)
    @test length(red_states) == 2   # u_r, u_i stay as constraints
    @test length(red_obs) == 1      # n solved from its equation
    @test isequal(only(red_obs).lhs, n.val)

    # G: PV-like case: output blocked from linear P-eq, nonlinear V-eq stays
    @parameters Vset Pset
    eqs_pv = [
        0 ~ Pset + u_r*i_r + i_i*u_i   # linear in u_r, depends on inputs → blocked
        0 ~ -(Vset^2) + u_r^2 + u_i^2  # nonlinear → can't be solved for either state
    ]
    red_eqs, red_obs, red_states = _reduce_equations(
        eqs_pv, Equation[], [u_r, u_i];
        outset=[u_r, u_i], ff_inputs=ff_set, verbose=false)
    @test length(red_states) == 2   # both stay
    @test isempty(red_obs)
    @test length(red_eqs) == 2      # both original equations preserved intact

    # H: 4-SCC FF chain: input→LC→LS→LS→LC→output
    # cn1, cn2 are dynamic states whose values appear as coefficients,
    # making SCCs 2 and 3 :linear_state respectively.
    # Both SCC2 and SCC3 are :linear_state with equal cost; the index tiebreaker
    # selects SCC2 (lower equation index) as the forbidden match.
    # → SCC2 must be forbidden; SCCs 1, 3 and 4 must still be solved.
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
    red_eqs, red_obs, red_states = _reduce_equations(
        ff_chain_eqs, Equation[], all_chain_states;
        outset=[c4], ff_inputs=Set([c_inp]), verbose=false)
    obs_lhs = Set(eq.lhs for eq in red_obs)

    @test length(red_states) == 3
    @test red_eqs[3] ∈ ff_chain_eqs[4:5] # one of the state coeff should stae
    @test topologicical_sorted(red_obs)
end

@testset "diff-state-aware alias removal keeps differential representative" begin
    # Regression for a PowerDynamics machine+governor higher-index failure
    # (PSSE GENSAL + HYGOV): a differential speed state `w` is connected through a
    # transitive alias chain (w ~ spd_out ~ spd_in ~ speed_input) into governor
    # signals, and a lead block differentiates the chained signal (D(gerr) in an eq).
    # The first-pass alias removal must keep the *differential* state `w` as the
    # representative of the whole group, so every chained member becomes an
    # observation of `w` directly. The previous per-equation handling routed the
    # chain through a non-differential intermediate (e.g. speed_input), which
    # stranded the differential equation and forced a higher-index residual
    # (D(gerr) leaking into the RHS) later in the pipeline → RHSDifferentialsError.
    @variables w(t) spd_out(t) spd_in(t) speed_input(t) gerr(t) lead(t)
    @parameters nref Tr
    eqs = [
        Dt(w) ~ -w,                          # w is the differential state
        0 ~ spd_out - w,                     # SPEED_out.u ~ w   (pure alias)
        0 ~ spd_in - spd_out,                # connect           (pure alias)
        0 ~ speed_input - spd_in,            # speed_input ~ SPEED_in.u (pure alias)
        0 ~ gerr - (nref - speed_input),     # governor error (NOT a pure alias)
        0 ~ lead - (gerr + Tr*Dt(gerr)),     # lead block differentiates the chain
    ]
    states = [w, spd_out, spd_in, speed_input, gerr, lead]
    red_eqs, red_obs, red_states = _reduce_equations(eqs, Equation[], states; verbose=false)

    # the differential equation for w must survive the reduction
    @test any(eq -> mtkext.isdifferential(eq.lhs) && isequal(only(eq.lhs.args), w.val), red_eqs)
    # every alias member observes the differential inner `w` *directly* (star, not chain)
    obs = Dict(eq.lhs => eq.rhs for eq in red_obs)
    for s in (spd_out, spd_in, speed_input)
        @test haskey(obs, s.val)
        @test isequal(obs[s.val], w.val)
    end
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

    red_eqs, red_obs, red_states = _reduce_equations(
        eqs, Equation[], all_states;
        outset=outputs, ff_inputs=ff_set, verbose=false)

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

    @test mtkext.get_alias(Dt(a) ~ Dt(b)) == nothing
    @test mtkext.get_alias(Dt(a) ~ b) == nothing
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

@testset "initf metadata to init_formulas" begin
    @testset "no nesting" begin
        @component function slack_initf(; name)
            @parameters begin
                u_init_r = 1
                u_init_i = 0
            end
            @variables begin
                u_r(t), [initf = u_init_r + 0.1]
                u_i(t), [initf = u_init_i]
                i_r(t), [guess=0]
                i_i(t), [guess=0]
            end
            eqs = [Dt(u_r) ~ 0, Dt(u_i) ~ 0]
            System(eqs, t, [u_r, u_i, i_r, i_i], [u_init_r, u_init_i]; name)
        end
        @named sys = slack_initf()
        vm = VertexModel(sys, [:i_r, :i_i], [:u_r, :u_i]; verbose=false)

        @test has_initformula(vm)
        formulas = collect(get_initformulas(vm))
        @test length(formulas) == 2

        f_ur = only(filter(f -> f.outsym == [:u_r], formulas))
        f_ui = only(filter(f -> f.outsym == [:u_i], formulas))
        @test f_ur.sym == [:u_init_r]
        @test f_ui.sym == [:u_init_i]

        out = NetworkDynamics.SymbolicView(zeros(1), f_ur.outsym)
        f_ur(out, NetworkDynamics.SymbolicView([2.0], f_ur.sym))
        @test out[:u_r] ≈ 2.1

        out = NetworkDynamics.SymbolicView(zeros(1), f_ui.outsym)
        f_ui(out, NetworkDynamics.SymbolicView([3.5], f_ui.sym))
        @test out[:u_i] ≈ 3.5
    end

    @testset "parameter target" begin
        # a parameter is a valid initf target -- this is the case MTK `bindings` cannot express
        @component function param_target(; name)
            @variables begin
                x(t), [guess=0]
                y(t), [guess=0]
            end
            @parameters begin
                setp, [guess=0, initf = 3*x + 1]
            end
            System([Dt(x) ~ -x + setp, Dt(y) ~ -y], t, [x, y], [setp]; name)
        end
        @named sys = param_target()
        vm = VertexModel(sys, [], [:y]; verbose=false)
        f = only(get_initformulas(vm))
        @test f.outsym == [:setp]
        @test f.sym == [:x]

        out = NetworkDynamics.SymbolicView(zeros(1), f.outsym)
        f(out, NetworkDynamics.SymbolicView([2.0], f.sym))
        @test out[:setp] ≈ 7.0
    end

    @testset "nested model: metadata expressions are namespaced" begin
        # The whole point of `collect_initf` recursing over the *hierarchical* system:
        # `renamespace` renames a symbol but does not descend into its metadata, so a naive
        # read off the flattened system would leave `a` bare instead of `c₊a`.
        @component function initf_inner(; name)
            @parameters a=1.0
            @variables begin
                x(t), [guess=0, initf = 2*a]
                y(t), [guess=0]
            end
            System([Dt(x) ~ -x, Dt(y) ~ -y], t, [x, y], [a]; name)
        end
        @component function initf_outer(; name)
            @named c = initf_inner()
            @variables z(t), [guess=0]
            System([Dt(z) ~ -z + c.x], t, [z], []; name, systems=[c])
        end
        @named outer = initf_outer()
        vm = VertexModel(outer, [], [:z]; verbose=false)

        f = only(get_initformulas(vm))
        @test f.outsym == [:c₊x]
        @test f.sym == [:c₊a]   # <- namespaced, NOT bare `:a`

        out = NetworkDynamics.SymbolicView(zeros(1), f.outsym)
        f(out, NetworkDynamics.SymbolicView([4.0], f.sym))
        @test out[:c₊x] ≈ 8.0
    end

    @testset "nested model: foreign symbols passed in stay in the parent namespace" begin
        # `@component` promotes symbolic kwargs to `ParentScope`, so an expression handed
        # down from the parent must NOT be namespaced into the child.
        @component function initf_child(; name, foreign)
            @parameters a=1.0
            @variables begin
                x(t), [guess=0, initf = 2*a + foreign]
                y(t), [guess=0]
            end
            System([Dt(x) ~ -x, Dt(y) ~ -y], t, [x, y], [a]; name)
        end
        @component function initf_parent(; name)
            @parameters q=5.0
            @named c = initf_child(; foreign=3*q)
            @variables z(t), [guess=0]
            System([Dt(z) ~ -z + c.x], t, [z], [q]; name, systems=[c])
        end
        @named p = initf_parent()
        vm = VertexModel(p, [], [:z]; verbose=false)

        f = only(get_initformulas(vm))
        @test f.outsym == [:c₊x]
        @test Set(f.sym) == Set([:c₊a, :q])   # child symbol namespaced, parent symbol bare

        vals = Dict(:c₊a => 4.0, :q => 2.0)
        out = NetworkDynamics.SymbolicView(zeros(1), f.outsym)
        f(out, NetworkDynamics.SymbolicView([vals[s] for s in f.sym], f.sym))
        @test out[:c₊x] ≈ 2*4.0 + 3*2.0
    end

    @testset "conflicting formula targets" begin
        # two initf entries forcing the same raw target to different values is a compile
        # error; two entries on different members of one alias class both compile and are
        # reported by the init-time duplicate-writer check, once normalization has collapsed
        # the class (formulas are ejected raw, all classification happens at init).
        @component function conflict_inner(; name)
            @parameters a=1.0
            @variables begin
                x(t), [guess=0, initf = 2*a]
                y(t), [guess=0]
            end
            System([Dt(x) ~ -x, Dt(y) ~ -y], t, [x, y], [a]; name)
        end
        @component function conflict_outer(; name)
            @parameters b=2.0
            @named c = conflict_inner()
            @variables begin
                z(t), [guess=0]
                w(t), [guess=0, initf = 5*b]   # w aliases c₊x -> two writers, different values
            end
            eqs = [w ~ c.x, Dt(z) ~ -z + w + b]
            System(eqs, t; name, systems=[c])
        end
        @named outer = conflict_outer()
        vm = VertexModel(outer, [], [:z]; verbose=false)
        @test length(get_initformulas(vm)) == 2
        err = try; initialize_component(vm; verbose=false); catch e; e; end
        @test err isa ArgumentError
        @test occursin("Multiple InitFormulas set the same symbol", sprint(showerror, err))
    end
end

@testset "set_initf: system-level init formulas" begin
    # the backward-flow showcase: the child knows only its own inverse (initf on x reads
    # the child's own output y), the parent knows what that output must be in steady state
    # and pins it via set_initf. Together they determine everything — no nonlinear solve.
    @component function pin_pi_block(; name)
        @parameters K_p=20.0 K_i=5.0
        @variables begin
            err(t)
            y(t)
            x(t), [guess=0, initf=(y - K_p*err)/K_i]
        end
        System([Dt(x) ~ err, y ~ K_p*err + K_i*x], t; name)
    end
    @component function pin_ctrl(; name)
        @named pi = pin_pi_block()
        @parameters K=2.0 vref=1.0
        @variables begin
            v(t) = 1.0
            i(t), [input=true]
            o(t), [output=true]
        end
        eqs = [pi.err ~ vref - v, Dt(v) ~ pi.y - K*v + i, o ~ v]
        sys = System(eqs, t; name, systems=[pi])
        set_initf(sys, pi.y => K*v - i)
    end

    @testset "parent pins the child's observable output" begin
        @named ctrl = pin_ctrl()
        vm = VertexModel(ctrl, [:i], [:o]; verbose=false)
        @test :pi₊y ∈ obssym(vm)
        @test NetworkDynamics.pinned_obssyms(vm) == Set([:pi₊y])

        io = IOBuffer()
        state = initialize_component(vm;
            default_overrides=Dict(:i => 0.5, :o => 1.0), verbose=true, io)
        @test occursin("No free variables!", String(take!(io)))
        @test state[:pi₊x] ≈ 0.3    # ((K*v - i) - K_p*0) / K_i
        @test !haskey(state, :pi₊y)
    end

    @testset "pairs survive a second nesting level" begin
        @component function pin_outer(; name)
            @named ctrl = pin_ctrl()
            @variables begin
                iin(t), [input=true]
                oout(t), [output=true]
            end
            System([ctrl.i ~ iin, oout ~ ctrl.o], t; name, systems=[ctrl])
        end
        @named outer = pin_outer()
        vm = VertexModel(outer, [:iin], [:oout]; verbose=false)
        @test NetworkDynamics.pinned_obssyms(vm) == Set([:ctrl₊pi₊y])
        state = initialize_component(vm;
            default_overrides=Dict(:iin => 0.5, :oout => 1.0), verbose=false)
        @test state[:ctrl₊pi₊x] ≈ 0.3
    end

    @testset "settable targets work like the variable option" begin
        @component function settable_target(; name)
            @parameters a=1.0
            @variables begin
                x(t), [guess=0]
                y(t), [guess=0]
            end
            sys = System([Dt(x) ~ -x, Dt(y) ~ -y], t, [x, y], [a]; name)
            set_initf(sys, x => 2a)
        end
        @named sys = settable_target()
        vm = VertexModel(sys, [], [:y]; verbose=false)
        f = only(get_initformulas(vm))
        @test f.outsym == [:x]
        @test f.sym == [:a]
    end

    @testset "conflict with a variable-option initf on the same target errors" begin
        @component function initf_conflict(; name)
            @parameters a=1.0
            @variables begin
                x(t), [guess=0, initf=2a]
                y(t), [guess=0]
            end
            sys = System([Dt(x) ~ -x, Dt(y) ~ -y], t, [x, y], [a]; name)
            set_initf(sys, x => 3a)   # different recipe for the same raw target
        end
        @named sys = initf_conflict()
        @test_throws "conflicting definitions" VertexModel(sys, [], [:y]; verbose=false)
    end

    @testset "eager validation and appending" begin
        @variables x(t) y(t)
        @named sys = System([Dt(x) ~ -x, Dt(y) ~ -y], t, [x, y], [])
        @test_throws ArgumentError set_initf(sys, x + y => 1.0)
        sys2 = set_initf(set_initf(sys, x => 1.0), y => 2.0)  # appends, non-mutating
        @test length(mtkext.collect_initf(sys2)) == 2
        @test isempty(mtkext.collect_initf(sys))
    end

    # The same pin can be spelled without set_initf, through a parent-local alias variable
    # carrying the initf. This works deterministically: alias-name selection prefers the
    # shallower name, so the parent-local variable always becomes the class representative
    # holding the defining equation, the child symbol becomes an alias leaf, and the child's
    # read expands onto the pinned representative and stops there.
    @testset "alternative spelling: parent-local alias variable" begin
        @component function ctrl_aliasvar(; name)
            @named pi = pin_pi_block()
            @parameters K=2.0 vref=1.0
            @variables begin
                v(t) = 1.0
                i(t), [input=true]
                o(t), [output=true]
            end
            @variables y_wish(t), [initf = K*v - i]
            eqs = [y_wish ~ pi.y, pi.err ~ vref - v, Dt(v) ~ pi.y - K*v + i, o ~ v]
            System(eqs, t; name, systems=[pi])
        end
        @named ctrl = ctrl_aliasvar()
        vm = VertexModel(ctrl, [:i], [:o]; verbose=false)
        @test NetworkDynamics.pinned_obssyms(vm) == Set([:y_wish])
        state = initialize_component(vm;
            default_overrides=Dict(:i => 0.5, :o => 1.0), verbose=false)
        @test state[:pi₊x] ≈ 0.3
    end
end

@testset "set_guessf: system-level guess formulas" begin
    @testset "subsystem-owned target is namespaced" begin
        @component function guessf_inner(; name)
            @parameters a=1.0
            @variables x(t) [guess=0.0] y(t) [guess=0.0]
            System([Dt(x) ~ -x + a, Dt(y) ~ -y], t, [x, y], [a]; name)
        end
        @component function guessf_outer(; name)
            @named c = guessf_inner()
            @variables z(t) [guess=0.0]
            sys = System([Dt(z) ~ -z + c.x], t, [z], []; name, systems=[c])
            set_guessf(sys, c.x => c.a)
        end
        @named outer = guessf_outer()
        vm = VertexModel(outer, [], [:z]; verbose=false)
        f = only(get_guessformulas(vm))
        @test f.outsym == [:c₊x]
        @test f.sym == [:c₊a]   # <- namespaced, NOT bare `:a`
    end

    @testset "conflicting guessf targets warn, never error" begin
        # unlike set_initf (a constraint → error), two guess recipes for one target are only
        # a hint clash: dedupe with a warning and keep one.
        @component function guessf_conflict(; name)
            @parameters a=1.0
            @variables x(t) [guess=0.0, guessf=2a] y(t) [guess=0.0]
            sys = System([Dt(x) ~ -x, Dt(y) ~ -y], t, [x, y], [a]; name)
            set_guessf(sys, x => 3a)   # different recipe for the same raw target
        end
        @named sys = guessf_conflict()
        local vm
        @test_logs (:warn, r"conflicting definitions") match_mode=:any begin
            vm = VertexModel(sys, [], [:y]; verbose=false)
        end
        @test length(get_guessformulas(vm)) == 1
    end

    @testset "eager validation and appending" begin
        @variables x(t) y(t)
        @named sys = System([Dt(x) ~ -x, Dt(y) ~ -y], t, [x, y], [])
        @test_throws ArgumentError set_guessf(sys, x + y => 1.0)
        sys2 = set_guessf(set_guessf(sys, x => 1.0), y => 2.0)  # appends, non-mutating
        @test length(mtkext.collect_guessf(sys2)) == 2
        @test isempty(mtkext.collect_guessf(sys))
    end
end

@testset "symbolic bindings on unknowns are rejected" begin
    # `@variables x(t) = <symbolic>` is a binding; on an unknown that is ambiguous and the
    # user must say `initf` instead.
    @component function state_binding(; name)
        @parameters u_init_r=1
        @variables begin
            u_r(t) = u_init_r + 0.1     # <- binding on an unknown
            u_i(t)
            i_r(t), [guess=0]
            i_i(t), [guess=0]
        end
        System([Dt(u_r) ~ 0, Dt(u_i) ~ 0], t, [u_r, u_i, i_r, i_i], [u_init_r]; name)
    end
    @named s = state_binding()
    err = try
        VertexModel(s, [:i_r, :i_i], [:u_r, :u_i]; verbose=false)
        nothing
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("binds unknown(s)", err.msg)
    @test occursin("initf", err.msg)

    # explicit `bindings=` kwarg is the same thing and must be caught too
    @component function state_binding_kw(; name)
        @parameters p=1.0
        @variables begin
            x(t)
            y(t), [guess=0]
        end
        System([Dt(x) ~ -x, Dt(y) ~ -y], t, [x, y], [p]; name, bindings=[x => 2p])
    end
    @named s2 = state_binding_kw()
    @test_throws ArgumentError VertexModel(s2, [], [:y]; verbose=false)
end

@testset "symbolic defaults on parameters still shadow" begin
    # A symbolic default on a *parameter* is a runtime dependency, not an init formula: MTK
    # lowers it to an observed equation and the parameter disappears from the compiled model.
    @component function param_shadow(; name)
        @parameters begin
            p1 = 2.0
            p2 = 3*p1        # <- parameter binding: legal, reduces the parameter count
        end
        @variables begin
            x(t), [guess=0]
            y(t), [guess=0]
        end
        System([Dt(x) ~ -x + p2, Dt(y) ~ -y], t, [x, y], [p1, p2]; name)
    end
    @named s = param_shadow()
    vm = VertexModel(s, [], [:y]; verbose=false)   # must not throw
    @test :p1 ∈ psym(vm)
    @test :p2 ∉ psym(vm)          # shadowed away
    @test :p2 ∈ obssym(vm)        # ...into an observable
    @test !has_initformula(vm)    # and it is NOT an init formula
end

@testset "guessf variable option to guess_formulas" begin
    @component function vertex_with_guessf(; name)
        @parameters u_init=1.0
        @variables x(t) [guess=0.5] y(t) [guessf=2*u_init]
        eqs = [Dt(x) ~ -x; Dt(y) ~ -y]
        System(eqs, t, [x, y], [u_init]; name)
    end
    @named sys = vertex_with_guessf()
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

    # a scalar guess coexists as plain metadata (the fallback), NOT promoted to a formula
    @test all(gf -> :x ∉ gf.outsym, get_guessformulas(vm))
    @test get_guess(vm, :x) == 0.5
end

@testset "guessf seeds the init solver" begin
    # end-to-end: a free variable with only a `guessf` (no scalar guess) gets its solver
    # starting value from the formula.
    @component function seedme(; name)
        @parameters p=5.0
        @variables x(t) [guessf=p] o(t) [guess=0.0]
        System([Dt(x) ~ p - x, o ~ x], t, [x, o], [p]; name)
    end
    @named sys = seedme()
    vm = VertexModel(sys, Symbol[], [:o]; verbose=false)
    @test has_guessformula(vm)
    @test !has_guess(vm, :x)
    state = initialize_component(vm; verbose=false)
    @test state[:x] ≈ 5.0
end

@testset "guessf for eliminated variables survives as a raw formula" begin
    # `z` is algebraically eliminated (z = 2x) into a scaled-alias observable. The guess
    # formula is ejected raw, targeting `z` as written; at init time `normalize` transports
    # it onto the surviving symbol through the aliasmap. (This used to be skipped with a
    # warning when formulas were resolved symbolically at compile time.)
    @component function vertex_with_eliminated_guessf(; name)
        @parameters u_init = 1.0
        @variables x(t) y(t) [guess=1.0] z(t) [guessf=2*u_init]
        eqs = [Dt(x) ~ -x + z, 0 ~ z - 2 * x, Dt(y) ~ -y]
        System(eqs, t, [x, y, z], [u_init]; name)
    end
    @named sys = vertex_with_eliminated_guessf()
    vm = VertexModel(sys, [], [:x, :y]; verbose=false)

    f = only(get_guessformulas(vm))
    @test f.outsym == [:z]      # raw target, the alias key
    @test f.sym == [:u_init]
    n = NetworkDynamics.normalize(f, get_aliasmap(vm), vm)
    @test n.outsym == [:x]      # ... lands on the settable survivor at init time
    guesses = Dict{Symbol,Float64}()
    NetworkDynamics.apply_guess_formulas!(guesses, Dict(:u_init => 1.0), [n])
    @test guesses[:x] ≈ 1.0     # z = 2x = 2*u_init  =>  x = u_init

    # the constant guess for `y` is preserved as plain metadata
    @test has_guess(vm, :y)
end

@testset "is_symbolic_constant classification" begin
    @parameters p
    @variables x(t)
    # plain numbers and free-variable-free symbolics are constants (kept as :guess metadata)
    @test mtkext.is_symbolic_constant(2.0)
    @test mtkext.symbolic_constant_value(2.0) === 2.0
    @test mtkext.symbolic_constant_value(3) === 3
    # an MTKv11 literal guess is a constant BasicSymbolic (not isa Number) and must still be
    # folded to a Float64 (else such guesses are silently lost on assembly)
    c = Symbolics.unwrap(Symbolics.Num(0.0))
    @test !(c isa Number)
    @test mtkext.is_symbolic_constant(c)
    @test mtkext.symbolic_constant_value(c) === 0.0
    @test mtkext.is_symbolic_constant(Symbolics.Num(sqrt(2)))
    @test mtkext.symbolic_constant_value(Symbolics.Num(sqrt(2))) ≈ sqrt(2)
    # values referencing other variables/parameters are NOT constant → rejected as a scalar
    # `:guess`, must be spelled with the `guessf` option instead
    @test !mtkext.is_symbolic_constant(x)
    @test !mtkext.is_symbolic_constant(2x)
    @test !mtkext.is_symbolic_constant(2p)
end

@testset "guessf referencing nonexistent symbol is skipped, not fatal" begin
    # a guessf whose RHS references a symbol the compiled component does not expose is
    # malformed; the attach-time validation catches it and the lenient MTK attach demotes
    # the error to a warn-and-skip.
    @component function vertex_with_bad_guessf(; name)
        @variables foo(t) x(t) [guessf=2*foo] y(t)
        System([Dt(x) ~ -x, Dt(y) ~ -y], t, [x, y], []; name)
    end
    @named badsys = vertex_with_bad_guessf()
    local vm
    @test_logs (:warn, r"Skipping an? \w*Formula") match_mode=:any begin
        vm = VertexModel(badsys, [], [:x, :y]; verbose=false)
    end
    @test !has_guessformula(vm)
end

@testset "symbolic guess (non-constant) is rejected, directing to guessf" begin
    # the `guess` option is scalar-only; a symbolic guess must use `guessf`/`set_guessf`.
    @component function vertex_symguess_option(; name)
        @parameters q=3.0
        @variables x(t) [guess=q] y(t)
        System([Dt(x) ~ -x + q, Dt(y) ~ -y], t, [x, y], [q]; name)
    end
    @named s1 = vertex_symguess_option()
    @test_throws "guessf" VertexModel(s1, [], [:y]; verbose=false)

    @component function vertex_symguess_kwarg(; name)
        @parameters q=3.0
        @variables x(t) y(t)
        System([Dt(x) ~ -x + q, Dt(y) ~ -y], t, [x, y], [q]; name, guesses=[x => 2*q])
    end
    @named s2 = vertex_symguess_kwarg()
    @test_throws "guessf" VertexModel(s2, [], [:y]; verbose=false)
end

@testset "_dedupe_resolved conflict policy" begin
    # `_dedupe_resolved` collapses resolved (target, rhs) entries that share a target.
    # Comparison is on the symbolic rhs, so identical definitions dedupe silently, while
    # genuinely conflicting ones warn (guesses) or error (bindings) per the `fail` kw.
    @variables a(t) b(t)
    mk(target, rhs; weak=false) = (; src=rhs, target, rhs, input_symbolic=Any[], input_names=Symbol[], weak)

    # distinct targets: all kept
    @test length(mtkext._dedupe_resolved([mk(:x, a), mk(:y, b)]; fail=:warn, kind="G")) == 2

    # identical definitions for the same target: deduped to one, and NOT an error even
    # under fail=:error (equal rhs is a harmless alias-merge duplicate, not a conflict)
    @test length(mtkext._dedupe_resolved([mk(:x, a), mk(:x, a)]; fail=:error, kind="I")) == 1

    # a weak writer yields to a strong one on the same target (order-independent); all-weak
    # duplicates stay weak
    @test only(mtkext._dedupe_resolved([mk(:x, a; weak=true), mk(:x, a)]; fail=:error, kind="I")).weak == false
    @test only(mtkext._dedupe_resolved([mk(:x, a), mk(:x, a; weak=true)]; fail=:error, kind="I")).weak == false
    @test only(mtkext._dedupe_resolved([mk(:x, a; weak=true), mk(:x, a; weak=true)]; fail=:error, kind="I")).weak == true

    # weak yields to a strong writer even with a *differing* rhs — no conflict, the strong wins
    let r = only(mtkext._dedupe_resolved([mk(:x, a), mk(:x, b; weak=true)]; fail=:error, kind="I"))
        @test !r.weak && isequal(r.rhs, a)
    end
    # two *weak* writers with differing rhs stay a genuine conflict
    @test_throws "conflicting definitions" mtkext._dedupe_resolved([mk(:x, a; weak=true), mk(:x, b; weak=true)]; fail=:error, kind="I")

    # conflicting definitions (same target, differing rhs), fail=:warn → keep one + warn
    local kept
    @test_logs (:warn, r"conflicting definitions target x") match_mode=:any begin
        kept = mtkext._dedupe_resolved([mk(:x, a), mk(:x, b)]; fail=:warn, kind="G")
    end
    @test length(kept) == 1

    # conflicting definitions, fail=:error → throw
    @test_throws "conflicting definitions" mtkext._dedupe_resolved([mk(:x, a), mk(:x, b)]; fail=:error, kind="I")
end

@testset "getproperty_symbolic resolves flat ₊-named leaf variable" begin
    # a flattened/transpiled leaf can carry `₊` *inside* a single variable name
    # (e.g. `washout₊derivative₊x` is one unknown, not the chain washout→derivative→x).
    # getproperty_symbolic must NOT descend into a nonexistent subsystem `a`, but
    # resolve `a₊b` as one flat variable.
    n = Symbol("a₊b")
    xab = only(@variables $n(t))
    @named foo = System([Dt(xab) ~ -xab], t, [xab], [])
    foo = mtkcompile(foo)
    r = mtkext.getproperty_symbolic(foo, n; might_contain_toplevel_ns=false)
    @test isequal(r, Symbolics.unwrap(xab))
end

@testset "bound parameters become observed" begin
    # in MTK11, bound parameters (functions of other parameters only) become observed (no parameter anymore)
    # this i helpful because we can actually "assign" parameter equivalents across models
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
    vm1 = VertexModel(node, [:u], [:x]; mtkcompile=true)
    @test :Sn  ∉ Set(NetworkDynamics.psym(vm1))
    @test :S_b ∈ Set(NetworkDynamics.psym(vm1))
    @test :Sn ∈ Set(NetworkDynamics.obssym(vm1)) # sn is observed
    vm2 = VertexModel(node, [:u], [:x]; mtkcompile=false)
    @test :Sn  ∉ Set(NetworkDynamics.psym(vm2))
    @test :S_b ∈ Set(NetworkDynamics.psym(vm2))
    @test :Sn ∈ Set(NetworkDynamics.obssym(vm2)) # sn is observed
end

@testset "SimpleLead: Dt(input) requires index reduction" begin
    # Minimal reproduction of the HYGOV SimpleLead pattern.
    #
    # A SimpleLead block has the equation  T*Dt(in) ~ K*out - in
    # which is an improper transfer function (1+sT)/K — it uses the
    # derivative of its input.
    #
    # When composed as:  network_input → Filter → Lead → network_output
    # the Lead's Dt(in) = Dt(filter_x) is already defined by the Filter's
    # diff eq.  After substitution the Lead equation becomes purely algebraic
    # in lead_out:
    #   lead_out = (filter_x + T_l * (P - filter_x)/T_f) / K_l
    #
    # Currently mtkcompile classifies T_l*Dt(lead_in) ~ K_l*lead_out - lead_in
    # as a second differential equation for lead_in, leaving lead_out undefined
    # → Nstates ≠ Neqs.

    # Variant A: pre-aliased (filter_x used directly in both equations)
    # After flattening this gives two Dt(filter_x) equations and no eq for lead_out.
    @mtkmodel FilterLeadNode_A begin
        @variables begin
            filter_x(t) = 0.0, [description="filter state"]
            lead_out(t),       [description="lead output"]
            u(t), [description="output signal", output=true]
            P(t), [description="input signal", input=true]
        end
        @parameters begin
            T_f = 1.0, [description="filter time constant"]
            T_l = 1.0, [description="lead time constant"]
            K_l = 2.0, [description="lead gain"]
        end
        @equations begin
            Dt(filter_x) ~ (P - filter_x) / T_f
            T_l * Dt(filter_x) ~ K_l * lead_out - filter_x
            u ~ lead_out
        end
    end
    @named node_a = FilterLeadNode_A()
    v = VertexModel(node_a, [:P], [:u]; verbose=false)
    @test length(sym(v)) == 2 && v.mass_matrix == Diagonal([1, 0])

    v = VertexModel(node_a, [:P], [:u], ff_to_constraint=false, verbose=false)
    @test dim(v) == 1

    # Variant B: separate variables with alias connection (closer to real @components)
    # filter_out and lead_in are distinct variables connected by an alias equation.
    # filter_x has one diff eq, lead_in has another, and filter_out ~ lead_in ties them.
    # After alias resolution, the lead's diff eq should become algebraic in lead_out.
    @mtkmodel FilterLeadNode_B begin
        @variables begin
            filter_x(t) = 0.0, [description="filter state"]
            filter_out(t),     [description="filter output = filter_x"]
            lead_in(t),        [description="lead input, aliased to filter_out"]
            lead_out(t),       [description="lead output"]
            u(t), [description="output signal", output=true]
            P(t), [description="input signal", input=true]
        end
        @parameters begin
            T_f = 1.0, [description="filter time constant"]
            T_l = 1.0, [description="lead time constant"]
            K_l = 2.0, [description="lead gain"]
        end
        @equations begin
            # Filter block (SimpleLag flattened)
            Dt(filter_x) ~ (P - filter_x) / T_f
            filter_out ~ filter_x
            # Lead block (SimpleLead flattened) — Dt(in) is the problematic term
            T_l * Dt(lead_in) ~ K_l * lead_out - lead_in
            # Connection: lead input = filter output
            lead_in ~ filter_out
            # Network output
            u ~ lead_out
        end
    end
    @named node_b = FilterLeadNode_B()
    v = VertexModel(node_b, [:P], [:u], verbose=false)
    # filter_x is a diff state, u is algebraic constraint (FF prevention)
    @test length(sym(v)) == 2 && v.mass_matrix == Diagonal([1, 0])

    v = VertexModel(node_b, [:P], [:u], ff_to_constraint=false, verbose=false)
    @test length(sym(v)) == 1 && v.mass_matrix == Diagonal([1])
end

@testset "Algebraic loop with clamp/min (K_G feedback regression)" begin
    # Reproduces the esst4b₊ loop seen in real power system models.
    #
    # Five algebraic variables form a genuine loop via K_G feedback:
    #   ce  = va_out - K_G*efd          ← efd feeds back here
    #   vmp = K_PM * ce
    #   vmo = clamp(x_int + vmp, ...)   ← clamp makes the chain nonlinear in ce
    #   vml = min(vmo, VOEL)
    #   efd = vml * vb_signal           ← closes the loop
    #
    # Without clamp/min the loop is fully reduced (all 5 states solved).
    # With clamp/min the loop is partially reduced: vmp, vmo, vml, efd are solved
    # but ce stays as the residual constraint because substitution through clamp/min
    # makes the loop equation nonlinear in ce.

    @variables ce(t) vmp(t) vmo(t) vml(t) efd(t)
    @parameters K_G K_PM va_out x_int vb_signal V_MMIN V_MMAX VOEL
    all_states = [ce, vmp, vmo, vml, efd]

    # Baseline: pure-linear loop (no clamp/min) → fully solved
    eqs_linear = [
        0 ~ ce  - va_out + K_G*efd,
        0 ~ vmp - K_PM*ce,
        0 ~ vmo - (x_int + vmp),
        0 ~ vml - vmo,
        0 ~ efd - vml*vb_signal,
    ]
    red_eqs, red_obs, red_states = _reduce_equations(eqs_linear, Equation[], all_states, verbose=true)
    @test isempty(red_states)
    @test length(red_obs) == 5

    # Same loop with clamp/min: ce remains as the residual algebraic constraint
    eqs_with_clamp = [
        0 ~ ce  - va_out + K_G*efd,
        0 ~ vmp - K_PM*ce,
        0 ~ vmo - clamp(x_int + vmp, V_MMIN, V_MMAX),
        0 ~ vml - min(vmo, VOEL),
        0 ~ efd - vml*vb_signal,
    ]
    red_eqs, red_obs, red_states = _reduce_equations(eqs_with_clamp, Equation[], all_states, verbose=true)
    @test length(red_states) == 1
    @test length(red_obs) == 4
end

@testset "selective expander" begin
    @variables a b c d e target1 target2
    obseqs = [
        a ~ target1
        b ~ 4
        c ~ a + b
        d ~ b + target2
        e ~ b + d
    ]
    ex = mtkext.selective_expander(obseqs, [target1, target2])
    @test isequal(ex(a), target1)
    @test isequal(ex(b), b)
    @test isequal(ex(c), target1 + b)
    @test isequal(ex(d), b + target2)
    @test isequal(ex(e), b + (b + target2))
end

@testset "RHS derivative substitution after diff-obs move-back" begin
    # Regression test: when two diff equations are solved and moved to obs, and one
    # depends on D(x) of the other (e.g. governor lead-lag: T3*D(x2) = x1+T2*D(x1)-x2),
    # the final step that moves all diff-obs back to eqs must substitute D(x1) in the
    # RHS of D(x2).  Without the fix this would leave D(x1) in the RHS and cause a
    # RHSDifferentialsError downstream.
    @variables x1(t) x2(t)
    @parameters T1 T2 T3 ref
    eqs = [
        Dt(x1) ~ (ref - x1) / T1,            # explicit ODE for x1
        0 ~ x2 - x1 - T2*Dt(x1) + T3*Dt(x2), # D(x2) depends on D(x1) in solution
    ]
    red_eqs, red_obs, red_states = _reduce_equations(eqs, Equation[], [x1, x2])
    # Both states keep their own ODE
    @test length(red_states) == 2
    @test length(red_eqs) == 2
    @test all(mtkext.isdifferential(eq.lhs) for eq in red_eqs)
    # No RHS should contain a derivative after the fix
    @test isempty(mtkext.rhs_differentials(red_eqs))
end

@testset "implicit_output algebraic loop should not crash" begin
    # Regression test: PQ-load connected to a busbar where the busbar uses
    # implicit_output(i_r + i_i) to signal an output coupling. This creates a
    # symbolic cycle in the observation dependency graph:
    #   pq_tu_r → b_tu_r → implicit_output(b_i_r + b_i_i) → b_i_r
    #           → b_ti_r → pq_ti_r → pq_tu_r  (via pq current equations)
    # The cycle is not a real algebraic loop at runtime (implicit_output = 0),
    # but the symbolic solver must not crash — it should leave the coupling
    # equations as residual constraints.
    @variables begin
        b_i_r(t), b_i_i(t)
        b_P(t), b_Q(t), b_u_mag(t), b_u_arg(t), b_i_mag(t), b_i_arg(t)
        b_tu_r(t), b_tu_i(t), b_ti_r(t), b_ti_i(t)
        pq_tu_r(t), pq_tu_i(t), pq_ti_r(t), pq_ti_i(t)
        b_u_r(t), b_u_i(t)  # voltage inputs (ff_inputs, not states)
    end
    @parameters P Q

    io = implicit_output
    eqs = Equation[
        0 ~ b_P + b_i_i*(io(b_i_i+b_i_r)+b_u_i) + (b_u_r+io(b_i_i+b_i_r))*b_i_r,
        0 ~ b_Q - b_i_i*(b_u_r+io(b_i_i+b_i_r)) + (io(b_i_i+b_i_r)+b_u_i)*b_i_r,
        0 ~ b_u_mag - sqrt((b_u_r+io(b_i_i+b_i_r))^2 + (io(b_i_i+b_i_r)+b_u_i)^2),
        0 ~ b_u_arg - atan(io(b_i_i+b_i_r)+b_u_i, b_u_r+io(b_i_i+b_i_r)),
        0 ~ b_i_mag - sqrt(b_i_i^2+b_i_r^2),
        0 ~ b_i_arg - atan(b_i_i, b_i_r),
        0 ~ -b_tu_r + b_u_r + io(b_i_i+b_i_r),
        0 ~ -b_tu_i + io(b_i_i+b_i_r) + b_u_i,
        0 ~ -b_ti_r + b_i_r,
        0 ~ b_i_i - b_ti_i,
        0 ~ (-P*pq_tu_r - Q*pq_tu_i)/(pq_tu_r^2+pq_tu_i^2) + pq_ti_r,
        0 ~ (-P*pq_tu_i + Q*pq_tu_r)/(pq_tu_r^2+pq_tu_i^2) + pq_ti_i,
        0 ~ -pq_tu_r + b_tu_r,
        0 ~ b_tu_i - pq_tu_i,
        0 ~ b_ti_r + pq_ti_r,
        0 ~ b_ti_i + pq_ti_i,
    ]
    states = [b_i_r, b_i_i, b_P, b_Q, b_u_mag, b_u_arg, b_i_mag, b_i_arg,
              b_tu_r, b_tu_i, b_ti_r, b_ti_i, pq_tu_r, pq_tu_i, pq_ti_r, pq_ti_i]

    # Should not throw; the coupling equations form an irresolvable symbolic loop
    # and must remain as residual constraints.
    red_eqs, red_obs, red_states = _reduce_equations(eqs, Equation[], states;
                                                      outset=[b_i_r, b_i_i],
                                                      ff_inputs=[b_u_r, b_u_i])
    @test length(red_eqs) == 2
end

@testset "Test pre alias reduction to avoid large algebraic cycles over trivial equations" begin
    @mtkmodel SauerPaiMachine begin
        @parameters begin
            R_s=0, [description="stator resistance"]
            X_d=2.0, [description="d-axis synchronous reactance"]
            X_q=2.0, [description="q-axis synchronous reactance"]
            X′_d=0.3, [description="d-axis transient reactance"]
            X′_q=0.3, [description="q-axis transient reactance"]
            X″_d=0.2, [description="d-axis subtransient reactance"]
            X″_q=0.2, [description="q-axis subtransient reactance"]
            X_ls=0.172, [description="stator leakage reactance"]
            T′_d0=6.66667, [description="d-axis transient time constant"]
            T″_d0=0.075, [description="d-axis subtransient time constant"]
            T′_q0=6.66667, [description="q-axis transient time constant"]
            T″_q0=0.075, [description="q-axis subtransient time constant"]
            H=1.3, [description="inertia constant"]
            D=0, [description="direct shaft damping"]
            # System and machine base
            S_b=100, [description="System power basis in MVA"]
            V_b=110, [description="System voltage basis in kV"]
            ω_b=2*pi*50, [description="System base frequency in rad/s"]
            Sn=200, [description="Machine power rating in MVA"]
            Vn=120, [description="Machine voltage rating in kV"]
            # input/parameter switches
            vf, [guess=1, bounds=(0,Inf), description="field voltage"]
            τ_m, [guess=1, bounds=(0, Inf), description="mechanical torque"]
        end
        @variables begin
            # inputs
            u_r(t)=1, [description="bus d-voltage", output=true]
            u_i(t)=0, [description="bus q-voltage", output=true]
            i_r(t)=-0.5, [description="bus d-current (flowing into bus)", input=true]
            i_i(t)=0, [description="bus d-current (flowing into bus)", input=true]
            ψ_d(t), [description="d-axis flux linkage"]
            ψ_q(t), [description="q-axis flux linkage"]
            ψ″_d(t), [guess=-1, description="flux linkage assosciated with X″_d"]
            ψ″_q(t), [guess=0, description="flux linkage assosciated with X″_q"]
            I_d(t), [guess=0, description="d-axis current"]
            I_q(t), [guess=0, description="q-axis current"]
            V_d(t), [guess=0, description="d-axis voltage"]
            V_q(t), [guess=1, description="q-axis voltage"]
            E′_d(t), [guess=0, description="transient voltage behind transient reactance in d-axis"]
            E′_q(t), [guess=-1, description="transient voltage behind transient reactance in q-axis"]
            δ(t), [guess=0, description="rotor angle"]
            ω(t), [guess=1, description="rotor speed"]
            τ_e(t), [bounds=(0, Inf), description="electrical torque"]
            i_mag(t), [description="terminal current magnitude"]
            i_arg(t), [description="terminal current angle"]
        end
        begin
            γ_d1 = (X″_d - X_ls)/(X′_d - X_ls)
            γ_q1 = (X″_q - X_ls)/(X′_q - X_ls)
            γ_d2 = (X′_d-X″_d)/(X′_d-X_ls)^2 # ~ (1 - γ_d1)/(X′_d - X_ls)
            γ_q2 = (X′_q-X″_q)/(X′_q-X_ls)^2 # ~ (1 - γ_q1)/(X′_q - X_ls)
            T_to_loc(α)  = [ sin(α) -cos(α);
                            cos(α)  sin(α)]
            T_to_glob(α) = [ sin(α)  cos(α);
                            -cos(α)  sin(α)]
        end
        @equations begin
            [u_r, u_i] .~ T_to_glob(δ)*[V_d, V_q] * (Vn/V_b)
            [I_d, I_q] .~ -T_to_loc(δ)*[i_r, i_i] * ((S_b/V_b)/(Sn/Vn))

            τ_e ~ ψ_d*I_q - ψ_q*I_d
            Dt(δ) ~ ω_b*(ω - 1)
            2*H * Dt(ω) ~ τ_m  - τ_e - D*(ω - 1)
            V_d ~ -R_s*I_d - ω * ψ_q
            V_q ~ -R_s*I_q + ω * ψ_d

            T′_d0 * Dt(E′_q) ~ -E′_q - (X_d - X′_d)*(I_d - γ_d2*ψ″_d - (1-γ_d1)*I_d + γ_d2*E′_q) + vf
            T′_q0 * Dt(E′_d) ~ -E′_d + (X_q - X′_q)*(I_q - γ_q2*ψ″_q - (1-γ_q1)*I_q - γ_q2*E′_d)
            T″_d0 * Dt(ψ″_d) ~ -ψ″_d + E′_q - (X′_d - X_ls)*I_d
            T″_q0 * Dt(ψ″_q) ~ -ψ″_q - E′_d - (X′_q - X_ls)*I_q

            ψ_d ~ -X″_d*I_d + γ_d1*E′_q + (1-γ_d1)*ψ″_d
            ψ_q ~ -X″_q*I_q - γ_q1*E′_d + (1-γ_q1)*ψ″_q
            i_mag ~ i_r^2 + i_i^2
            i_arg ~ atan(i_i, i_r)
        end
    end

    @variables u_r(t) u_i(t) u_r_inter(t) u_i_inter(t)
    @variables i_r(t) i_i(t)

    # test compilation of 2 gen model
    @named gen1 = SauerPaiMachine()
    @named gen2 = SauerPaiMachine()
    eqs = [
        u_r ~ u_r_inter
        u_i ~ u_i_inter
        u_r_inter ~ gen1.u_r
        u_r_inter ~ gen2.u_r
        u_i_inter ~ gen1.u_i
        u_i_inter ~ gen2.u_i
        i_r ~ gen1.i_r + gen2.i_r
        i_i ~ gen1.i_i + gen2.i_i
    ]
    bus = System(eqs, t, [u_r, u_i, u_r_inter, u_i_inter, i_r, i_i], []; name=:bus, systems=[gen1, gen2])
    vm = VertexModel(bus, [:i_r, :i_i], [:u_r, :u_i]; verbose=false, check=false)
    @test dim(vm) == sum(vm.mass_matrix) + 2

    # test compilation of 3 gen model (same busbar equations, but more parallel machines)
    @named gen1 = SauerPaiMachine()
    @named gen2 = SauerPaiMachine()
    @named gen3 = SauerPaiMachine()
    eqs = [
        u_r ~ u_r_inter
        u_i ~ u_i_inter
        u_r_inter ~ gen1.u_r
        u_r_inter ~ gen2.u_r
        u_r_inter ~ gen3.u_r
        u_i_inter ~ gen1.u_i
        u_i_inter ~ gen2.u_i
        u_i_inter ~ gen3.u_i
        i_r ~ gen1.i_r + gen2.i_r + gen3.i_r
        i_i ~ gen1.i_i + gen2.i_i + gen3.i_i
    ]
    bus = System(eqs, t, [u_r, u_i, u_r_inter, u_i_inter, i_r, i_i], []; name=:bus, systems=[gen1, gen2, gen3])
    vm = VertexModel(bus, [:i_r, :i_i], [:u_r, :u_i]; verbose=false, check=false)
    @test dim(vm) == sum(vm.mass_matrix) + 2

    @testset "test effectiveness of scaling pseudoalg" begin
        genp = (S_b=100, V_b=1, ω_b=2π*50, R_s=0.000125, T″_d0=0.01, T″_q0=0.01, X_ls=0.01460, X_d=0.1460, X′_d=0.0608, X″_d=0.06, X_q=0.1000, X′_q=0.0969, X″_q=0.06, T′_d0=8.96, T′_q0=0.310, H=23.64, Sn=100, Vn=1)
        @named gen = SauerPaiMachine(; genp...)

        # init problem as current osurce model
        vm = VertexModel(gen, [:u_r, :u_i], [:i_r, :i_i], ff_to_constraint=false) # test that the machine model itself is consistent
        set_default!(vm, :u_r, 1.0)
        set_default!(vm, :u_i, 0.0)
        set_default!(vm, :i_r, -0.45)
        set_default!(vm, :i_i, 0.1)
        initialize_component(vm);
    end
end

@testset "AliasMap extraction" begin
    @testset "_match_scaled_var acceptance" begin
        @variables x(t) y(t)
        @parameters k V
        # pure scaled aliases of a single variable are accepted...
        @test mtkext._match_scaled_var(x)      == (1.0, :x)
        @test mtkext._match_scaled_var(-x)     == (-1.0, :x)
        @test mtkext._match_scaled_var(2.5x)   == (2.5, :x)
        @test mtkext._match_scaled_var(-2.5x)  == (-2.5, :x)
        @test mtkext._match_scaled_var(x/2)    == (0.5, :x)
        @test mtkext._match_scaled_var(V)      == (1.0, :V) # parameters are settable, so they alias
        # ...everything else stays an ordinary observable
        @test mtkext._match_scaled_var(x + 1)  === nothing # affine
        @test mtkext._match_scaled_var(x + y)  === nothing # sum
        @test mtkext._match_scaled_var(2x + 3y) === nothing
        @test mtkext._match_scaled_var(x^2)    === nothing # nonlinear
        @test mtkext._match_scaled_var(x*y)    === nothing
        @test mtkext._match_scaled_var(sin(x)) === nothing
        @test mtkext._match_scaled_var(k*x)    === nothing # symbolic factor
        @test mtkext._match_scaled_var(0*x)    === nothing # degenerate
        @test mtkext._match_scaled_var(Num(1.0)) === nothing # no variable at all
    end

    @testset "extraction from compiled component" begin
        @mtkmodel AliasTestBus begin
            @variables begin
                u_r(t) = 1.0
                θ(t)
                scaled(t)
                nl(t)
                Vmeas(t)
                i_r(t), [input=true]
                P(t), [output=true]
            end
            @parameters begin
                Vset = 1.0
            end
            @equations begin
                Dt(u_r) ~ -u_r + i_r
                θ ~ -u_r        # sign flipped alias of a state
                scaled ~ 2*θ    # chain: scaled = 2*θ = -2*u_r
                nl ~ u_r^2      # not an alias
                Vmeas ~ Vset    # alias of a parameter
                P ~ u_r + nl
            end
        end
        @named atb = AliasTestBus()
        vm = VertexModel(atb, [:i_r], [:P])

        am = get_aliasmap(vm)
        @test am == AliasMap(:θ      => (-1.0, :u_r),
                             :scaled => (-2.0, :u_r), # factors multiply along the chain
                             :Vmeas  => (1.0, :Vset))
        @test !haskey(am, :nl)

        # structural invariants: aliases are observables, roots are settable
        settable = NetworkDynamics.settable_symbols(vm)
        @test all(k -> k ∈ NetworkDynamics.obssym(vm), keys(am))
        @test all(k -> k ∉ settable, keys(am))
        @test all(v -> v[2] ∈ settable, values(am))

        # the extracted factors must agree with what the obs function actually computes
        set_default!(vm, :i_r, 0.5)
        set_default!(vm, :u_r, 0.8)
        for (alias, (factor, canonical)) in am
            @test get_initial_state(vm, alias) ≈ factor * get_initial_state(vm, canonical)
        end
    end

    @testset "no aliases in ordinary components" begin
        # observables of the library components are all genuinely computed quantities;
        # this pins the empty-aliasmap (backwards compatible) path
        for c in [Lib.dqbus_pv(), Lib.dqbus_swing(), Lib.swing_mtk(), Lib.line_mtk()]
            @test isempty(get_aliasmap(c))
        end
    end
end
