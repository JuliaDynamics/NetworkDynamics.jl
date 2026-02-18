using NetworkDynamics, Graphs
using SteadyStateDiffEq, OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqNonlinearSolve
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as Dt
using SymbolicIndexingInterface: SymbolicIndexingInterface as SII
using Chairmarks
using NonlinearSolve
using SciMLLogging
using NetworkDynamics: SubtaskError

@__MODULE__()==Main ? includet(joinpath(pkgdir(NetworkDynamics), "test", "ComponentLibrary.jl")) : (const Lib = Main.Lib)

@testset "test find_fixpoint" begin
    function swing_equation!(dv, v, esum, (M, P, D), t)
        dv[1] = v[2]
        dv[2] = 1/M *(P - D * v[2] + esum[1])
        nothing
    end
    swing_vertex = VertexModel(f=swing_equation!, g=1, sym=[:θ, :ω], psym=[:M=>1, :P, :D=>0.1])

    function powerflow!(e, v_s, v_d, (K,), t)
        e[1] = K * sin(v_s[1] - v_d[1])
    end
    powerflow_edge = EdgeModel(g=AntiSymmetric(powerflow!); outdim=1, psym=[:K=>6])
    g = watts_strogatz(4, 2, 0.0)

    nd = Network(g, swing_vertex, powerflow_edge)
    p = NWParameter(nd)
    p.v[1:4, :P] = [1, -1, 1, -1]

    s0 = find_fixpoint(nd, p)
    @test pflat(s0) == pflat(p)

    s1 = find_fixpoint(nd, p; alg=DynamicSS(Rodas5P()))
    @test isapprox(diff(s0.v[1:nv(g),:θ]), diff(s1.v[1:nv(g),:θ]); atol=1e-8)
end

@testset "test component initialization" begin
    @mtkmodel InitSwing begin
        @variables begin
            u_r(t)=1, [description="bus d-voltage", output=true]
            u_i(t)=0.1, [description="bus q-voltage", output=true]
            i_r(t)=1, [description="bus d-current (flowing into bus)", input=true]
            i_i(t)=0.1, [description="bus d-current (flowing into bus)", input=true]
            ω(t), [guess=0.0, description="Rotor frequency"]
            θ(t), [guess=0.0, bounds=[-π, π], description="Rotor angle"]
            Pel(t), [guess=1, description="Electrical Power injected into the grid"]
        end
        @parameters begin
            M=0.005, [description="Inertia"]
            D=0.1, [description="Damping"]
            V, [description="Voltage magnitude", guess=1, bounds=(0,Inf)]
            ω_ref=0, [description="Reference frequency"]
            Pm, [guess=0.1,description="Mechanical Power"]
            aux, [description="Auxiliary, unused parameter"]
        end
        @equations begin
            Dt(θ) ~ ω - ω_ref
            Dt(ω) ~ 1/M * (Pm - D*ω - Pel)
            Pel ~ u_r*i_r + u_i*i_i
            u_r ~ V*cos(θ)
            u_i ~ V*sin(θ)
        end
    end
    sys = InitSwing(name=:swing)
    vf = VertexModel(sys, [:i_r, :i_i], [:u_r, :u_i])

    @test vf.symmetadata[:u_r][:default] == 1
    @test vf.symmetadata[:u_i][:default] == 0.1
    @test vf.symmetadata[:Pm][:guess] == 0.1
    @test vf.symmetadata[:θ][:bounds] == [-π, π]
    @test vf.symmetadata[:i_r][:default] == 1
    @test vf.symmetadata[:i_i][:default] == 0.1
    @test filter(p->NetworkDynamics.is_unused(vf, p), psym(vf)) == [:aux]
    @test NetworkDynamics.has_default_or_init(vf, :aux) == false

    NetworkDynamics.initialize_component!(vf; verbose=true)
    @test NetworkDynamics.init_residual(vf) < 1e-8

    @test_throws ComponentInitError NetworkDynamics.initialize_component!(vf, tol=0.0; verbose=true)

    # make empty problem (no free vars)
    vf_comp = copy(vf)
    set_default!(vf_comp, :Pm, get_init(vf_comp, :Pm))
    set_default!(vf_comp, :θ, get_init(vf_comp, :θ))
    set_default!(vf_comp, :ω, get_init(vf_comp, :ω))
    set_default!(vf_comp, :V, get_init(vf_comp, :V))

    NetworkDynamics.initialize_component!(vf_comp; verbose=true)
    dump_initial_state(vf)
    dump_initial_state(vf_comp)

    # check with wrong default/ on observed
    vf_def = copy(vf)
    set_default!(vf_def, :Pel, -100)
    set_bounds!(vf_def, :Pel, (-Inf, 0))

    buf = IOBuffer()
    NetworkDynamics.initialize_component!(vf_def; verbose=false, io=buf)
    out = String(take!(buf))
    @test contains(out, r"has broken bounds")
    @test contains(out, r"has observables that differ")

    vf_conflict = copy(vf)
    defaults = NetworkDynamics.get_defaults_dict(vf_conflict)
    # Modify the defaults dictionary to create conflicting values
    defaults[:u_r] = 0.0
    defaults[:u_i] = 1.0

    # Test that the metadata synchronization works with conflicting defaults
    NetworkDynamics.initialize_component!(vf_conflict;
        defaults=defaults,
        verbose=true)

    # Verify the metadata was updated and values were applied correctly
    @test get_default(vf_conflict, :u_r) == 0.0
    @test get_default(vf_conflict, :u_i) == 1.0
    @test get_init(vf_conflict, :V)  ≈ 1
    @test get_init(vf_conflict, :θ) % 2pi ≈ pi/2

    # change metadtaa by providing custom input
    vf_sync = copy(vf)
    custom_defaults = NetworkDynamics.get_defaults_dict(vf_sync)
    custom_guesses = NetworkDynamics.get_guesses_dict(vf_sync)
    custom_bounds = NetworkDynamics.get_bounds_dict(vf_sync)

    delete!(custom_defaults, :u_r) # removed
    custom_defaults[:θ] = 0.0      # additional
    custom_defaults[:u_i] = 0.0      # changed

    delete!(custom_guesses, :θ)    # removed
    custom_guesses[:u_r] = 0.1     # additional
    custom_guesses[:ω] = 1.1       # changed

    delete!(custom_bounds, :θ)       # removed
    custom_bounds[:u_r] = (0.0, 1.1) # additional

    # Initialize with custom metadata
    NetworkDynamics.initialize_component!(vf_sync;
        defaults=custom_defaults,
        guesses=custom_guesses,
        bounds=custom_bounds,
        verbose=true
    )

    # Direct comparison to verify metadata synchronization
    @test NetworkDynamics.get_defaults_dict(vf_sync) == custom_defaults
    @test NetworkDynamics.get_guesses_dict(vf_sync) == custom_guesses
    @test NetworkDynamics.get_bounds_dict(vf_sync) == custom_bounds
end

@testset "test component initialization with bounds" begin
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
            [u_r, u_i] .~ T_to_glob(δ)*[V_d, V_q] * Vn/V_b
            [I_d, I_q] .~ -T_to_loc(δ)*[i_r, i_i] * (S_b/V_b)/(Sn/Vn)

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

    sys = SauerPaiMachine(name=:swing)
    vf = VertexModel(sys, [:i_r, :i_i], [:u_r, :u_i])

    @test get_initial_state(vf, sym(vf)) == [nothing, nothing, nothing, nothing, nothing, nothing, 1.0, 0.0]
    @test get_initial_state(vf, insym(vf)) == [-0.5, 0.0]
    @test get_initial_state(vf, outsym(vf)) == [1.0, 0.0]
    get_initial_state(vf, psym(vf)) # evaluates
    get_initial_state(vf, obssym(vf)) # no error

    NetworkDynamics.initialize_component!(vf; verbose=true)
    @test get_initial_state(vf, :vf) > 1
    @test !any(isnothing, get_initial_state(vf, obssym(vf)))

    # known bad seed which validates the constraint
    set_guess!(vf, :E′_d,  -0.31205)
    set_guess!(vf, :E′_q,   1.3396)
    set_guess!(vf, :δ,      0.4516)
    set_guess!(vf, :ψ″_d,   0.50574)
    set_guess!(vf, :ψ″_q,   1.353)
    set_guess!(vf, :ω,     -1.55)

    buf = IOBuffer()
    NetworkDynamics.initialize_component!(vf; verbose=false, apply_bound_transformation=false, io=buf)
    out = String(take!(buf))
    @test contains(out, r"has broken bounds")

    @test get_initial_state(vf, :vf) < 1 # does not conserve
    NetworkDynamics.initialize_component!(vf; verbose=true, apply_bound_transformation=true)
    @test get_initial_state(vf, :vf) > 1 # conserves

    # no lets try to force it negative
    set_bounds!(vf, :vf, (-Inf, 0))
    set_guess!(vf, :vf, -1)
    NetworkDynamics.initialize_component!(vf; verbose=true)
    @test get_initial_state(vf, :vf) < 1

    # check performance
    prob, _ = NetworkDynamics.initialization_problem(vf,
        get_defaults_dict(vf), get_guesses_dict(vf), get_bounds_dict(vf);
        apply_bound_transformation=true);
    du = zeros(length(prob.f.resid_prototype));
    b = @b $(prob.f)($du, $(prob.u0), nothing)
    @test iszero(b.allocs)
end

@testset "Test edge initialization" begin
    @mtkmodel StaticPowerLine begin
        @variables begin
            srcθ(t), [description = "voltage angle at src end", input=true]
            dstθ(t), [description = "voltage angle at dst end", input=true]
            P(t), [description = "flow towards node at dst end", output=true]
            Δθ(t)
        end
        @parameters begin
            active = 1, [description = "line active"]
            K = 1.63, [description = "Line conductance"]
            limit = 1, [description = "Line limit"]
        end
        @equations begin
            Δθ ~ srcθ - dstθ
            P ~ active*K*sin(Δθ)
        end
    end
    em = EdgeModel(StaticPowerLine(name=:line_mtk), [:srcθ], [:dstθ], AntiSymmetric([:P]))

    state = initialize_component(em;
        defaults=Dict(:K=>1.63, :active=>1, :srcθ=>0.1, :dstθ=>0),
        guesses=Dict(:P=>0.0, :₋P=>0.0),
        verbose=true)

    @test iszero(sum(get_initial_state(em, state, [:P, :₋P])))
    @test iszero(init_residual(em, state))

    # test init_residual with missing state:
    state_minimal = Dict(
        :P => state[:P],
        :₋P =>state[:₋P],
        :K => state[:K],
        :active => state[:active],
        :dstθ => state[:dstθ],
        :srcθ=> state[:srcθ],
    )
    @test init_residual(em, state_minimal) == 0
    delete!(state_minimal, :P)
    @test_throws ArgumentError init_residual(em, state_minimal)

    state2 = initialize_component(em;
        default_overrides=Dict(:srcθ=>0.0, :dstθ=>0.1),
        guess_overrides=Dict(:P=>0.0, :₋P=>0.0),
        verbose=true)
    @test state2[:P]  == -state[:P]
    @test state2[:₋P] == -state[:₋P]

    em_mut = copy(em)
    initialize_component!(em_mut;
        default_overrides=Dict(:srcθ=>0.0, :dstθ=>0.1),
        guess_overrides=Dict(:P=>0.0, :₋P=>0.0),
        verbose=true)
    @test get_guess(em_mut, :P) == 0.0
    @test get_guess(em_mut, :₋P) == 0.0
    @test get_default(em_mut, :srcθ) == 0.0
    @test get_default(em_mut, :dstθ) == 0.1
end

@testset "Initialization constraint construction" begin
    using NetworkDynamics: @initconstraint, InitConstraint
    ic1 = @initconstraint :x + :y
    ic2 = InitConstraint([:x, :y], 1) do out, u
        out[1] = u[:x] + u[:y]
    end
    out1 = [0.0]
    out2 = [0.0]
    u = rand(2)
    ic1(out1, u)
    ic2(out2, u)
    @test out1 == out2

    ic1 = @initconstraint begin
        :x + :y
        :z^2
    end
    ic2 = InitConstraint([:x, :y, :z], 2) do out, u
        out[1] = u[:x] + u[:y]
        out[2] = u[:z]^2
    end
    out1 = [0.0, 0.0]
    out2 = [0.0, 0.0]
    u = rand(3)
    ic1(out1, u)
    ic2(out2, u)
    @test out1 == out2
end

@testset "InitConstraint combining constructor" begin
    using NetworkDynamics: @initconstraint, InitConstraint, dim

    # Test basic combining of two constraints
    c1 = @initconstraint :x + :y
    c2 = @initconstraint :z^2 - 1
    combined = InitConstraint(c1, c2)

    @test dim(combined) == 2
    @test Set(combined.sym) == Set([:x, :y, :z])

    # Test that combined constraint works correctly
    u = [1.0, 2.0, 3.0]  # x=1, y=2, z=3
    res_combined = zeros(2)
    res_individual = zeros(2)

    combined(res_combined, u)
    c1(view(res_individual, 1:1), u[1:2])
    c2(view(res_individual, 2:2), u[3])

    @test res_combined == res_individual
    @test res_combined[1] ≈ 3.0  # x + y = 1 + 2
    @test res_combined[2] ≈ 8.0  # z^2 - 1 = 9 - 1

    # check perfomance
    u_sview = SymbolicView(u, combined.sym)
    b = @b $(combined.f)($res_combined, $u_sview)
    @test b.allocs == 0

    # Test combining constraints with overlapping symbols
    c3 = @initconstraint :x - :y
    c4 = @initconstraint :x^2 + :w
    combined2 = InitConstraint(c3, c4)

    @test dim(combined2) == 2
    @test Set(combined2.sym) == Set([:x, :y, :w])

    # Test combining constraints with different dimensions
    c5 = @initconstraint begin
        :a + :b
        :c - :d
    end
    c6 = @initconstraint :e^2
    combined3 = InitConstraint(c5, c6)

    @test dim(combined3) == 3
    @test Set(combined3.sym) == Set([:a, :b, :c, :d, :e])

    u2 = [1.0, 2.0, 3.0, 4.0, 5.0]  # a=1, b=2, c=3, d=4, e=5
    res = zeros(3)
    combined3(res, u2)
    @test res[1] ≈ 3.0   # a + b = 1 + 2
    @test res[2] ≈ -1.0  # c - d = 3 - 4
    @test res[3] ≈ 25.0  # e^2 = 25

    # Test error case: empty varargs
    @test_throws ArgumentError InitConstraint()

    # Test single constraint (should work)
    single = InitConstraint(c1)
    @test single === c1

    # Test combining multiple constraints
    c7 = @initconstraint :p
    c8 = @initconstraint :q
    c9 = @initconstraint :r
    multi_combined = InitConstraint(c7, c8, c9)

    @test dim(multi_combined) == 3
    @test Set(multi_combined.sym) == Set([:p, :q, :r])
    u = [-1, 2, 5]
    res = zeros(3)
    multi_combined(res, u)
    @test res == u

    # test error on int indices
    c10 = InitConstraint([:z, :x, :y], 1) do out, u
        out[1] = u[1] + 2*u[2] + 3*u[3]  # u[1]=z, u[2]=x, u[3]=y
    end
    c11 = @initconstraint :z
    @test_throws ArgumentError InitConstraint(c10, c11)
end

_recursive_replace(el::Symbol, dict) = haskey(dict, el) ? dict[el] : NaN
_recursive_replace(containter, dict) = map(el -> _recursive_replace(el, dict), containter)
@testset "test input mapping" begin
    using NetworkDynamics: insym_normalized, outsym_normalized
    vm = Lib.swing_mtk()
    c = @initconstraint begin
        :Pdamping # observable
        :P # input
        :θ # output
        :ω # state
        :Pmech # parameter
    end
    dict = Dict(:Pdamping => 1, :P=> 2, :θ => 3, :ω => 4, :Pmech => 5)

    outdata = _recursive_replace(outsym_normalized(vm), dict)
    udata   = _recursive_replace(sym(vm), dict)
    indata  = _recursive_replace(insym_normalized(vm), dict)
    pdata   = _recursive_replace(psym(vm), dict)
    obsdata = _recursive_replace(obssym(vm), dict)

    mapping! = NetworkDynamics.generate_init_input_mapping(vm, c)

    syms = zeros(5)
    mapping!(syms, outdata, udata, indata, pdata, obsdata)
    @test syms == 1:5

    em = Lib.line_mtk()
    c = @initconstraint begin
        :P # output
        :dstθ # input
        :active #param
        :Δθ # obs
    end
    dict = Dict(:P => 1, :dstθ => 2, :active => 3, :Δθ => 4)

    outdata = _recursive_replace(outsym_normalized(em), dict)
    udata   = _recursive_replace(sym(em), dict)
    indata  = _recursive_replace(insym_normalized(em), dict)
    pdata   = _recursive_replace(psym(em), dict)
    obsdata = _recursive_replace(obssym(em), dict)

    mapping! = NetworkDynamics.generate_init_input_mapping(em, c)
    syms = zeros(4)
    mapping!(syms, outdata, udata, indata, pdata, obsdata)
    @test syms == 1:4
end

@testset "init with additonal contraints" begin
    @mtkmodel InitSwing begin
        @variables begin
            u_r(t), [description="bus d-voltage", output=true, guess=1]
            u_i(t)=1, [description="bus q-voltage", output=true, guess=1]
            i_r(t)=1, [description="bus d-current (flowing into bus)", input=true]
            i_i(t)=0.1, [description="bus d-current (flowing into bus)", input=true]
            u_mag(t), [description="bus voltage magnitude"]
            ω(t), [guess=0.0, description="Rotor frequency"]
            θ(t), [guess=0.0, bounds=[-π, π], description="Rotor angle"]
            Pel(t), [guess=1, description="Electrical Power injected into the grid"]
        end
        @parameters begin
            M=0.005, [description="Inertia"]
            D=0.1, [description="Damping"]
            V, [description="Voltage magnitude", guess=1]
            ω_ref=0, [description="Reference frequency"]
            Pm, [guess=0.1,description="Mechanical Power"]
            aux, [description="Auxiliary, unused parameter"]
        end
        @equations begin
            Dt(θ) ~ ω - ω_ref
            Dt(ω) ~ 1/M * (Pm - D*ω - Pel)
            Pel ~ u_r*i_r + u_i*i_i
            u_r ~ V*cos(θ)
            u_i ~ V*sin(θ)
            u_mag ~ sqrt(u_r^2 + u_i^2)
        end
    end
    vm = VertexModel(InitSwing(name=:swing), [:i_r, :i_i], [:u_r, :u_i])
    ic = @initconstraint begin
        :u_mag - 1
        :Pel - 0.1
    end
    res1 = initialize_component(vm, additional_initconstraint=ic, verbose=false)
    set_initconstraint!(vm, ic)
    res2 = initialize_component(vm, verbose=false)
    @test res1 == res2

    initialize_component!(vm)
    @test get_initial_state(vm, :u_mag) ≈ 1
    @test get_initial_state(vm, :Pel) ≈ 0.1

    # test if the whole observable stuff is allocation free
    prob, _ = NetworkDynamics.initialization_problem(
        vm, get_guesses_dict(vm), get_defaults_dict(vm), get_bounds_dict(vm),
        InitConstraint(get_initconstraints(vm)...);
        apply_bound_transformation=true);
    du = zeros(length(prob.f.resid_prototype));
    b = @b $(prob.f)($du, $(prob.u0), nothing)
    @test iszero(b.allocs)

    @test_throws ArgumentError set_initconstraint!(vm, @initconstraint :wrong_symbol)

    em = Lib.line_mtk()
    delete_default!(em, :active)
    set_initconstraint!(em, @initconstraint begin
        :Δθ - 1 # observable
        :srcθ
        :P # force P to zero
    end)
    initialize_component!(em, guess_overrides=Dict(:P=>0,:₋P=>0,:srcθ=>π,:dstθ=>-π,:active=>1.0), verbose=true)
    @test get_initial_state(em, :Δθ) ≈ 1
    @test get_initial_state(em, :dstθ) ≈ -1
    @test get_initial_state(em, :P) ≈ 0 atol=1e-10
    @test get_initial_state(em, :active) ≈ 0 atol=1e-10
end

@testset "multistep init of powergrid like network" begin
    # Create a simple network with Kuramoto oscillators
    g = cycle_graph(5) # 5-node cycle graph
    v1s = Lib.dqbus_slack()
    v2s = Lib.dqbus_pv(Pset=1.5, Vset=1.0)
    v3s = Lib.dqbus_pq(Pset=-1.0, Qset=-0.1)
    v4s = Lib.dqbus_pq(Pset=-1.0, Qset=-0.1)
    v5s = Lib.dqbus_pq(Pset=-1.0, Qset=-0.1)
    e = Lib.dqline(X=0.1, R=0.01)

    nws = Network(g, [v1s, v2s, v3s, v4s, v5s], e)
    pf = find_fixpoint(nws)

    v1 = Lib.dqbus_swing_and_load()
    set_initconstraint!(v1, @initconstraint begin
        :load₊Pinj + 1.0
        - :u_r^2 - :u_i^2 + :swing₊V^2
    end)
    v2 = Lib.dqbus_swing()
    v3 = Lib.dqbus_pq()
    v4 = Lib.dqbus_pq()
    v5 = Lib.dqbus_pq()
    nw = Network(g, [v1, v2, v3, v4, v5], e; dealias=true)

    default_overrides = merge(
        interface_values(pf),
        Dict(eidxs(1:5, :active) .=> nothing))

    s_nmut = initialize_componentwise(nw; subverbose=true, verbose=true, default_overrides)
    s_mut = initialize_componentwise!(nw; subverbose=true, verbose=true, default_overrides)
    s_parallel = initialize_componentwise(nw; subverbose=true, verbose=true, default_overrides, parallel=true)
    s_meta = NWState(nw)

    for (k, v) in interface_values(pf)
        @test s_nmut[k] ≈ v atol=1e-10
        @test s_parallel[k] ≈ v atol=1e-10
    end
    @test s_meta[VIndex(1, :swing₊V)] ≈ s_meta[VIndex(1, :u_mag)]
    @test s_meta[VIndex(1, :load₊Pset)] ≈ -1.0

    @testset "test passing of subalg and subsolve_kwargs" begin
        @test_throws SciMLBase.NonSolverError initialize_componentwise(nw; subalg=:foo)
        @test_throws SciMLBase.NonSolverError initialize_componentwise(nw; subalg=Dict(VIndex(:swing_and_load) => :foo))
        alg = NonlinearSolve.GaussNewton()
        @test_throws ComponentInitError initialize_componentwise(nw; subalg=alg) # does not work with that solver
        initialize_componentwise(nw; subsolve_kwargs=(;verbose=Detailed())) # prints solver faield ouput

        # inject high tolerances to make the tol test fail
        @test_throws ComponentInitError initialize_componentwise(nw; verbose=true, subsolve_kwargs=Dict(VIndex(2) => (;reltol=100, abstol=100)))
        @test_throws SubtaskError initialize_componentwise(nw; verbose=true, subsolve_kwargs=Dict(VIndex(2) => (;reltol=100, abstol=100)), parallel=true)
        # initialize_componentwise(nw; subsolve_kwargs=(;verbose=Detailed()))

        # test parallel versions
        @test_throws SciMLBase.NonSolverError initialize_componentwise(nw; subalg=:foo, parallel=false)
        @test_throws SubtaskError initialize_componentwise(nw; subalg=Dict(VIndex(:swing_and_load) => :foo), parallel=true)
    end
end

@testset "InitFormula tests" begin
    @testset "InitFormula construction and basic functionality" begin
        # Test basic macro functionality
        if1 = @initformula begin
            :Vset = sqrt(:u_r^2 + :u_i^2)
        end

        @test if1.outsym == [:Vset]
        @test if1.sym == [:u_r, :u_i]
        @test NetworkDynamics.dim(if1) == 1

        # Test multiple assignments
        if2 = @initformula begin
            :Vset = sqrt(:u_r^2 + :u_i^2)
            :Pset = :u_r * :i_r + :u_i * :i_i
        end

        @test if2.outsym == [:Vset, :Pset]
        @test if2.sym == [:u_r, :u_i, :i_r, :i_i]
        @test NetworkDynamics.dim(if2) == 2

        # Test function call
        out = Dict{Symbol,Float64}()
        u_vals = [1.0, 2.0, 3.0, 4.0]  # values for u_r, u_i, i_r, i_i
        u_view = NetworkDynamics.SymbolicView(u_vals, [:u_r, :u_i, :i_r, :i_i])
        if2(out, u_view)

        @test out[:Vset] ≈ sqrt(1.0^2 + 2.0^2)
        @test out[:Pset] ≈ 1.0 * 3.0 + 2.0 * 4.0
    end

    @testset "InitFormula with complex expressions" begin
        # Test with more complex mathematical expressions
        if_complex = @initformula begin
            :magnitude = sqrt(:x^2 + :y^2)
            :angle = atan(:y, :x)
            :power = :voltage * :current * cos(:phase_diff)
        end

        @test if_complex.outsym == [:magnitude, :angle, :power]
        @test Set(if_complex.sym) == Set([:x, :y, :voltage, :current, :phase_diff])
        @test NetworkDynamics.dim(if_complex) == 3

        # Test evaluation
        out = Dict{Symbol,Float64}()
        vals = [3.0, 4.0, 230.0, 10.0, π/6]  # x, y, voltage, current, phase_diff
        u_view = NetworkDynamics.SymbolicView(vals, [:x, :y, :voltage, :current, :phase_diff])
        if_complex(out, u_view)

        @test out[:magnitude] ≈ 5.0
        @test out[:angle] ≈ atan(4.0, 3.0)
        @test out[:power] ≈ 230.0 * 10.0 * cos(π/6)
    end

    @testset "InitFormula validation" begin
        # Use ComponentLibrary model for testing
        swing_model = Lib.swing_mtk()

        # Test valid formula - overriding existing state symbol
        valid_state_formula = @initformula :θ = :ω + 1
        @test NetworkDynamics.assert_initformula_compat(swing_model, valid_state_formula) == valid_state_formula

        # Test another valid formula - overriding existing parameter
        valid_param_formula = @initformula :Pmech = 2 * :M + :D
        @test NetworkDynamics.assert_initformula_compat(swing_model, valid_param_formula) == valid_param_formula

        # Test invalid input symbols (symbol doesn't exist in component)
        invalid_input = @initformula :θ = :nonexistent_symbol + 1
        @test_throws ArgumentError NetworkDynamics.assert_initformula_compat(swing_model, invalid_input)

        # Test invalid output symbols (symbol doesn't exist in component)
        invalid_output = @initformula :nonexistent_output = :ω + 1
        @test_throws ArgumentError NetworkDynamics.assert_initformula_compat(swing_model, invalid_output)

        # Test setting observable symbols (not allowed)
        observable_conflict = @initformula :Pdamping = 1.0  # Pdamping is observable
        @test_throws ArgumentError NetworkDynamics.assert_initformula_compat(swing_model, observable_conflict)

        # Test multiple valid overrides
        multi_valid = @initformula begin
            :θ =  2     # θ override is valid
            :M = :D + 1     # M override is valid
            :ω =  0.1   # ω override is valid
        end
        @test NetworkDynamics.assert_initformula_compat(swing_model, multi_valid) == multi_valid
    end
end

@testset "GuessFormula tests" begin
    @testset "GuessFormula construction and basic functionality" begin
        # Test basic macro functionality
        gf1 = @guessformula begin
            :Vset = sqrt(:u_r^2 + :u_i^2)
        end

        @test gf1.outsym == [:Vset]
        @test gf1.sym == [:u_r, :u_i]
        @test NetworkDynamics.dim(gf1) == 1

        # Test multiple assignments
        gf2 = @guessformula begin
            :Vset = sqrt(:u_r^2 + :u_i^2)
            :Pset = :u_r * :i_r + :u_i * :i_i
        end

        @test gf2.outsym == [:Vset, :Pset]
        @test gf2.sym == [:u_r, :u_i, :i_r, :i_i]
        @test NetworkDynamics.dim(gf2) == 2

        # Test function call
        out = Dict{Symbol,Float64}()
        u_vals = [1.0, 2.0, 3.0, 4.0]  # values for u_r, u_i, i_r, i_i
        u_view = NetworkDynamics.SymbolicView(u_vals, [:u_r, :u_i, :i_r, :i_i])
        gf2(out, u_view)

        @test out[:Vset] ≈ sqrt(1.0^2 + 2.0^2)
        @test out[:Pset] ≈ 1.0 * 3.0 + 2.0 * 4.0
    end

    @testset "GuessFormula validation - KEY DIFFERENCE from InitFormula" begin
        # Use ComponentLibrary model for testing
        swing_model = Lib.swing_mtk()

        # Test valid formula - overriding existing state symbol
        valid_state_formula = @guessformula :θ = :ω + 1
        @test NetworkDynamics.assert_guessformula_compat(swing_model, valid_state_formula) == valid_state_formula

        # Test another valid formula - overriding existing parameter
        valid_param_formula = @guessformula :Pmech = 2 * :M + :D
        @test NetworkDynamics.assert_guessformula_compat(swing_model, valid_param_formula) == valid_param_formula

        # Test invalid input symbols (symbol doesn't exist in component)
        invalid_input = @guessformula :θ = :nonexistent_symbol + 1
        @test_throws ArgumentError NetworkDynamics.assert_guessformula_compat(swing_model, invalid_input)

        # Test invalid output symbols (symbol doesn't exist in component)
        invalid_output = @guessformula :nonexistent_output = :ω + 1
        @test_throws ArgumentError NetworkDynamics.assert_guessformula_compat(swing_model, invalid_output)

        # GuessFormulas cant read observables
        valid_observable_input = @guessformula :θ = :Pdamping / 10  # Pdamping is observable - allowed as input
        @test_throws ArgumentError NetworkDynamics.assert_guessformula_compat(swing_model, valid_observable_input)

        # But writing to observables is still not allowed
        invalid_observable_output = @guessformula :Pdamping = 1.0  # Pdamping is observable - not allowed as output
        @test_throws ArgumentError NetworkDynamics.assert_guessformula_compat(swing_model, invalid_observable_output)

        # Test multiple valid overrides
        multi_valid = @guessformula begin
            :θ = 2     # θ override is valid
            :M = :D + 1     # M override is valid
            :ω = 0.1   # ω override is valid
        end
        @test NetworkDynamics.assert_guessformula_compat(swing_model, multi_valid) == multi_valid
    end
end

@testset "apply_guess_formulas! tests - KEY DIFFERENCE: Layered lookup" begin
    using NetworkDynamics: apply_guess_formulas!

    @testset "Basic formula application with layered lookup" begin
        # Test that defaults take priority over guesses
        gf1 = @guessformula :voltage_mag = sqrt(:u_r^2 + :u_i^2)

        defaults = Dict(:u_r => 3.0, :u_i => 4.0)
        guesses = Dict(:u_r => 100.0, :u_i => 200.0)  # Should be ignored

        apply_guess_formulas!(guesses, defaults, [gf1]; verbose=true)

        # Formula should use defaults (3.0, 4.0), not guesses (100.0, 200.0)
        @test guesses[:voltage_mag] ≈ 5.0  # sqrt(3^2 + 4^2)
        # Defaults remain unchanged
        @test defaults[:u_r] == 3.0
        @test defaults[:u_i] == 4.0
        # Original guesses remain unchanged
        @test guesses[:u_r] == 100.0
        @test guesses[:u_i] == 200.0
    end

    @testset "Layered lookup: defaults prioritized, guesses fallback" begin
        # Test that formulas read from defaults first, then guesses
        gf1 = @guessformula :result = :from_default + :from_guess

        defaults = Dict(:from_default => 10.0)  # Only one value in defaults
        guesses = Dict(:from_guess => 5.0)      # Other value in guesses

        apply_guess_formulas!(guesses, defaults, [gf1]; verbose=true)

        @test guesses[:result] ≈ 15.0  # 10 (from defaults) + 5 (from guesses)
        # Original dicts unchanged except for new result
        @test defaults[:from_default] == 10.0
        @test guesses[:from_guess] == 5.0
    end

    @testset "GuessFormulas write to guesses, not defaults" begin
        # Test that GuessFormulas only modify guesses dict
        gf1 = @guessformula :new_value = :input * 3

        defaults = Dict(:input => 7.0)
        guesses = Dict{Symbol,Float64}()

        apply_guess_formulas!(guesses, defaults, [gf1]; verbose=true)

        # New value goes to guesses
        @test guesses[:new_value] ≈ 21.0
        # Defaults unchanged
        @test defaults == Dict(:input => 7.0)
    end

    @testset "Overwriting existing guesses" begin
        # Test that formulas can overwrite existing guesses
        gf1 = @guessformula :existing_guess = :new_value * 3

        defaults = Dict(:new_value => 7.0)
        guesses = Dict(:existing_guess => 100.0)  # Will be overwritten

        apply_guess_formulas!(guesses, defaults, [gf1]; verbose=true)

        # updating while also having a default
        defaults[:existing_guess] = 42
        apply_guess_formulas!(guesses, defaults, [gf1]; verbose=true)

        @test guesses[:existing_guess] ≈ 21.0  # 7 * 3, overwrites original 100.0
    end

    @testset "Error handling" begin
        # Test missing input symbol (not in defaults or guesses)
        gf_missing = @guessformula :output = :missing_symbol + 1
        defaults = Dict(:other => 5.0)
        guesses = Dict{Symbol,Float64}()

        @test_throws ArgumentError apply_guess_formulas!(guesses, defaults, [gf_missing]; verbose=false)

        # Test NaN input in defaults
        gf_nan = @guessformula :output = :input + 1
        defaults_nan = Dict(:input => NaN)
        guesses_nan = Dict{Symbol,Float64}()

        @test_throws ArgumentError apply_guess_formulas!(defaults_nan, guesses_nan, [gf_nan]; verbose=false)
    end

    @testset "Complex dependency chain with layered lookup" begin
        # Test dependencies reading from both defaults and guesses
        gf_root = @guessformula :shared = :from_default * :from_guess
        gf_branch1 = @guessformula :result1 = :shared + 1
        gf_branch2 = @guessformula :temp = :shared - 1
        gf_final = @guessformula :result2 = :temp * 3

        defaults = Dict(:from_default => 4.0)
        guesses = Dict(:from_guess => 2.0)

        apply_guess_formulas!(guesses, defaults, [gf_final, gf_branch1, gf_root, gf_branch2]; verbose=false)

        @test guesses[:shared] ≈ 8.0      # 4 * 2
        @test guesses[:result1] ≈ 9.0     # 8 + 1
        @test guesses[:temp] ≈ 7.0        # 8 - 1
        @test guesses[:result2] ≈ 21.0    # 7 * 3
        # Defaults unchanged
        @test defaults == Dict(:from_default => 4.0)
    end

    @testset "Circular dependency detection" begin
        # This should fail during topological sorting
        gf1 = @guessformula :a = :b + 1
        gf2 = @guessformula :b = :a + 1
        defaults = Dict(:start => 1.0)
        guesses = Dict{Symbol,Float64}()

        @test_throws ArgumentError apply_guess_formulas!(guesses, defaults, [gf1, gf2]; verbose=false)
    end

    @testset "has_guessformula, set_guessformula!, and get_guessformulas" begin
        swing_model = Lib.swing_mtk()

        # Initially no guess formula
        @test !has_guessformula(swing_model)

        # Set a guess formula
        gf1 = @guessformula :θ = :ω + 1
        set_guessformula!(swing_model, gf1)

        @test has_guessformula(swing_model)
        formulas = get_guessformulas(swing_model)
        @test length(formulas) == 1
        @test formulas[1] === gf1

        # Overwrite with multiple formulas
        gf2 = @guessformula :Pmech = 2 * :M
        gf_tuple = (gf1, gf2)
        set_guessformula!(swing_model, gf_tuple)

        formulas = get_guessformulas(swing_model)
        @test length(formulas) == 2
        @test formulas[1] === gf1
        @test formulas[2] === gf2
    end

    @testset "add_guessformula! - additive behavior and duplicate detection" begin
        swing_model = Lib.swing_mtk()

        gf1 = @guessformula :θ = :ω + 1
        gf2 = @guessformula :Pmech = 2 * :M

        # Add first formula
        @test add_guessformula!(swing_model, gf1) == true
        @test length(get_guessformulas(swing_model)) == 1

        # Add second formula
        @test add_guessformula!(swing_model, gf2) == true
        @test length(get_guessformulas(swing_model)) == 2

        # Try to add duplicate - should return false
        @test add_guessformula!(swing_model, gf1) == false
        @test length(get_guessformulas(swing_model)) == 2  # Still 2

        # Another duplicate
        @test add_guessformula!(swing_model, gf2) == false
        @test length(get_guessformulas(swing_model)) == 2
    end

    @testset "delete_guessformulas!" begin
        swing_model = Lib.swing_mtk()

        gf1 = @guessformula :θ = :ω + 1
        set_guessformula!(swing_model, gf1)
        @test has_guessformula(swing_model)

        delete_guessformulas!(swing_model)
        @test !has_guessformula(swing_model)
    end

    @testset "GuessFormula validation during set/add" begin
        swing_model = Lib.swing_mtk()

        # Valid formula - should work
        valid_gf = @guessformula :θ = :ω + 1
        set_guessformula!(swing_model, valid_gf; check=true)
        @test has_guessformula(swing_model)

        # Invalid formula - should error
        invalid_gf = @guessformula :θ = :nonexistent + 1
        @test_throws ArgumentError set_guessformula!(swing_model, invalid_gf; check=true)

        # With check=false, should work
        set_guessformula!(swing_model, invalid_gf; check=false)
        @test has_guessformula(swing_model)
    end

    @testset "GuessFormula combined with InitFormula" begin
        @named sal = Lib.SwingAndLoadDQ()
        vm = VertexModel(sal, [:i_r, :i_i], [:u_r, :u_i])
        add_initformula!(vm, @initformula begin
            :u_r = 1.0
            :u_i = 0.1
            :i_r = 0.5
            :i_i = 0.1
        end)
        add_guessformula!(vm, @guessformula(:swing₊V = sqrt(:u_r^2 + :u_i^2)))
        add_guessformula!(vm, @guessformula(:swing₊θ = atan(:u_i, :u_r)))
        add_initconstraint!(vm, @initconstraint begin
            :swing₊Pmech - :load₊Pset
            :load₊Qset
        end)

        initialize_component!(vm)

        @test get_guess(vm, :swing₊V) ≈ sqrt(1.0^2 + 0.1^2)
        @test get_guess(vm, :swing₊θ) ≈ atan(0.1, 1.0)
    end
end

@testset "Topological sorting of InitFormulas" begin
    using NetworkDynamics: topological_sort_formulas

    @testset "Independent formulas" begin
        # Test formulas with no dependencies - order should be preserved
        f1 = @initformula :a = :x + 1
        f2 = @initformula :b = :y + 2
        f3 = @initformula :c = :z * 3

        formulas = [f1, f2, f3]
        sorted = topological_sort_formulas(formulas)

        @test [f.outsym for f in sorted] == [[:c], [:b], [:a]]
    end

    @testset "Linear dependency chain" begin
        # A depends on B, B depends on C: C → B → A
        f_a = @initformula :final = :intermediate + 10    # depends on B
        f_b = @initformula :intermediate = :base * 2      # depends on C
        f_c = @initformula :base = :input + 1             # independent

        formulas = [f_a, f_b, f_c]  # wrong order
        sorted = topological_sort_formulas(formulas)

        # Should be reordered as C, B, A
        @test [f.outsym[1] for f in sorted] == [:base, :intermediate, :final]
    end

    @testset "Tree dependencies" begin
        # Multiple formulas depend on root formula
        f_root = @initformula :shared = :input * 2
        f_branch1 = @initformula :result1 = :shared + 1
        f_branch2 = @initformula :result2 = :shared - 1
        f_branch3 = @initformula :result3 = :shared * 3

        formulas = [f_branch2, f_branch1, f_root, f_branch3]  # mixed order
        sorted = topological_sort_formulas(formulas)

        # Root should come first
        @test sorted[1].outsym == [:shared]
        # Other three can be in any order after root
        rest_outputs = Set([f.outsym[1] for f in sorted[2:end]])
        @test rest_outputs == Set([:result1, :result2, :result3])
    end

    @testset "Complex dependency graph" begin
        # More complex: A→B→D, A→C→D
        f_a = @initformula :start = :input
        f_b = @initformula :path1 = :start + 1
        f_c = @initformula :path2 = :start * 2
        f_d = @initformula :end = :path1 + :path2

        formulas = [f_d, f_c, f_b, f_a]  # reverse order
        sorted = topological_sort_formulas(formulas)

        # A must come first, D must come last
        @test sorted[1].outsym == [:start]
        @test sorted[end].outsym == [:end]

        # B and C must come after A but before D
        middle_outputs = Set([f.outsym[1] for f in sorted[2:end-1]])
        @test middle_outputs == Set([:path1, :path2])
    end

    @testset "Error cases" begin
        # Test output symbol conflicts
        f1 = @initformula :conflict = :x + 1
        f2 = @initformula :conflict = :y + 2  # Same output symbol
        @test_throws ArgumentError topological_sort_formulas([f1, f2])

        # Test self-dependency - should fail at construction time
        @test_throws ArgumentError @initformula :self = :self + 1

        # Test circular dependency
        f_cycle1 = @initformula :a = :b + 1
        f_cycle2 = @initformula :b = :a + 1
        @test_throws ArgumentError topological_sort_formulas([f_cycle1, f_cycle2])

        # Test longer cycle A→B→C→A
        f_a = @initformula :a = :c + 1
        f_b = @initformula :b = :a + 1
        f_c = @initformula :c = :b + 1
        @test_throws ArgumentError topological_sort_formulas([f_a, f_b, f_c])
    end

    @testset "Edge cases" begin
        # Empty vector
        @test topological_sort_formulas(InitFormula[]) == InitFormula[]

        # Single formula
        f = @initformula :single = :input
        @test topological_sort_formulas([f]) == [f]

        # Two independent formulas
        f1 = @initformula :first = :x
        f2 = @initformula :second = :y
        sorted = topological_sort_formulas([f1, f2])
        @test length(sorted) == 2
        @test sorted[1].outsym ∈ [[:first], [:second]]
        @test sorted[2].outsym ∈ [[:first], [:second]]
        @test sorted[1].outsym != sorted[2].outsym
    end
end

@testset "apply_init_formulas! tests" begin
    using NetworkDynamics: apply_init_formulas!, topological_sort_formulas

    @testset "Basic formula application" begin
        # Test basic functionality with independent formulas
        f1 = @initformula :voltage_mag = sqrt(:u_r^2 + :u_i^2)
        f2 = @initformula :power = :u_r * :i_r + :u_i * :i_i

        defaults = Dict(:u_r => 3.0, :u_i => 4.0, :i_r => 2.0, :i_i => 1.0)

        result = apply_init_formulas!(defaults, [f1, f2]; verbose=true)

        @test result[:voltage_mag] ≈ 5.0  # sqrt(3^2 + 4^2)
        @test result[:power] ≈ 10.0       # 3*2 + 4*1
        @test result[:u_r] == 3.0         # original values preserved
        @test result[:u_i] == 4.0
        @test result[:i_r] == 2.0
        @test result[:i_i] == 1.0
    end

    @testset "Overwriting existing defaults" begin
        # Test that formulas can overwrite existing defaults
        f1 = @initformula :existing = :new_value * 3

        defaults = Dict(:existing => 100.0, :new_value => 7.0)

        # Test with verbose output
        result = apply_init_formulas!(defaults, [f1]; verbose=true)

        @test result[:existing] ≈ 21.0  # 7 * 3, overwrites original 100.0
        @test result[:new_value] == 7.0
    end

    @testset "Multiple output formula" begin
        # Test formula with multiple outputs
        f_multi = @initformula begin
            :mag = sqrt(:x^2 + :y^2)
            :angle = atan(:y, :x)
        end

        defaults = Dict(:x => 1.0, :y => 1.0)

        result = apply_init_formulas!(defaults, [f_multi]; verbose=true)

        @test result[:mag] ≈ sqrt(2.0)
        @test result[:angle] ≈ π/4
    end

    @testset "Error handling" begin
        # Test missing input symbol
        f_missing = @initformula :output = :missing_symbol + 1
        defaults = Dict(:other => 5.0)

        @test_throws ArgumentError apply_init_formulas!(defaults, [f_missing]; verbose=false)

        # Test NaN input
        f_nan = @initformula :output = :input + 1
        defaults_nan = Dict(:input => NaN)

        @test_throws ArgumentError apply_init_formulas!(defaults_nan, [f_nan]; verbose=false)

        # Test missing input
        f_missing_val = @initformula :output = :input + 1
        defaults_missing = Dict(:input => missing)

        @test_throws ArgumentError apply_init_formulas!(defaults_missing, [f_missing_val]; verbose=false)

        # Test nothing input
        f_nothing = @initformula :output = :input + 1
        defaults_nothing = Dict(:input => nothing)

        @test_throws ArgumentError apply_init_formulas!(defaults_nothing, [f_nothing]; verbose=false)
    end

    @testset "Complex dependency chain" begin
        # Test more complex dependencies: root → branch1, root → branch2 → final
        f_root = @initformula :shared = :input * 2
        f_branch1 = @initformula :result1 = :shared + 1
        f_branch2 = @initformula :temp = :shared - 1
        f_final = @initformula :result2 = :temp * 3

        defaults = Dict(:input => 4.0)

        # Mix up the order to test topological sorting
        result = apply_init_formulas!(defaults, [f_final, f_branch1, f_root, f_branch2]; verbose=false)

        @test result[:shared] ≈ 8.0      # 4 * 2
        @test result[:result1] ≈ 9.0     # 8 + 1
        @test result[:temp] ≈ 7.0        # 8 - 1
        @test result[:result2] ≈ 21.0    # 7 * 3
    end

    @testset "Empty formulas list" begin
        # Test with empty formulas list
        defaults = Dict(:x => 1.0, :y => 2.0)
        original = copy(defaults)

        result = apply_init_formulas!(defaults, []; verbose=false)

        @test result == original  # Should be unchanged
    end

    @testset "Circular dependency detection" begin
        # This should fail during topological sorting
        f1 = @initformula :a = :b + 1
        f2 = @initformula :b = :a + 1
        defaults = Dict(:start => 1.0)

        @test_throws ArgumentError apply_init_formulas!(defaults, [f1, f2]; verbose=false)
    end
end

@testset "test initialization with additional initformula and initconstraint" begin
    # Create a component model for testing both mutating and non-mutating versions
    @mtkmodel TestMergeModel begin
        @variables begin
            u_r(t)=1, [description="d-voltage", output=true]
            u_i(t)=0, [description="q-voltage", output=true]
            i_r(t)=1, [description="d-current", input=true]
            i_i(t)=0.1, [description="q-current", input=true]
            θ(t), [guess=0.0, description="angle"]
            ω(t), [guess=0.0, description="frequency"]
            V_calc(t), [description="calculated voltage magnitude"]
            P_calc(t), [description="calculated power"]
        end
        @parameters begin
            M=0.005, [description="inertia"]
            D=0.1, [description="damping"]
            P_target=0.5, [guess=0.5, description="target power"]
            V_target=1.0, [description="target voltage"]
        end
        @equations begin
            Dt(θ) ~ ω
            Dt(ω) ~ 1/M * (P_target - D*ω)
            u_r ~ cos(θ)
            u_i ~ sin(θ)
            V_calc ~ sqrt(u_r^2 + u_i^2)
            P_calc ~ u_r*i_r + u_i*i_i
        end
    end

    # Set up test data
    initial_formula = @initformula :V_target = 2.0
    initial_constraint = @initconstraint :u_i
    additional_formula = @initformula :P_target = 1
    additional_constraint = @initconstraint :θ - 0.1

    # Test 1: Non-mutating version - additional parameters work same as metadata
    vm_nonmut = VertexModel(TestMergeModel(name=:test_nonmut), [:i_r, :i_i], [:u_r, :u_i])
    set_initformula!(vm_nonmut, initial_formula)
    set_initconstraint!(vm_nonmut, initial_constraint)

    result1 = initialize_component(vm_nonmut;
        additional_initformula=additional_formula,
        additional_initconstraint=additional_constraint,
        verbose=false, tol=Inf)

    add_initformula!(vm_nonmut, additional_formula)
    add_initconstraint!(vm_nonmut, additional_constraint)
    result2 = initialize_component(vm_nonmut; verbose=false, tol=Inf)

    @test result1 == result2

    # Test 2: Both versions produce identical results with additional parameters
    vm_mut = VertexModel(TestMergeModel(name=:test_mut), [:i_r, :i_i], [:u_r, :u_i])
    vm_nonmut_compare = VertexModel(TestMergeModel(name=:test_nonmut_compare), [:i_r, :i_i], [:u_r, :u_i])

    for vm in [vm_mut, vm_nonmut_compare]
        set_initformula!(vm, initial_formula)
        set_initconstraint!(vm, initial_constraint)
    end

    result_nonmut = initialize_component(vm_nonmut_compare;
        additional_initformula=additional_formula,
        additional_initconstraint=additional_constraint,
        verbose=false, tol=Inf)

    initialize_component!(vm_mut;
        additional_initformula=additional_formula,
        additional_initconstraint=additional_constraint,
        verbose=false, tol=Inf)

    # Compare results - validates both functionality and bug fix
    @test result_nonmut[:V_target] ≈ get_initial_state(vm_mut, :V_target)
    @test result_nonmut[:P_target] ≈ get_initial_state(vm_mut, :P_target)
    @test result_nonmut[:θ] ≈ get_initial_state(vm_mut, :θ)
    @test result_nonmut[:u_i] ≈ get_initial_state(vm_mut, :u_i)
end

@testset "test duplicate detection for add_initconstraint! and add_initformula!" begin
    # Create a component model for testing
    swing_model = Lib.swing_mtk()

    @testset "add_initconstraint! duplicate detection" begin
        # Create test constraints
        c1 = @initconstraint :θ - 0.1
        c2 = @initconstraint :ω + 0.5

        # Test adding first constraint returns true
        @test add_initconstraint!(swing_model, c1)
        @test has_initconstraint(swing_model)

        # Test adding second constraint returns true
        @test add_initconstraint!(swing_model, c2)
        constraints = get_initconstraints(swing_model)
        @test length(constraints) == 2

        # Test adding duplicate constraint returns false
        @test !add_initconstraint!(swing_model, c1)
        @test !add_initconstraint!(swing_model, c2)
        constraints_after = get_initconstraints(swing_model)
        @test length(constraints_after) == 2  # Should still be 2
    end

    @testset "add_initformula! duplicate detection" begin
        # Create test formulas
        f1 = @initformula :θ = :ω + 1
        f2 = @initformula :Pmech = 2 * :M + :D

        # Test adding first formula returns true
        @test add_initformula!(swing_model, f1)
        @test has_initformula(swing_model)

        # Test adding second formula returns true
        @test add_initformula!(swing_model, f2)
        formulas = get_initformulas(swing_model)
        @test length(formulas) == 2

        # Test adding duplicate formula returns false
        @test !add_initformula!(swing_model, f1)
        @test !add_initformula!(swing_model, f2)
        formulas_after = get_initformulas(swing_model)
        @test length(formulas_after) == 2  # Should still be 2
    end
end

@testset "Test combining of formulas" begin
    using NetworkDynamics: collect_formulas
    f1 = @initformula :foo = :bar
    f2 = @initformula :baz = :qux
    tup = (f1, f2)
    arr = [f1, f2]

    # test all combinations
    @test collect_formulas(nothing, nothing) === nothing
    @test collect_formulas(nothing, f1) == [f1]
    @test collect_formulas(nothing, tup) == [f1, f2]
    @test collect_formulas(nothing, arr) == [f1, f2]
    @test collect_formulas(f1, nothing) == [f1]
    @test collect_formulas(f1, f1) == [f1, f1]
    @test collect_formulas(f1, tup) == [f1, f1, f2]
    @test collect_formulas(f1, arr) == [f1, f1, f2]
    @test collect_formulas(tup, nothing) == [f1, f2]
    @test collect_formulas(tup, f1) == [f1, f2, f1]
    @test collect_formulas(tup, tup) == [f1, f2, f1, f2]
    @test collect_formulas(tup, arr) == [f1, f2, f1, f2]
    @test collect_formulas(arr, nothing) == [f1, f2]
    @test collect_formulas(arr, f1) == [f1, f2, f1]
    @test collect_formulas(arr, tup) == [f1, f2, f1, f2]
    @test collect_formulas(arr, arr) == [f1, f2, f1, f2]
end

@testset "Test cobining of constraints" begin
    using NetworkDynamics: merge_initconstraints
    c1 = @initconstraint :x + :y
    c2 = @initconstraint :y - :z
    tup = (c1, c2)
    arr = [c1, c2]

    similar_enough(a::InitConstraint, b::InitConstraint) =
        a.sym == b.sym && a.dim == b.dim && a.prettyprint == b.prettyprint

    @test merge_initconstraints(nothing, nothing) === nothing
    @test merge_initconstraints(nothing, c1) == InitConstraint(c1)
    @test similar_enough(merge_initconstraints(nothing, arr), InitConstraint(c1, c2))
    @test similar_enough(merge_initconstraints(c1, nothing), InitConstraint(c1))
    @test similar_enough(merge_initconstraints(c1, c1), InitConstraint(c1, c1))
    @test similar_enough(merge_initconstraints(c1, tup), InitConstraint(c1, c1, c2))
    @test similar_enough(merge_initconstraints(c1, arr), InitConstraint(c1, c1, c2))
    @test similar_enough(merge_initconstraints(tup, nothing), InitConstraint(c1, c2))
    @test similar_enough(merge_initconstraints(tup, c1), InitConstraint(c1, c2, c1))
    @test similar_enough(merge_initconstraints(tup, tup), InitConstraint(c1, c2, c1, c2))
    @test similar_enough(merge_initconstraints(tup, arr), InitConstraint(c1, c2, c1, c2))
    @test similar_enough(merge_initconstraints(arr, nothing), InitConstraint(c1, c2))
    @test similar_enough(merge_initconstraints(arr, c1), InitConstraint(c1, c2, c1))
    @test similar_enough(merge_initconstraints(arr, tup), InitConstraint(c1, c2, c1, c2))
    @test similar_enough(merge_initconstraints(arr, arr), InitConstraint(c1, c2, c1, c2))
end

@testset "test error messages while initialization" begin
    # Create a simple network with Kuramoto oscillators
    g = cycle_graph(5) # 5-node cycle graph
    v1s = Lib.dqbus_slack()
    v2s = Lib.dqbus_pv(Pset=1.5, Vset=1.0)
    v3s = Lib.dqbus_pq(Pset=-1.0, Qset=-0.1)
    v4s = Lib.dqbus_pq(Pset=-1.0, Qset=-0.1)
    v5s = Lib.dqbus_timedeppq(Pfun=t->-1.0+0.1*sin(t), Qset=-0.1)
    e = Lib.dqline(X=0.1, R=0.01)
    nws = Network(g, [v1s, v2s, v3s, v4s, v5s], e)

    @test_throws NetworkInitError find_fixpoint(nws; alg=SSRootfind())

    s1 = find_fixpoint(nws; alg=SSRootfind(), t=0)
    s2 = find_fixpoint(nws; alg=OrdinaryDiffEqNonlinearSolve.NewtonRaphson(), t=0)
    @test s1 ≈ s2

    u0nan = uflat(NWState(nws))
    u0nan[1] = NaN
    @test_throws NetworkInitError find_fixpoint(nws, u0nan)

    Pfun_error(t) = t==1.0 ? error("GOTCHA!") : 0.0
    @register_symbolic Pfun_error(t)
    v5s_error = Lib.dqbus_timedeppq(Pfun=Pfun_error, Qset=-0.1)
    nws_error = Network(g, [v1s, v2s, v3s, v4s, v5s_error], e)
    find_fixpoint(nws_error) # no error for t=NaN
    @test_throws NetworkInitError find_fixpoint(nws_error; t=1.0)

    #### Dynamic errors
    v1 = Lib.dqbus_swing()
    v2 = Lib.dqbus_swing()
    nw = Network(g, [v1, v2, v3s, v4s, v5s], e; dealias=true)

    default_overrides = interface_values(s1)
    # if no free variables, errors in init_residual
    @test_throws ComponentInitError initialize_componentwise(nw; subverbose=false, verbose=false, default_overrides)
    # goes through
    initialize_componentwise(nw; subverbose=false, verbose=false, default_overrides, t=0)

    delete_default!(nw[VIndex(5)], :Qset)
    set_guess!(nw[VIndex(5)], :Qset, -0.2)
    @test_throws ComponentInitError initialize_componentwise(nw; subverbose=false, verbose=false, default_overrides)
    initialize_componentwise(nw; subverbose=false, verbose=false, default_overrides, t=0)

    # test parallel versions of error handling
    @test_throws SubtaskError initialize_componentwise(nw; subverbose=false, verbose=false, default_overrides, parallel=true)
    initialize_componentwise(nw; subverbose=false, verbose=false, default_overrides, t=0, parallel=true)

    delete_default!(nw[VIndex(5)], :Qset)
    set_guess!(nw[VIndex(5)], :Qset, -0.2)
    @test_throws SubtaskError initialize_componentwise(nw; subverbose=false, verbose=false, default_overrides, parallel=true)
    initialize_componentwise(nw; subverbose=false, verbose=false, default_overrides, t=0, parallel=true)
end
