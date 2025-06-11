using NetworkDynamics, Graphs
using SteadyStateDiffEq, OrdinaryDiffEqRosenbrock
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as Dt
using SymbolicIndexingInterface: SymbolicIndexingInterface as SII
using Chairmarks


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
            V=sqrt(u_r^2 + u_i^2), [description="Voltage magnitude"]
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

    @test_throws ErrorException NetworkDynamics.initialize_component!(vf, tol=0.0; verbose=true)

    # make empty problem
    vf_comp = copy(vf)
    set_default!(vf_comp, :Pm, get_init(vf_comp, :Pm))
    set_default!(vf_comp, :θ, get_init(vf_comp, :θ))
    set_default!(vf_comp, :ω, get_init(vf_comp, :ω))

    NetworkDynamics.initialize_component!(vf_comp; verbose=true)
    dump_initial_state(vf)
    dump_initial_state(vf_comp)

    # check with wrong default/ on observed
    vf_def = copy(vf)
    set_default!(vf_def, :Pel, -100)
    set_bounds!(vf_def, :Pel, (-Inf, 0))
    @test_logs (:warn, r"has broken bounds") match_mode=:any begin
        NetworkDynamics.initialize_component!(vf_def; verbose=false)
    end
    @test_logs (:warn, r"has observables that differ") match_mode=:any begin
        NetworkDynamics.initialize_component!(vf_def; verbose=false)
    end

    vf_conflict = copy(vf)
    defaults = NetworkDynamics.get_defaults_dict(vf_conflict)
    # Modify the defaults dictionary to create conflicting values
    defaults[:u_r] = 0.0
    defaults[:u_i] = 1.0

    # Test that the metadata synchronization works with conflicting defaults
    NetworkDynamics.initialize_component!(vf_conflict;
        defaults=defaults,
        verbose=false)

    # Verify the metadata was updated and values were applied correctly
    @test get_default(vf_conflict, :u_r) == 0.0
    @test get_default(vf_conflict, :u_i) == 1.0

    # change metadtaa by providing custom input
    vf_sync = copy(vf)
    custom_defaults = NetworkDynamics.get_defaults_dict(vf_sync)
    custom_guesses = NetworkDynamics.get_guesses_dict(vf_sync)
    custom_bounds = NetworkDynamics.get_bounds_dict(vf_sync)

    delete!(custom_defaults, :u_r) # removed
    custom_defaults[:θ] = 0.0      # additional
    custom_defaults[:M] = 0.1      # changed

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
    @test_logs (:warn, r"broken bounds") match_mode=:any begin
        NetworkDynamics.initialize_component!(vf; verbose=true, apply_bound_transformation=false)
    end
    @test get_initial_state(vf, :vf) < 1 # does not conserve
    NetworkDynamics.initialize_component!(vf; verbose=true, apply_bound_transformation=true)
    @test get_initial_state(vf, :vf) > 1 # conserves

    # no lets try to force it negative
    set_bounds!(vf, :vf, (-Inf, 0))
    set_guess!(vf, :vf, -1)
    NetworkDynamics.initialize_component!(vf; verbose=true)
    @test get_initial_state(vf, :vf) < 1

    # check performance
    prob, _ = NetworkDynamics.initialization_problem(vf; apply_bound_transformation=true);
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

    state2 = initialize_component(em;
        additional_defaults=Dict(:srcθ=>0.0, :dstθ=>0.1),
        additional_guesses=Dict(:P=>0.0, :₋P=>0.0),
        verbose=true)
    @test state2[:P]  == -state[:P]
    @test state2[:₋P] == -state[:₋P]

    em_mut = copy(em)
    initialize_component!(em_mut;
        additional_defaults=Dict(:srcθ=>0.0, :dstθ=>0.1),
        additional_guesses=Dict(:P=>0.0, :₋P=>0.0),
        verbose=true)
    @test get_guess(em_mut, :P) == 0.0
    @test get_guess(em_mut, :₋P) == 0.0
    @test get_default(em_mut, :srcθ) == 0.0
    @test get_default(em_mut, :dstθ) == 0.1
end
