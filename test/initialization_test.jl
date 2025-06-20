using NetworkDynamics, Graphs
using SteadyStateDiffEq, OrdinaryDiffEqRosenbrock
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as Dt
using SymbolicIndexingInterface: SymbolicIndexingInterface as SII
using Chairmarks

(isinteractive() && @__MODULE__()==Main ? includet : include)("ComponentLibrary.jl")

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
        verbose=true)

    # Verify the metadata was updated and values were applied correctly
    @test get_default(vf_conflict, :u_r) == 0.0
    @test get_default(vf_conflict, :u_i) == 1.0
    @test get_init(vf_conflict, :θ) ≈ pi/2

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
    u_sview = SymbolicView(u, combined_alloc.sym)
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

@testset "test input mapping" begin
    vm = Lib.swing_mtk()
    c = @initconstraint begin
        :Pdamping # observable
        :P # input
        :θ # output
        :ω # state
        :Pmech # parameter
    end
    mapping! = NetworkDynamics.generate_init_input_mapping(vm, c)
    syms = zeros(5)
    mapping!(syms, ([3],), [NaN, 4], ([2],), [NaN, 5, NaN], [1])
    @test syms == 1:5

    em = Lib.line_mtk()
    c = @initconstraint begin
        :P # output
        :dstθ # input
        :active #param
        :Δθ # obs
    end
    mapping! = NetworkDynamics.generate_init_input_mapping(em, c)
    syms = zeros(4)
    mapping!(syms, ([NaN],[1]), [], ([NaN],[2]), [NaN, NaN, 3], [4])
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
    set_initconstraint!(vm, @initconstraint begin
        :u_mag - 1
        :Pel - 0.1
    end)
    initialize_component!(vm)
    @test get_initial_state(vm, :u_mag) ≈ 1
    @test get_initial_state(vm, :Pel) ≈ 0.1

    # test if the whole observable stuff is allocation free
    prob, _ = NetworkDynamics.initialization_problem(vm; apply_bound_transformation=true);
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
    s_meta = NWState(nw)

    for (k, v) in interface_values(pf)
        @test s_nmut[k] ≈ v atol=1e-10
    end
    @test s_meta[VIndex(1, :swing₊V)] ≈ s_meta[VIndex(1, :u_mag)]
    @test s_meta[VIndex(1, :load₊Pset)] ≈ -1.0
end
