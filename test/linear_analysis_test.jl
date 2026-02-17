using NetworkDynamics, Graphs
using LinearAlgebra
using SteadyStateDiffEq
using Test
using NetworkDynamics: ForwardDiff

@__MODULE__()==Main ? includet(joinpath(pkgdir(NetworkDynamics), "test", "ComponentLibrary.jl")) : (const Lib = Main.Lib)

@testset "Linear Stability Tests" begin

    @testset "Diffusion network tests" begin
        # Create diffusion network once and test all functions
        g = path_graph(3)
        vertex_model = Lib.diffusion_vertex()
        edge_model = Lib.diffusion_edge()
        nw = Network(g, vertex_model, edge_model)
        s0 = find_fixpoint(nw, NWParameter(nw))

        # Test isfixpoint function
        @test isfixpoint(s0)
        @test isfixpoint(s0; tol=1e-12)
        @test isfixpoint(s0; tol=1e-15)

        # Test with non-fixpoint state
        s_non_fix = NWState(s0)
        s_non_fix.v[1, :s] = 1.0
        @test !isfixpoint(s_non_fix)

        # Test jacobian_eigenvals function
        λ = jacobian_eigenvals(nw, s0)
        @test length(λ) == 3  # Should have 3 eigenvalues for 3 vertices
        @test all(real.(λ) .≤ 0)  # All eigenvalues should have non-positive real parts
        @test nw.mass_matrix == I  # Verify unconstrained system

        # Test custom eigenvalue function
        λ_custom = jacobian_eigenvals(s0; eigvalf=eigvals)
        @test λ ≈ λ_custom

        # Test is_linear_stable function
        @test is_linear_stable(s0)
        @test_throws "state s0 is not a fixpoint" is_linear_stable(nw, s_non_fix)
    end

    @testset "Kuramoto system test" begin
        # Test with a more complex system: Kuramoto oscillators
        g = path_graph(4)
        vertex_model = Lib.kuramoto_second()
        edge_model = Lib.kuramoto_edge()

        nw = Network(g, vertex_model, edge_model)

        # Set up parameters
        s = NWState(nw)
        s.p.v[:, :Pm] .= [1, -0.5, -0.5, 0]  # Power injections
        s.p.v[:, :M] .= 1.0  # Inertia
        s.p.v[:, :D] .= 0.2  # Damping
        s.p.e[:, :K] .= 2.0  # Coupling strength

        # Find fixpoint
        s0 = find_fixpoint(nw, s)

        # Test fixpoint check
        @test isfixpoint(s0)

        # Test eigenvalue computation
        λ = jacobian_eigenvals(s0)
        @test length(λ) == 8  # 4 nodes × 2 states each

        # Test stability (should be stable for this configuration)
        @test is_linear_stable(s0; marginally_stable=true)
        show_participation_factors(s0)

        # Test that eigenvalues can be complex
        @test any(isa.(λ, Complex))
    end

    @testset "Eigenvalue sensitivity" begin
        nw, s0 = Lib.powergridlike_network()
        @test nw.mass_matrix != I  # DAE system

        sens = eigenvalue_sensitivity(s0, 1)
        n_p = length(pflat(s0))
        @test sens.mode_idx == 1
        @test length(sens.sensitivities) == n_p
        @test length(sens.scaled_sensitivities) == n_p
        @test length(sens.param_syms) == n_p

        # Cross-check: verify a sensitivity against finite-difference eigenvalue shift
        # Uses eigvals() which returns the same default ordering as eigen()
        sys0 = linearize_network(s0)
        red0 = NetworkDynamics.reduce_dae(sys0)
        ev_base = eigvals(red0.A)

        for k in 1:n_p
            abs(sens.sensitivities[k]) > 1e-10 || continue
            ε = 1e-7 * max(abs(sens.param_vals[k]), 1.0)
            p_pert = copy(pflat(s0))
            pidx = NetworkDynamics.SII.parameter_index(nw, sens.param_syms[k])
            p_pert[pidx] += ε
            h!(dx, x) = nw(dx, x, p_pert, s0.t)
            A_pert = ForwardDiff.jacobian(h!, similar(uflat(s0)), uflat(s0))
            ev_pert = eigvals(NetworkDynamics.reduce_dae(NetworkDescriptorSystem(; M=sys0.M, A=A_pert)).A)
            fd_sensitivity = (ev_pert[1] - ev_base[1]) / ε
            @test sens.sensitivities[k] ≈ fd_sensitivity rtol=1e-4
        end

        # Scaled sensitivity = p * raw sensitivity (for finite params)
        for k in eachindex(sens.param_vals)
            isfinite(sens.param_vals[k]) || continue
            @test sens.scaled_sensitivities[k] ≈ sens.param_vals[k] * sens.sensitivities[k]
        end

        # show function should not error
        io = IOBuffer()
        show_eigenvalue_sensitivity(io, sens)
        output = String(take!(io))
        @test contains(output, "Eigenvalue Sensitivity")

        # params subset
        subset = sens.param_syms[1:3]
        sens_sub = eigenvalue_sensitivity(s0, 1; params=subset)
        @test length(sens_sub.sensitivities) == 3
        @test length(sens_sub.param_syms) == 3
        @test sens_sub.sensitivities ≈ sens.sensitivities[1:3]
        @test sens_sub.param_vals ≈ sens.param_vals[1:3]

        # mode_idx out of range
        @test_throws ArgumentError eigenvalue_sensitivity(s0, 0)
        @test_throws ArgumentError eigenvalue_sensitivity(s0, 1000)
    end

    @testset "DAE system with mass matrix" begin
        nw, s0 = Lib.powergridlike_network();

        # This system has a non-identity mass matrix due to algebraic constraints
        @test nw.mass_matrix != I
        @test isfixpoint(nw, s0)

        # Test eigenvalue computation for DAE system (should use reduced Jacobian)
        λ = jacobian_eigenvals(nw, s0)
        @test length(λ) == sum(nw.mass_matrix) # reduce first

        # Test stability
        @test is_linear_stable(nw, s0; marginally_stable=true)

        show_participation_factors(s0)

        # Verify this system uses the DAE branch (reduced Jacobian approach)
        c_idx = findall(LinearAlgebra.diag(nw.mass_matrix) .== 0)
        d_idx = findall(LinearAlgebra.diag(nw.mass_matrix) .== 1)
        @test length(c_idx) > 0  # Has algebraic constraints
        @test length(d_idx) > 0  # Has differential equations
    end

    @testset "input/output linearization" begin
        nw, s0 = Lib.powergridlike_network()

        # SISO: single in/single out (not wrapped in vector) → scalar transfer function value
        sys_siso = linearize_network(s0; in=VIndex(1, :u_r), out=VIndex(1, :i_r))
        @test sys_siso.insym isa VIndex
        @test sys_siso.outsym isa VIndex
        val_siso = sys_siso(10im * 2π)
        @test val_siso isa Number
        @test !isa(val_siso, AbstractMatrix)

        # MIMO: in=[single]/out=[single] (wrapped in vector) → matrix transfer function value
        sys_mimo = linearize_network(s0; in=[VIndex(1, :u_r)], out=[VIndex(1, :i_r)])
        @test sys_mimo.insym isa Vector
        @test sys_mimo.outsym isa Vector
        val_mimo = sys_mimo(10im * 2π)
        @test val_mimo isa AbstractMatrix
        @test size(val_mimo) == (1, 1)
        @test val_mimo[1, 1] ≈ val_siso

        # One representative per perturbation class; verify classification then linearize once
        all_ins = [
            VIndex(1, :u_r),     # vo_map: vertex output
            VIndex(1, :i_r),     # vi_map: vertex input
            EIndex(1, :src_i_r), # eo_map: edge output src side
            EIndex(1, :dst_i_r), # eo_map: edge output dst side
            EIndex(1, :src_u_r), # ei_map: edge input src side
            EIndex(1, :dst_u_r), # ei_map: edge input dst side
            VIndex(2, :M),       # p_map: vertex parameter
            EIndex(1, :R),       # p_map: edge parameter
        ]
        _, maps = NetworkDynamics._classify_perturbation_channels(nw, all_ins)
        @test !isempty(maps.vo_map)
        @test !isempty(maps.vi_map)
        @test !isempty(maps.eo_map)
        @test !isempty(maps.ei_map)
        @test !isempty(maps.p_map)

        sys_all = linearize_network(s0; in=all_ins, out=VIndex(2, :u_r))
        @test size(sys_all.B, 2) == length(all_ins)
        @test sys_all(10im * 2π) isa AbstractMatrix

        # reduce_dae consistency: G_unreduced(s) ≈ G_reduced(s) at multiple frequencies
        sys_full = linearize_network(s0;
            in=[VIndex(1, :u_r), VIndex(1, :u_i)],
            out=[VIndex(1, :i_r), VIndex(1, :i_i)])
        @test sys_full.M != LinearAlgebra.I  # has algebraic constraints
        sys_red = NetworkDynamics.reduce_dae(sys_full)
        @test sys_red.M == LinearAlgebra.I   # reduced to ODE
        @test NetworkDynamics.dim(sys_red) < NetworkDynamics.dim(sys_full)
        for s in [0.1im, 1.0im, 10im, 100im, -1+1im, -0.5+5im, -1-3im, 3im, 7im, 50im] .* 2π
            @test sys_full(s) ≈ sys_red(s) rtol=1e-8
        end
    end
end
