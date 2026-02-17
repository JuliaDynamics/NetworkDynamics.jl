using NetworkDynamics, Graphs
using LinearAlgebra
using SteadyStateDiffEq
using Test
using NetworkDynamics: ForwardDiff

(isinteractive() && @__MODULE__()==Main ? includet : include)("ComponentLibrary.jl")

@testset "Linear Stability Tests" begin

    @testset "Diffusion network tests" begin
        # Create diffusion network once and test all functions
        g = path_graph(3)
        vertex_model = Lib.diffusion_vertex()
        edge_model = Lib.diffusion_edge()
        nw = Network(g, vertex_model, edge_model)
        s0 = find_fixpoint(nw, NWParameter(nw))

        # Test isfixpoint function
        @test isfixpoint(nw, s0)
        @test isfixpoint(nw, s0; tol=1e-12)
        @test isfixpoint(nw, s0; tol=1e-15)

        # Test with non-fixpoint state
        s_non_fix = NWState(s0)
        s_non_fix.v[1, :s] = 1.0
        @test !isfixpoint(nw, s_non_fix)

        # Test jacobian_eigenvals function
        λ = jacobian_eigenvals(nw, s0)
        @test length(λ) == 3  # Should have 3 eigenvalues for 3 vertices
        @test all(real.(λ) .≤ 0)  # All eigenvalues should have non-positive real parts
        @test nw.mass_matrix == I  # Verify unconstrained system

        # Test custom eigenvalue function
        λ_custom = jacobian_eigenvals(nw, s0; eigvalf=eigvals)
        @test λ ≈ λ_custom

        # Test is_linear_stable function
        @test is_linear_stable(nw, s0)
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
        @test isfixpoint(nw, s0)

        # Test eigenvalue computation
        λ = jacobian_eigenvals(nw, s0)
        @test length(λ) == 8  # 4 nodes × 2 states each

        # Test stability (should be stable for this configuration)
        @test is_linear_stable(nw, s0; marginally_stable=true)
        show_participation_factors(s0)

        # Test that eigenvalues can be complex
        @test any(isa.(λ, Complex))
    end

    @testset "Eigenvalue sensitivity" begin
        nw, s0 = Lib.powergridlike_network()
        @test nw.mass_matrix != I  # DAE system

        sens = eigenvalue_sensitivity(nw, s0, 1)
        n_p = length(pflat(s0))
        @test sens.mode_idx == 1
        @test length(sens.sensitivities) == n_p
        @test length(sens.scaled_sensitivities) == n_p
        @test length(sens.param_syms) == n_p

        # Cross-check: verify a sensitivity against finite-difference eigenvalue shift
        # Uses eigvals() which returns the same default ordering as eigen()
        sys0 = linearize_network(nw, s0)
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

        # NWState convenience dispatch
        sens2 = eigenvalue_sensitivity(s0, 1)
        @test sens.eigenvalue ≈ sens2.eigenvalue
        @test sens.sensitivities ≈ sens2.sensitivities

        # params subset
        subset = sens.param_syms[1:3]
        sens_sub = eigenvalue_sensitivity(nw, s0, 1; params=subset)
        @test length(sens_sub.sensitivities) == 3
        @test length(sens_sub.param_syms) == 3
        @test sens_sub.sensitivities ≈ sens.sensitivities[1:3]
        @test sens_sub.param_vals ≈ sens.param_vals[1:3]

        # mode_idx out of range
        @test_throws ArgumentError eigenvalue_sensitivity(nw, s0, 0)
        @test_throws ArgumentError eigenvalue_sensitivity(nw, s0, 1000)
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
end
