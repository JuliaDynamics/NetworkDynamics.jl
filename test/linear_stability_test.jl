using NetworkDynamics, Graphs
using LinearAlgebra
using SteadyStateDiffEq
using Test

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

        # Test linear_eigenvals function
        λ = linear_eigenvals(nw, s0)
        @test length(λ) == 3  # Should have 3 eigenvalues for 3 vertices
        @test all(real.(λ) .≤ 0)  # All eigenvalues should have non-positive real parts
        @test nw.mass_matrix == I  # Verify unconstrained system

        # Test custom eigenvalue function
        λ_custom = linear_eigenvals(nw, s0; eigvalf=eigvals)
        @test λ ≈ λ_custom

        # Test is_linear_stable function
        @test is_linear_stable(nw, s0)
        @test_throws ErrorException is_linear_stable(nw, s_non_fix)
    end

    @testset "Kuramoto system test" begin
        # Test with a more complex system: Kuramoto oscillators
        g = path_graph(4)
        vertex_model = Lib.kuramoto_second()
        edge_model = Lib.kuramoto_edge()

        nw = Network(g, vertex_model, edge_model)

        # Set up parameters
        p = NWParameter(nw)
        p.v[:, :Pm] .= [1, -0.5, -0.5, 0]  # Power injections
        p.v[:, :M] .= 1.0  # Inertia
        p.v[:, :D] .= 0.1  # Damping
        p.e[:, :K] .= 2.0  # Coupling strength

        # Find fixpoint
        s0 = find_fixpoint(nw, p)

        # Test fixpoint check
        @test isfixpoint(nw, s0)

        # Test eigenvalue computation
        λ = linear_eigenvals(nw, s0)
        @test length(λ) == 8  # 4 nodes × 2 states each

        # Test stability (should be stable for this configuration)
        @test is_linear_stable(nw, s0)

        # Test that eigenvalues can be complex
        @test any(isa.(λ, Complex))
    end

    @testset "DAE system with mass matrix" begin
        # Test DAE system with mixed differential/algebraic equations
        g = path_graph(3)
        # Create vertex with mixed mass matrix: some differential, some algebraic
        mixed_vertex = VertexModel(f=(dx, x, edges, p, t) -> (dx[1] = -x[1] + edges[1]; dx[2] = 0),
                                  g=(y, x, p, t) -> y[1] = x[1],
                                  dim=2, outdim=1,
                                  mass_matrix=LinearAlgebra.Diagonal([1, 0]))  # Mixed ODE/DAE

        simple_edge = EdgeModel(g=AntiSymmetric((e, v_s, v_d, p, t) -> e[1] = v_s[1] - v_d[1]),
                               outdim=1)

        nw = Network(g, mixed_vertex, simple_edge)

        # This system has a non-identity mass matrix due to algebraic constraints
        @test nw.mass_matrix != I
        s0 = find_fixpoint(nw)
        @test isfixpoint(nw, s0)

        # Test eigenvalue computation for DAE system (should use reduced Jacobian)
        λ = linear_eigenvals(nw, s0)
        @test length(λ) == 3  # Should have 3 eigenvalues (differential variables only)
        @test all(real.(λ) .< 0)  # Should be stable

        # Test stability
        @test is_linear_stable(nw, s0)

        # Verify this system uses the DAE branch (reduced Jacobian approach)
        c_idx = findall(LinearAlgebra.diag(nw.mass_matrix) .== 0)
        d_idx = findall(LinearAlgebra.diag(nw.mass_matrix) .== 1)
        @test length(c_idx) > 0  # Has algebraic constraints
        @test length(d_idx) > 0  # Has differential equations
    end

end
