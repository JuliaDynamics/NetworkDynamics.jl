using NetworkDynamics, Graphs
using LinearAlgebra
using SteadyStateDiffEq
using Test
using NetworkDynamics: ForwardDiff

@__MODULE__()==Main ? includet(joinpath(pkgdir(NetworkDynamics), "test", "ComponentLibrary.jl")) : (const Lib = Main.Lib)

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
    λ = jacobian_eigenvals(s0)
    @test length(λ) == 3  # Should have 3 eigenvalues for 3 vertices
    @test all(real.(λ) .≤ 0)  # All eigenvalues should have non-positive real parts
    @test nw.mass_matrix == I  # Verify unconstrained system

    # Test custom eigenvalue function
    λ_custom = jacobian_eigenvals(s0; eigvalf=eigvals)
    @test λ ≈ λ_custom

    # Test is_linear_stable function
    @test is_linear_stable(s0)
    @test_throws "state s0 is not a fixpoint" is_linear_stable(s_non_fix)
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
    @test isfixpoint(s0)

    # Test eigenvalue computation for DAE system (should use reduced Jacobian)
    λ = jacobian_eigenvals(s0)
    @test length(λ) == sum(nw.mass_matrix) # reduce first

    # Test stability
    @test is_linear_stable(s0; marginally_stable=true)

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

@testset "Test _blkdiag" begin
    using NetworkDynamics: _blkdiag
    A = [1;;]
    B = nothing
    C = [1.0 2; 3 4]
    D = [5 6; 7 8]
    diag = _blkdiag(A, B, C, D)
    @test diag[1:1,1:1] == A
    @test diag[2:3,2:3] == C
    @test diag[4:5,4:5] == D

    diag2 = _blkdiag(diag, C)
    @test diag2[1:1,1:1] == A
    @test diag2[2:3,2:3] == C
    @test diag2[4:5,4:5] == D
    @test diag2[6:7,6:7] == C

    M = LinearAlgebra.I
    @test_throws MethodError _blkdiag(M, A)
    @test_throws AssertionError _blkdiag((M,(2,2)), (B, (0,0)), (C, (3,3)))  # wrong size for C
    diag3 = _blkdiag((M,(3,3)), (B, (0,0)), (C, (2,2)))
    @test diag3[1:3,1:3] == M
    @test diag3[4:5,4:5] == C

    # rectangular matrices (B: n×m, not square)
    Brect = [1.0 2; 3 4; 5 6]  # 3×2
    diag4 = _blkdiag(Brect, Brect)
    @test size(diag4) == (6, 4)
    @test diag4[1:3, 1:2] == Brect
    @test diag4[4:6, 3:4] == Brect

    # all nothing gives nothing
    @test _blkdiag(nothing, nothing) == nothing
end

@testset "_cat_syms" begin
    using NetworkDynamics: _cat_syms
    vec = [VIndex(1,:a), VIndex(2,:b)]
    sing = VIndex(3,:c)
    @test _cat_syms(vec, sing) == vcat(vec, sing)
    @test _cat_syms((vec,2), (sing, 1), (nothing, 3)) == vcat(vec, sing, nothing, nothing, nothing)

    # all nothing gives nothing
    @test _cat_syms((nothing,0), (nothing,0)) == nothing
    @test _cat_syms((nothing,0), (sing,1), (nothing,0)) == [sing]
end

@testset "LTI composition: append, series (*), parallel (+)" begin
    # two simple SISO systems
    # s1: 2-state, 1 in, 1 out
    s1 = NetworkDescriptorSystem(;
        A=[0.0 1; -1 -0.1],
        B=reshape([0.0; 1.0],2,1),
        C=reshape([1.0, 0.0],1,2),
        D=fill(0.0,1,1),
        sym=[VIndex(1,:x), VIndex(1,:v)],
        insym=VIndex(1,:u),
        outsym=VIndex(1,:y)
    )
    # s2: 1-state, 1 in, 1 out
    s2 = NetworkDescriptorSystem(;
        A=fill(-2.0,1,1),
        B=fill(1.0,1,1),
        C=fill(1.0,1,1),
        D=fill(0.0,1,1),
        sym=nothing,
        insym=VIndex(2,:u),
        outsym=VIndex(2,:y)
    )

    # --- append ---
    sa = append(s1, s2)
    @test size(sa.A) == (3,3)
    @test size(sa.B) == (3,2)
    @test size(sa.C) == (2,3)
    @test size(sa.D) == (2,2)
    # block structure: top-left A1, bottom-right A2, zero off-diagonal
    @test sa.A[1:2,1:2] == s1.A
    @test sa.A[3:3,3:3] == s2.A
    @test sa.A[1:2,3:3] == zeros(2,1)
    @test sa.A[3:3,1:2] == zeros(1,2)
    # sym: s1 has syms, s2 has nothing → padded with nothing
    @test sa.sym == [VIndex(1,:x), VIndex(1,:v), nothing]
    @test sa.insym == [VIndex(1,:u), VIndex(2,:u)]
    @test sa.outsym == [VIndex(1,:y), VIndex(2,:y)]
    # variadic: append of 3 systems equals two nested appends
    @test Matrix(append(s1, s2, s1).A) ≈ Matrix(append(append(s1,s2),s1).A)

    # --- series: s2 * s1 (output of s1 → input of s2) ---
    ss = s2 * s1
    @test size(ss.A) == (3,3)
    @test size(ss.B) == (3,1)
    @test size(ss.C) == (1,3)
    @test size(ss.D) == (1,1)
    @test ss.insym == s1.insym   # input is s1's input
    @test ss.outsym == s2.outsym # output is s2's output
    # verify transfer function at one frequency: G2(s)*G1(s)
    s_val = 2.0im
    @test ss(s_val) ≈ s2(s_val)*s1(s_val) rtol=1e-10

    # --- parallel: s1 + s2 (same 1-in/1-out, outputs summed) ---
    # need matching dims; use s2+s2
    sp = s2 + s2
    @test size(sp.A) == (2,2)
    @test size(sp.B) == (2,1)
    @test size(sp.C) == (1,2)
    @test size(sp.D) == (1,1)
    @test sp.D ≈ s2.D + s2.D
    # transfer function: G(s) = G2(s) + G2(s)
    @test sp(s_val) ≈ 2 .* s2(s_val) rtol=1e-10

    # dimension mismatch errors
    s_2out = NetworkDescriptorSystem(;
        A=fill(-1.0,1,1),
        B=fill(1.0,1,1),
        C=fill(1.0,2,1),
        D=fill(0.0,2,1),
        insym=[VIndex(1,:u)],
        outsym=[VIndex(1,:y1), VIndex(1,:y2)]
    )
    @test_throws DimensionMismatch s2 + s_2out   # mismatched outdim
    @test_throws DimensionMismatch s_2out * s_2out  # outdim=2 ≠ indim=1

    # --- MIMO: 2 inputs, 2 outputs ---
    # sm: 2-state, 2 in, 2 out  (full B, C, D matrices)
    sm = NetworkDescriptorSystem(;
        A = [-1.0 0.5; -0.5 -2.0],
        B = [1.0 0; 0 1.0],
        C = [1.0 0; 0 1.0],
        D = zeros(2,2),
        insym = [VIndex(1,:u1), VIndex(1,:u2)],
        outsym = [VIndex(1,:y1), VIndex(1,:y2)],
    )

    # append(MIMO, SISO): 2+1 states, 2+1 in, 2+1 out
    sa_ms = append(sm, s2)
    @test size(sa_ms.A) == (3,3)
    @test size(sa_ms.B) == (3,3)
    @test size(sa_ms.C) == (3,3)
    @test size(sa_ms.D) == (3,3)
    @test sa_ms.insym  == [VIndex(1,:u1), VIndex(1,:u2), VIndex(2,:u)]
    @test sa_ms.outsym == [VIndex(1,:y1), VIndex(1,:y2), VIndex(2,:y)]

    # insym/outsym distinction: SISO stores scalars, MIMO-with-1-channel stores vectors
    s_mimo1 = NetworkDescriptorSystem(;
        A=fill(-3.0,1,1),
        B=fill(1.0,1,1),
        C=fill(1.0,1,1),
        D=fill(0.0,1,1),
        insym=[VIndex(3,:u)],
        outsym=[VIndex(3,:y)]
    )  # length-1 vectors
    @test s2.insym isa VIndex          # SISO: scalar
    @test s_mimo1.insym isa Vector     # MIMO-1: vector
    # append of SISO and MIMO-1 should concatenate into a 2-vector
    sa2 = append(s2, s_mimo1)
    @test sa2.insym isa Vector
    @test length(sa2.insym) == 2

    # series MIMO→MIMO: sm has 2 out, need a system with 2 in feeding into it
    # build s_pre: 1-state, 2 in, 2 out (wide B, D)
    s_pre = NetworkDescriptorSystem(;
        A = fill(-0.5,1,1),
        B = [1.0 0.5],          # 1×2
        C = [1.0; -1.0],        # 2×1
        D = [0.5 0; 0 0.5],     # 2×2
        insym = [VIndex(4,:u1), VIndex(4,:u2)],
        outsym = [VIndex(4,:y1), VIndex(4,:y2)],
    )
    ss_mm = sm * s_pre   # output of s_pre (2) feeds input of sm (2)
    @test size(ss_mm.A) == (3,3)
    @test size(ss_mm.B) == (3,2)
    @test size(ss_mm.C) == (2,3)
    @test size(ss_mm.D) == (2,2)
    @test ss_mm.insym  == s_pre.insym
    @test ss_mm.outsym == sm.outsym
    s_val = 2.0im
    @test ss_mm(s_val) ≈ sm(s_val) * s_pre(s_val) rtol=1e-10

    # parallel MIMO+MIMO: same 2-in/2-out
    sp_mm = sm + sm
    @test size(sp_mm.D) == (2,2)
    @test sp_mm(s_val) ≈ 2 .* sm(s_val) rtol=1e-10

    # --- scalar multiplication ---
    s3 = 3.0 * s2
    @test s3(s_val) ≈ 3.0 .* s2(s_val) rtol=1e-10
    @test s3.A == s2.A && s3.B == s2.B  # A, B unchanged
    @test (s2 * 3.0).C == s3.C && (s2 * 3.0).D == s3.D  # commutative

    # unary negation and subtraction
    @test (-s2)(s_val) ≈ -s2(s_val) rtol=1e-10
    @test (s2 - s2)(s_val) ≈ zero(s2(s_val)) atol=1e-14
    @test (sm - sm)(s_val) ≈ zeros(2,2) atol=1e-14

    # --- feedback: SISO negative feedback ---
    # T(s) = s1(s) * inv(I + s2(s)*s1(s))
    sf = feedback(s1, s2)
    @test size(sf.A) == (3,3)
    @test size(sf.B) == (3,1)
    @test size(sf.C) == (1,3)
    @test size(sf.D) == (1,1)
    @test sf.insym == s1.insym
    @test sf.outsym == s1.outsym
    T_expected = s1(s_val) * inv(I + s2(s_val) * s1(s_val))
    @test sf(s_val) ≈ T_expected rtol=1e-10

    # feedback: positive feedback
    sf_pos = feedback(s1, s2; pos=true)
    T_pos = s1(s_val) * inv(I - s2(s_val) * s1(s_val))
    @test sf_pos(s_val) ≈ T_pos rtol=1e-10

    # feedback: MIMO
    # use sm (2×2) in forward, sm in feedback
    sf_mm = feedback(sm, sm)
    @test size(sf_mm.A) == (4,4)
    @test size(sf_mm.B) == (4,2)
    @test size(sf_mm.C) == (2,4)
    T_mm = sm(s_val) * inv(I + sm(s_val) * sm(s_val))
    @test sf_mm(s_val) ≈ T_mm rtol=1e-10

    # feedback: dimension mismatch
    @test_throws DimensionMismatch feedback(s_2out, s2) # outdim=2 ≠ indim=1 for s2→s_2out
    @test_throws DimensionMismatch feedback(s2, s_2out) # outdim=2 ≠ indim=1

    # feedback: multiple frequencies
    for s_test in [0.1im, 1.0+1.0im, 10.0im]
        @test sf(s_test) ≈ s1(s_test) * inv(I + s2(s_test) * s1(s_test)) rtol=1e-10
        @test sf_mm(s_test) ≈ sm(s_test) * inv(I + sm(s_test) * sm(s_test)) rtol=1e-10
    end
end

@testset "Test open loop linearization" begin
    nw, s0 = Lib.powergridlike_network();
    vidx = VIndex(:swing_and_load)
    v1 = linearize_component(s0, vidx)
    @test dim(v1) == dim(nw[vidx])
    @test [s.subidx for s in sym(v1)] == NetworkDynamics.sym(nw[vidx])
    @test [s.subidx for s in insym(v1)] == NetworkDynamics.insym_flat(nw[vidx])
    @test [s.subidx for s in outsym(v1)] == NetworkDynamics.outsym_flat(nw[vidx])

    eidx = EIndex(1)
    e1 = linearize_component(s0, eidx)
    @test dim(e1) == dim(nw[eidx])
    @test [s.subidx for s in sym(e1)] == NetworkDynamics.sym(nw[eidx])
    @test [s.subidx for s in insym(e1)] == NetworkDynamics.insym_flat(nw[eidx])
    @test [s.subidx for s in outsym(e1)] == NetworkDynamics.outsym_flat(nw[eidx])

    # now we can go for the full system initialization
    Ynw, Zbus, Yinj = open_loop_linearization(s0)
    @test dim(Ynw) == 0

    # Ynw.D should be the nodal admittance matrix (with sign convention: current into node from network)
    # Construct it manually from edge parameters: Y_line = active/(R + jX)
    # In real 2x2 block form: y_block = [[G, -B], [B, G]] where G + jB = active/(R + jX)
    # Nodal admittance: Y_nodal[s,s] += y_block, Y_nodal[d,d] += y_block,
    #                   Y_nodal[s,d] -= y_block, Y_nodal[d,s] -= y_block
    # Ynw.D = -Y_nodal (current into node = negative of standard injection convention)
    edgelist = collect(edges(nw.im.g))
    n = nv(nw)
    Y_nodal = zeros(2n, 2n)
    for (k, e) in enumerate(edgelist)
        s, d = src(e), dst(e)
        R = s0.p[EIndex(k, :R)]
        X = s0.p[EIndex(k, :X)]
        denom = R^2 + X^2
        G = R / denom
        B = (-X) / denom
        yblock = [G -B; B G]
        si, di = 2(s-1), 2(d-1)
        Y_nodal[si+1:si+2, si+1:si+2] .+= yblock
        Y_nodal[di+1:di+2, di+1:di+2] .+= yblock
        Y_nodal[si+1:si+2, di+1:di+2] .-= yblock
        Y_nodal[di+1:di+2, si+1:si+2] .-= yblock
    end
    @test -Y_nodal ≈ Matrix(Ynw.D) rtol=1e-10

    # closed-loop eigenvalues should match full nonlinear jacobian eigenvalues
    @test Yinj === nothing  # no injector nodes in this network
    cl = feedback(Zbus, Ynw; pos=true)
    λ_cl = sort(jacobian_eigenvals(cl), by=λ->(real(λ), imag(λ)))
    λ_full = sort(jacobian_eigenvals(s0), by=λ->(real(λ), imag(λ)))
    @test λ_cl ≈ λ_full rtol=1e-8
end

@testset "open loop linearization with injectors" begin
    nw, s0 = Lib.powergridlike_injector_network()
    Ynw, Zbus, Yinj = open_loop_linearization(s0)

    # manually construct the admittance matrix from R and X values of the
    # edge models and comparing them
    _inj_i = Int[]
    for eidx in 1:ne(nw)
        NetworkDynamics.is_loopback(nw.im.edgem[eidx]) || continue
        push!(_inj_i, nw.im.edgevec[eidx].src)
    end
    unique!(sort!(_inj_i))
    bus_i = setdiff(1:nv(nw), _inj_i)
    bus_rank = Dict(v => i for (i, v) in enumerate(bus_i))  # global vertex → bus index

    # Ynw.D = -Y_nodal over bus nodes only, skipping loopback edges
    n_bus = length(bus_i)
    Y_nodal = zeros(2n_bus, 2n_bus)
    for (k, e) in enumerate(collect(edges(nw.im.g)))
        NetworkDynamics.is_loopback(nw.im.edgem[k]) && continue
        s, d = src(e), dst(e)
        R = s0.p[EIndex(k, :R)]
        X = s0.p[EIndex(k, :X)]
        denom = R^2 + X^2
        G = R / denom
        B = (-X) / denom
        yblock = [G -B; B G]
        si, di = 2(bus_rank[s]-1), 2(bus_rank[d]-1)
        Y_nodal[si+1:si+2, si+1:si+2] .+= yblock
        Y_nodal[di+1:di+2, di+1:di+2] .+= yblock
        Y_nodal[si+1:si+2, di+1:di+2] .-= yblock
        Y_nodal[di+1:di+2, si+1:si+2] .-= yblock
    end
    @test -Y_nodal ≈ Matrix(Ynw.D) rtol=1e-10

    # closed-loop eigenvalues should match full nonlinear jacobian eigenvalues
    cl = feedback(Zbus, Ynw + Yinj; pos=true)
    λ_cl = sort(jacobian_eigenvals(cl), by=λ->(real(λ), imag(λ)))
    λ_full = sort(jacobian_eigenvals(s0), by=λ->(real(λ), imag(λ)))
    @test λ_cl ≈ λ_full rtol=1e-8
end
