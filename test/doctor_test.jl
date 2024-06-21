using NetworkDynamics

@testset "Test AccessTracker" begin
    using NetworkDynamics: AccessTracker, reads, writes, has_reads, has_writes, has_uninit_reads, oob_reads, oob_writes, has_oob, has_similars
    a = AccessTracker(rand(5))
    a[1]
    a[2] = 1
    @test reads(a) == [1]
    @test writes(a) == [2]

    a = AccessTracker(rand(5))
    a .= 1
    @test has_reads(a) == false
    @test writes(a) == 1:5
    @test has_writes(a) == true

    a = AccessTracker(rand(5));
    has_uninit_reads(a)
    a[1] = 1
    _ = a[1]
    has_uninit_reads(a)
    _ = a[2]
    has_uninit_reads(a)
    _ = a[-4]
    @test oob_reads(a) == [-4]
    a[17] = :foo
    @test oob_writes(a) == [17]
    @test has_oob(a) == true

    a = AccessTracker(rand(5))
    v = view(a, 2:3)
    v[1] = 1
    @test writes(a) == [2]

    p = AccessTracker(rand(1))
    u = AccessTracker(rand(1))
    du = AccessTracker(rand(1))
    du .= p .* u
    @test !has_similars(du)
    @test !has_similars(p)
    @test !has_similars(u)

    du .= p[1] * u
    @test !has_similars(du)
    @test !has_similars(p)
    @test has_similars(u)

    e = AccessTracker(rand(1))
    p = AccessTracker(rand(1))
    v_s = AccessTracker(rand(1))
    v_d = AccessTracker(rand(1))
    e .= p .* (v_s .- v_d) # * σ
    @test !has_similars(e)
    @test !has_similars(p)
    @test !has_similars(v_s)
    @test !has_similars(v_d)
    e .= p .* (v_s - v_d) # * σ
    @test !has_similars(e)
    @test !has_similars(p)
    @test has_similars(v_s)
    @test !has_similars(v_d)
end

@testset "chk_component" begin
    using Logging
    @test_logs min_level=Logging.Warn ODEVertex(3,2) do du, u, edges, p, t
        du[1:2] .= p
        du[3] = 4
    end
    # don't warn on faulty broadcast (DimensionMismatch)
    @test_logs min_level=Logging.Warn ODEVertex(3,2) do du, u, edges, p, t
        du[1:2] .= edges
        du[3] = 4
    end
end
