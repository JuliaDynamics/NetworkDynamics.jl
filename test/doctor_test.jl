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
    # don't warn on correct component
    fv = (du, u, edges, p, t) -> begin
        du[1:2] .= p
        du[3] = 4
    end
    @test_logs min_level=Logging.Warn VertexFunction(;f=fv, g=1, dim=3, pdim=2)

    # don't warn on faulty broadcast (DimensionMismatch)
    f = (du, u, edges, p, t) -> begin
        du[1:2] .= edges
        du[3] = 4
    end
    @test_logs min_level=Logging.Warn VertexFunction(;f,g=1:3,dim=3,pdim=2)
    # but error if we know the in dim
    @test_logs (:warn, ) min_level=Logging.Warn VertexFunction(;f,g=1:3,dim=3,pdim=2,indim=3)
end
