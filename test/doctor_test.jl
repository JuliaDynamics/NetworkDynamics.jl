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
    @test_logs min_level=Logging.Warn VertexModel(;f=fv, g=1, dim=3, pdim=2)

    # don't warn on faulty broadcast (DimensionMismatch)
    f = (du, u, edges, p, t) -> begin
        du[1:2] .= edges
        du[3] = 4
    end
    @test_logs min_level=Logging.Warn VertexModel(;f,g=1:3,dim=3,pdim=2)
    # but error if we know the in dim
    @test_logs (:warn, ) min_level=Logging.Warn VertexModel(;f,g=1:3,dim=3,pdim=2,indim=3)

    # warn on allocating f or g
    falloc = (du, u, in, p, t) -> begin
        du[1] = length(string(u[1]))
        du[2] = u[2]
    end
    @test_logs (:warn, r"allocates on every call \(f: [1-9]") min_level=Logging.Warn VertexModel(;f=falloc, g=1, dim=2, name=:alloc)
    galloc = (odst, vsrc, vdst, p, t) -> begin
        odst[1] = length(string(vdst[1]))
    end
    @test_logs (:warn, r"allocates on every call \(f: 0 bytes, g: [1-9]") min_level=Logging.Warn EdgeModel(;g=AntiSymmetric(galloc), outdim=1, indim=1, ff=PureFeedForward(), name=:alloc)

    # other argument layouts: external inputs and p=nothing for pdim=0
    fext = (dv, v, ein, ext, p, t) -> (dv[1] = length(string(ext[1])); nothing)
    @test_logs (:warn, r"f: [1-9]") min_level=Logging.Warn VertexModel(f=fext, g=1, dim=1, extin=[VIndex(2,:a)])
    feext = (de, e, vsrc, vdst, ext, p, t) -> (de[1] = ext[1]; nothing)
    geext = (os, od, e, vsrc, vdst, ext, p, t) -> (os[1] = length(string(ext[1])); od[1] = 0.0; nothing)
    @test_logs (:warn, r"f: 0 bytes, g: [1-9]") min_level=Logging.Warn EdgeModel(f=feext, g=geext, outdim=1, dim=1, extin=[EIndex(2,:a)])
    fnop = (du, u, ein, p, t) -> (du[1] = isnothing(p) ? length(string(u[1])) : 0.0; nothing)
    @test_logs (:warn, r"f: [1-9]") min_level=Logging.Warn VertexModel(f=fnop, g=1, dim=1, pdim=0)
end
