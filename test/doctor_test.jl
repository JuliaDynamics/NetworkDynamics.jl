using Test
using NetworkDynamics
using ModelingToolkitBase
using ModelingToolkitBase: t_nounits as t, D_nounits as Dt
using IOCapture

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

    # report allocating f or g
    falloc = (du, u, in, p, t) -> begin
        du[1] = length(string(u[1]))
        du[2] = u[2]
    end
    @test_logs (:warn, r"Component :alloc\n- f allocates") min_level=Logging.Warn VertexModel(;f=falloc, g=1, dim=2, name=:alloc)
    galloc = (odst, vsrc, vdst, p, t) -> begin
        odst[1] = length(string(vdst[1]))
    end
    @test_logs (:warn, r"^Component :alloc\n- g allocates") min_level=Logging.Warn EdgeModel(;g=AntiSymmetric(galloc), outdim=1, indim=1, ff=PureFeedForward(), name=:alloc)

    # other argument layouts: external inputs and p=nothing for pdim=0
    fext = (dv, v, ein, ext, p, t) -> (dv[1] = length(string(ext[1])); nothing)
    @test_logs (:warn, r"- f allocates") min_level=Logging.Warn VertexModel(f=fext, g=1, dim=1, extin=[VIndex(2,:a)])
    feext = (de, e, vsrc, vdst, ext, p, t) -> (de[1] = ext[1]; nothing)
    geext = (os, od, e, vsrc, vdst, ext, p, t) -> (os[1] = length(string(ext[1])); od[1] = 0.0; nothing)
    @test_logs (:warn, r"^Component :\w+\n- g allocates") min_level=Logging.Warn EdgeModel(f=feext, g=geext, outdim=1, dim=1, extin=[EIndex(2,:a)])
    fnop = (du, u, ein, p, t) -> (du[1] = isnothing(p) ? length(string(u[1])) : 0.0; nothing)
    @test_logs (:warn, r"- f allocates") min_level=Logging.Warn VertexModel(f=fnop, g=1, dim=1, pdim=0)
end

@testset "chk_component with Duals" begin
    using Logging
    # For Duals `hi` and `lo` are Union{Float64,Int64,Dual}. With both as bounds `clamp` has too
    # many signatures to union-split, so its result is `Any` and gets boxed. For Float64 the
    # unions stay small enough to split.
    lookup(x) = x > 1 ? 1.1 : x
    funion = (du, u, in, p, t) -> begin
        lim = lookup(u[1])
        hi = u[2] > 0 ? lim : 0
        lo = u[3] > 0 ? lim : 0
        du[1] = clamp(u[4], -lo, hi)
        du[2:4] .= 0
        nothing
    end
    # construction skips the Dual pass
    v = @test_logs min_level=Logging.Warn VertexModel(f=funion, g=1, dim=4, name=:union)
    @test_logs (:warn, r"Component :union\n- f allocates .* per call, but only with ForwardDiff Duals") chk_component(v)

    # a Float64 buffer can't take the Duals
    buf = zeros(1)
    fbuf = (du, u, in, p, t) -> (buf[1] = u[1]; du[1] = buf[1]; nothing)
    v = VertexModel(f=fbuf, g=1, dim=1)
    @test_logs (:warn, r"- f fails with ForwardDiff Duals") chk_component(v)
end

@testset "chk_network" begin
    using Graphs: path_graph
    e = EdgeModel(; g=AntiSymmetric((o, vs, vd, p, t) -> o[1] = vs[1] - vd[1]), outdim=1)
    vclean = VertexModel(f=(du, u, in, p, t) -> (du[1] = in[1] - u[1]; nothing), g=1, dim=1)
    # clean also means the network code itself doesn't allocate for Duals
    cap = IOCapture.capture() do
        chk_network(Network(path_graph(3), vclean, e))
    end
    @test contains(cap.output, "✓ The network rhs doesn't allocate")
    # threaded execution allocates by design, the check uses a sequential copy
    nw = Network(path_graph(3), vclean, e; execution=ThreadedExecution{true}(), aggregator=ThreadedAggregator(+))
    cap = IOCapture.capture() do
        chk_network(nw)
    end
    @test contains(cap.output, "✓ The network rhs doesn't allocate")

    lookup(x) = x > 1 ? 1.1 : x
    funion = (du, u, in, p, t) -> begin
        lim = lookup(u[1])
        hi = in[1] > 0 ? lim : 0
        lo = u[2] > 0 ? lim : 0
        du[1] = clamp(u[1], -lo, hi)
        du[2] = 0
        nothing
    end
    vunion = VertexModel(f=funion, g=1, dim=2, name=:union)
    nw = Network(path_graph(3), [vclean, vunion, vunion], e)
    cap = IOCapture.capture() do
        chk_network(nw)
    end
    @test contains(cap.output, r"Dual states +[1-9]\d* allocations")
    @test contains(cap.output, r":union \(vertices 2, 3\)\n +✗ f allocates .* per call, but only with ForwardDiff Duals")
    @test !contains(cap.output, ":VertexM")
end

@testset "chk_component on MTK model with Duals" begin
    # the same type instability as above, written as MTK equations
    @variables x(t)=0.5 y(t)=0.3 lim(t) hi(t) lo(t) i(t) [input=true] o(t) [output=true]
    @parameters a=1 b=1
    eqs = [lim ~ ifelse(x > 1, 1.1, x),
           hi ~ ifelse(a > 0, lim, 0),
           lo ~ ifelse(b > 0, lim, 0),
           Dt(x) ~ i - x,
           Dt(y) ~ clamp(x, -lo, hi) - y,
           o ~ y]
    @named sys = System(eqs, t)
    v = VertexModel(sys, [:i], [:o])
    @test_logs (:warn, r"- f allocates .* per call, but only with ForwardDiff Duals") chk_component(v)
end
