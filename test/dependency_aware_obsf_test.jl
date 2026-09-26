using NetworkDynamics
using NetworkDynamics: DependencyAwareObsf, obsf, assignment_mask, requires_input
using ModelingToolkitBase
using ModelingToolkitBase: t_nounits as t, D_nounits as Dt
using Graphs
using Chairmarks
using Test
import SymbolicIndexingInterface as SII

@__MODULE__()==Main ? includet(joinpath(pkgdir(NetworkDynamics), "test", "ComponentLibrary.jl")) : (const Lib = Main.Lib)

@component function DAOTestVertex(; name)
    @variables begin
        x(t)=1
        i(t), [input=true]
        e(t)
        o(t), [output=true]
        a(t)
        b(t)
        c(t)
        d(t)
        w(t)
        z(t)
    end
    @parameters K=1
    eqs = [
        Dt(x) ~ -x + i,
        o ~ x^2,
        a ~ 2x,    # state only
        b ~ a + K, # through another observable
        c ~ i*x,   # reads the input
        d ~ c + a, # reads the input through c
        w ~ o + 1, # through the output
        z ~ e + 1, # reads the external input
    ]
    System(eqs, t; name)
end
@named daov = DAOTestVertex()
vm = VertexModel(daov, [:i], [:o]; extin=[:e => VIndex(2, :x)])

@testset "DependencyAwareObsf of an MTK component" begin
    o = obsf(vm)
    @test o isa DependencyAwareObsf
    @test obssym(vm) == [:a, :b, :c, :d, :w, :z]

    needs = Dict(s => requires_input(o, i) for (i, s) in enumerate(obssym(vm)))
    @test needs == Dict(:a=>false, :b=>false, :c=>true, :d=>true, :w=>false, :z=>true)

    x, i, e, K = 1.5, 0.7, 0.3, 2.0
    args = ([x], [i], [e], [K], 0.0)
    full = zeros(6)
    o(full, args...)
    @test full ≈ [2x, 2x+K, i*x, i*x+2x, x^2+1, e+1]

    # observables each observable needs, besides itself
    obsneeds = Dict(1=>Int[], 2=>[1], 3=>Int[], 4=>[1, 3], 5=>Int[], 6=>Int[])
    for idxs in [[1], [2], [3], [4], [5], [6], [2, 3], [1, 4, 6], 1:6]
        out = zeros(6)
        o(out, args...; required=idxs)
        @test out[idxs] == full[idxs]
        # everything else is NaN
        computed = union(idxs, (obsneeds[i] for i in idxs)...)
        @test findall(!isnan, out) == sort(computed)

        mout = zeros(6)
        o(mout, args...; mask=assignment_mask(o, idxs))
        @test isequal(mout, out)
    end
    @test_throws ArgumentError o(zeros(6), args...; required=[1], mask=assignment_mask(o, [1]))

    @test requires_input(o, [1, 2, 5]) == false
    @test requires_input(o, [1, 3]) == true
    # Bool masks work like indices
    boolmask = BitVector([1, 0, 1, 0, 0, 0])
    @test assignment_mask(o, boolmask) == assignment_mask(o, [1, 3])
    @test requires_input(o, boolmask) == true

    # `required=` builds its mask on every call, the other styles must not allocate. This is
    # measured behind a function barrier, `@b` does not see an unspecialized varargs call.
    allocs_plain(o, out, args) = @allocated o(out, args...)
    allocs_mask(o, out, args, mask) = @allocated o(out, args...; mask)
    out = zeros(6)
    mask = assignment_mask(o, [2, 3])
    allocs_plain(o, out, args); allocs_mask(o, out, args, mask)
    @test allocs_plain(o, out, args) == 0
    @test allocs_mask(o, out, args, mask) == 0

    # copies share the obsf, components from the same system are equal
    vm2 = VertexModel(daov, [:i], [:o]; extin=[:e => VIndex(2, :x)])
    @test obsf(copy(vm)) === o
    @test obsf(vm2) == o
    @test hash(obsf(vm2)) == hash(o)
end

@testset "input check follows ff_to_constraint" begin
    @component function DAOFFVertex(; name)
        @variables begin
            x(t)=1
            i(t), [input=true]
            o(t), [output=true]
            w(t)
        end
        System([Dt(x) ~ -x, o ~ x + i, w ~ 2o], t; name)
    end
    @named daoff = DAOFFVertex()
    # with the output promoted to a state, w only reads states
    vc = VertexModel(daoff, [:i], [:o]; ff_to_constraint=true)
    @test !requires_input(obsf(vc), findfirst(==(:w), obssym(vc)))
    vf = VertexModel(daoff, [:i], [:o]; ff_to_constraint=false)
    @test requires_input(obsf(vf), findfirst(==(:w), obssym(vf)))
end

@testset "SII.observed skips the buffer fill if possible" begin
    nw = Network(path_graph(3), vm, Lib.diffusion_edge())
    u1 = [1.0, 2.0, 3.0]
    p = pflat(NWParameter(nw))
    # same batch, but different observables per member
    freeidxs = [VIndex(1,:a), VIndex(3,:b), VIndex(2,:w)]
    inidxs = [VIndex(1,:c), VIndex(2,:d), VIndex(3,:z)]

    ref = SII.observed(nw, vcat(freeidxs, inidxs))(u1, p, 0.0)
    @test ref[1:3] == [2*1.0, 2*3.0+1, 2.0^2+1]
    @test ref[6] == 2.0 + 1 # external input is x of vertex 2

    # a buffer which was not filled is all NaN
    outbuf() = NetworkDynamics.PreallocationTools.get_tmp(nw.caches.output, Float64)
    @test SII.observed(nw, freeidxs)(u1, p, 0.0) == ref[1:3]
    @test all(isnan, outbuf())
    @test SII.observed(nw, VIndex(3,:b))(u1, p, 0.0) == ref[2]
    @test all(isnan, outbuf())

    @test SII.observed(nw, inidxs)(u1, p, 0.0) == ref[4:6]
    @test !any(isnan, outbuf())

    # the masked batches run allocation free, with and without the buffer fill
    res = zeros(3)
    for idxs in (freeidxs, inidxs)
        getter = SII.observed(nw, idxs)
        @test (@b $getter($u1, $p, 0.0, $res)).allocs == 0
    end
end

@testset "plain obsf keeps the full evaluation" begin
    vplain = VertexModel(f=(dx, x, i, p, t) -> (dx[1] = -x[1] + i[1]), g=1, dim=1, sym=[:x],
                         insym=[:i], obssym=[:q, :r],
                         obsf=(out, x, i, p, t) -> (out[1] = 2x[1]; out[2] = x[1] + i[1]))
    @test requires_input(obsf(vplain), 1)
    @test isnothing(assignment_mask(obsf(vplain), [1]))
    @test requires_input(obsf(vplain), [1, 2])

    nw = Network(path_graph(2), vplain, Lib.diffusion_edge())
    u1 = [1.0, 3.0]
    p = pflat(NWParameter(nw))
    nw(zeros(2), [10.0, -10.0], p, 0.0)
    # vertex 1 receives x2 - x1 from the diffusion edge
    @test SII.observed(nw, [VIndex(1,:q), VIndex(1,:r)])(u1, p, 0.0) == [2.0, 1.0 + 2.0]
end
