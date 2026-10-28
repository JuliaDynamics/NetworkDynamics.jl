using ModelingToolkitBase
using ModelingToolkitBase: t_nounits as t, D_nounits as Dt
using NetworkDynamics
using NetworkDynamics: get_initformulas, get_guessformulas, has_metadata, get_metadata
using Symbolics: SymbolicUtils
using Symbolics: Symbolics
using StableRNGs: StableRNG
using Test
mtkext = Base.get_extension(NetworkDynamics, :NetworkDynamicsMTKExt)

# Field-by-field comparison of two component models, returns the names of the differing fields.
# Generated functions are compared by type, which contains the hash of their body.
_fmt(r::NetworkDynamics.ResolutionRule) = (r.outsym, r.sym, r.label, r.provenance, r.optional)
_fmt(f::NetworkDynamics.InitFormula) = (f.outsym, f.sym, f.prettyprint, f.weak, f.optional)
_fmt(f::NetworkDynamics.GuessFormula) = (f.outsym, f.sym, f.prettyprint)
_fmt(x::AbstractVector) = map(_fmt, x)
_fmt(x::Tuple) = map(_fmt, x)
_fmt(x) = x
function _mdval(k, v)
    if k in (:observed, :equations, :full_equations, :outputeqs, :full_outputeqs)
        repr.(v)
    elseif k === :odesystem
        (repr.(equations(v)), sort(repr.(collect(initial_conditions(v)))))
    else
        _fmt(v)
    end
end
function model_diff(a, b)
    d = Symbol[]
    for fn in (:name, :sym, :psym, :obssym, :insym, :outsym, :mass_matrix, :ff, :aliasmap,
               :extin, :_outsym_flat, :_obssym_all, :symmetadata)
        isequal(getfield(a, fn), getfield(b, fn)) || push!(d, fn)
    end
    for fn in (:f, :g, :obsf)
        typeof(getfield(a, fn)) == typeof(getfield(b, fn)) || push!(d, fn)
    end
    keys(a.metadata) == keys(b.metadata) || push!(d, :metadata_keys)
    for k in intersect(keys(a.metadata), keys(b.metadata))
        isequal(_mdval(k, a.metadata[k]), _mdval(k, b.metadata[k])) || push!(d, Symbol("md_", k))
    end
    d
end

# outputs and state derivatives at random but reproducible inputs
function fg_output(c)
    outs, du, u, ins, p, _t = NetworkDynamics.rand_inputs_fg(StableRNG(1), c)
    NetworkDynamics.compfg(c)(outs, du, u, ins, p, 0.0)
    (outs, du)
end

compile_counts() = (mtkext.COMPILE_CACHE.misses, mtkext.COMPILE_CACHE.hits)

# Builds `ctor()` with the cache switched off, and twice into an empty cache (miss, then hit).
function fresh_miss_hit(ctor)
    fresh = NetworkDynamics.with_mtk_model_cache(ctor, false)
    miss, hit = NetworkDynamics.with_mtk_model_cache() do
        NetworkDynamics.wipe_mtk_model_cache!()
        ctor(), ctor()
    end
    @test compile_counts() == (1, 1)
    fresh, miss, hit
end

@testset "cache is off by default" begin
    @test mtkext.MODEL_CACHE_MODE[] === :off
    @test !NetworkDynamics.mtk_model_cache_enabled()
    @test_throws ArgumentError NetworkDynamics.set_mtk_model_cache!(:yes)
end

@testset "cached build equals fresh build" begin
    function swing()
        @variables θ(t)=0 ω(t)=0 P(t) Pel(t)
        @parameters M=1 [description = "Inertia"] D=0.1 Pm=1 [guess = 0.5]
        System([Dt(θ) ~ ω, Dt(ω) ~ 1/M * (Pm - D*ω + P), Pel ~ -P], t; name=:swing)
    end
    function line()
        @variables srcθ(t) dstθ(t) P(t) Δθ(t)
        @parameters K=1
        System([Δθ ~ dstθ - srcθ, P ~ K*sin(Δθ)], t; name=:line)
    end
    function hierarchical()
        @variables x(t)=1 u(t)
        sub = System([Dt(x) ~ -x + u], t; name=:sub)
        @variables i(t) y(t)
        System([sub.u ~ i, y ~ 2*sub.x], t; systems=[sub], name=:outer)
    end
    function extsys()
        @variables x(t)=1 i(t) ext(t) y(t)
        System([Dt(x) ~ -x + i + ext, y ~ x], t; name=:extsys)
    end
    ctors = [
        () -> VertexModel(swing(), [:P], [:θ]),
        () -> EdgeModel(line(), [:srcθ], [:dstθ], AntiSymmetric([:P])),
        () -> VertexModel(hierarchical(), [:i], [:y]),
        () -> VertexModel(extsys(), [:i], [:y]; extin=[:ext => VIndex(1, :x)]),
    ]
    for ctor in ctors
        fresh, miss, hit = fresh_miss_hit(ctor)
        @test isempty(model_diff(fresh, miss))
        @test isempty(model_diff(fresh, hit))
        @test isequal(fg_output(fresh), fg_output(hit))
        # per-instance parts are not shared between instances
        @test get_metadata(hit, :odesystem) !== get_metadata(miss, :odesystem)
        @test get_metadata(hit, :equations) !== get_metadata(miss, :equations)
    end
end

@testset "instances differing only in init/guess formulas" begin
    function formulasys(k; guess=false)
        @variables x(t)=1 y(t) i(t)
        @parameters a=1 b
        sys = System([Dt(x) ~ -a*x + b + i, y ~ x], t; name=:foo)
        guess ? NetworkDynamics.set_guessf(sys, b => k * x) : NetworkDynamics.set_initf(sys, b => k * x)
    end
    for guess in (false, true)
        getf = guess ? get_guessformulas : get_initformulas
        ctor(k) = VertexModel(formulasys(k; guess), [:i], [:y])
        fresh = NetworkDynamics.with_mtk_model_cache(false) do
            [f.prettyprint for f in getf(ctor(2))], [f.prettyprint for f in getf(ctor(3))]
        end
        cached = NetworkDynamics.with_mtk_model_cache() do
            [f.prettyprint for f in getf(ctor(2))], [f.prettyprint for f in getf(ctor(3))]
        end
        @test fresh[1] != fresh[2]
        @test cached == fresh
    end
end

# callable structs as operations, one of them printing without its field
struct CacheTestGain
    k::Float64
end
(g::CacheTestGain)(x) = g.k * x
struct CacheTestQuietGain
    k::Float64
end
(g::CacheTestQuietGain)(x) = g.k * x
Base.show(io::IO, ::CacheTestQuietGain) = print(io, "gain")

@testset "operations are compared by identity" begin
    function opsys(op)
        @variables x(t)=1 i(t) y(t)
        rhs = Symbolics.wrap(SymbolicUtils.term(op, Symbolics.unwrap(i); type=Real))
        System([Dt(x) ~ -x + rhs, y ~ x], t; name=:opsys)
    end
    build(sys) = fg_output(VertexModel(sys, [:i], [:y]))

    s1, s2 = opsys(CacheTestGain(1.0)), opsys(CacheTestGain(2.0))
    fresh = NetworkDynamics.with_mtk_model_cache(() -> (build(s1), build(s2)), false)
    cached = NetworkDynamics.with_mtk_model_cache(() -> (build(s1), build(s2)))
    @test compile_counts() == (2, 0)
    @test !isequal(fresh[1], fresh[2])
    @test isequal(cached, fresh)

    # Closures and callables which print alike. RuntimeGeneratedFunctions identifies a function
    # body by its printed form, so their generated code collides even without the cache. Only
    # the cache key can be tested here.
    key(sys) = mtkext._structural_key(sys, ([:i],), ([:y],), pairs((; ff_to_constraint=true)))
    closure(k) = x -> k * x
    q1, q2 = opsys(CacheTestQuietGain(1.0)), opsys(CacheTestQuietGain(2.0))
    c1, c2 = opsys(closure(1.0)), opsys(closure(2.0))
    @test repr(equations(q1)) == repr(equations(q2))
    @test !isequal(key(q1), key(q2))
    @test isequal(key(q1), key(opsys(CacheTestQuietGain(1.0))))
    @test !isequal(key(c1), key(c2))
    fresh = NetworkDynamics.with_mtk_model_cache(() -> (build(c1), build(c2)), false)
    @test_broken !isequal(fresh[1], fresh[2])
end

@testset "defaults stay per instance, bindings are structural" begin
    function psys(adef, bdef)
        @variables x(t)=1 i(t) y(t)
        @parameters a=adef b=bdef
        System([Dt(x) ~ -b*x + a*i, y ~ x], t; name=:psys)
    end
    @parameters a
    variants = [(1.0, 2.0), (3.0, 4.0), (1.0, 2a), (1.0, 3a), (3.0, 2a)]
    fresh = NetworkDynamics.with_mtk_model_cache(false) do
        [VertexModel(psys(v...), [:i], [:y]) for v in variants]
    end
    cached = NetworkDynamics.with_mtk_model_cache() do
        [VertexModel(psys(v...), [:i], [:y]) for v in variants]
    end
    # numeric defaults share an entry, each binding needs its own
    @test compile_counts() == (3, 2)
    for (f, c) in zip(fresh, cached)
        @test isempty(model_diff(f, c))
        @test isequal(fg_output(f), fg_output(c))
    end
    @test get_metadata(cached[1], :a, :default) == 1.0
    @test get_metadata(cached[2], :a, :default) == 3.0
end

@testset "warnings show on the filling build, verbose builds bypass the cache" begin
    function conflicting()
        @variables x(t) i(t) y(t)
        @parameters a=1
        sys = System([0 ~ -x + a + i, y ~ x], t; name=:conflicting)
        sys = NetworkDynamics.set_guessf(sys, x => 2a)
        NetworkDynamics.set_guessf(sys, x => 3a)
    end
    NetworkDynamics.with_mtk_model_cache() do
        @test_logs (:warn, r"conflicting definitions") VertexModel(conflicting(), [:i], [:y])
        @test_logs VertexModel(conflicting(), [:i], [:y])
        @test compile_counts() == (1, 1)
        @test_logs (:info,) match_mode=:any VertexModel(conflicting(), [:i], [:y]; verbose=true)
        @test compile_counts() == (1, 1)
    end
end

@testset "mtkcompile builds bypass the cache" begin
    function mtksys()
        @variables x(t)=1 i(t) y(t) z(t)
        @parameters a=1
        System([Dt(x) ~ -a*x + i, z ~ 2x, y ~ z], t; name=:mtksys)
    end
    NetworkDynamics.with_mtk_model_cache() do
        NetworkDynamics.wipe_mtk_model_cache!()
        VertexModel(mtksys(), [:i], [:y]; mtkcompile=true)
        VertexModel(mtksys(), [:i], [:y]; mtkcompile=true)
        @test compile_counts() == (0, 0)
        @test isempty(mtkext.COMPILE_CACHE.entries)
    end
end

@testset "cache lifetime" begin
    function sizedsys(n)
        @variables x(t)=1 i(t) y(t)
        System([Dt(x) ~ -n*x + i, y ~ x], t; name=:sized)
    end
    oldsize = mtkext.COMPILE_CACHE.maxsize
    try
        NetworkDynamics.with_mtk_model_cache() do
            NetworkDynamics.set_mtk_model_cache!(true; maxsize=2)
            NetworkDynamics.wipe_mtk_model_cache!()
            VertexModel(sizedsys(1), [:i], [:y])
            VertexModel(sizedsys(2), [:i], [:y])
            VertexModel(sizedsys(1), [:i], [:y]) # 1 is now the most recently used
            VertexModel(sizedsys(3), [:i], [:y]) # evicts 2
            @test length(mtkext.COMPILE_CACHE.entries) == 2
            @test compile_counts() == (3, 1)
            VertexModel(sizedsys(1), [:i], [:y])
            @test compile_counts() == (3, 2)
            VertexModel(sizedsys(2), [:i], [:y])
            @test compile_counts() == (4, 2)
        end
    finally
        NetworkDynamics.set_mtk_model_cache!(mtkext.MODEL_CACHE_MODE[]; maxsize=oldsize)
    end

    # with_mtk_model_cache restores the previous setting, switching off empties the cache
    NetworkDynamics.with_mtk_model_cache(:on) do
        VertexModel(sizedsys(1), [:i], [:y])
        NetworkDynamics.with_mtk_model_cache(false) do
            @test !NetworkDynamics.mtk_model_cache_enabled()
            @test isempty(mtkext.COMPILE_CACHE.entries)
        end
        @test mtkext.MODEL_CACHE_MODE[] === :on
    end
    @test mtkext.MODEL_CACHE_MODE[] === :off

    # a failing build leaves no entry behind
    NetworkDynamics.with_mtk_model_cache() do
        NetworkDynamics.wipe_mtk_model_cache!()
        for _ in 1:2
            @test_throws Exception VertexModel(sizedsys(1), [:i], [:nonexistent])
        end
        @test isempty(mtkext.COMPILE_CACHE.entries)
    end
end

@testset "check mode" begin
    function checksys()
        @variables x(t)=1 i(t) y(t)
        @parameters a=1
        System([Dt(x) ~ -a*x + i, y ~ x], t; name=:checksys)
    end
    fresh = NetworkDynamics.with_mtk_model_cache(() -> VertexModel(checksys(), [:i], [:y]), false)
    NetworkDynamics.with_mtk_model_cache(:check) do
        NetworkDynamics.wipe_mtk_model_cache!()
        VertexModel(checksys(), [:i], [:y])
        hit = VertexModel(checksys(), [:i], [:y])
        @test compile_counts() == (1, 1)
        @test isempty(model_diff(fresh, hit))

        # a corrupted entry is caught on the next hit
        entry = only(values(mtkext.COMPILE_CACHE.entries))
        @variables z(t)
        entry.value = (; entry.value..., observed=[entry.value.observed; z ~ 0])
        @test_throws r"different build.*observed" VertexModel(checksys(), [:i], [:y])
    end
end

@testset "concurrent builds" begin
    function concsys()
        @variables x(t)=1 i(t) y(t)
        @parameters a=1
        System([Dt(x) ~ -a*x + i, y ~ x^2], t; name=:conc)
    end
    fresh = NetworkDynamics.with_mtk_model_cache(() -> VertexModel(concsys(), [:i], [:y]), false)
    built = NetworkDynamics.with_mtk_model_cache() do
        NetworkDynamics.wipe_mtk_model_cache!()
        sys = concsys()
        fetch.([Threads.@spawn VertexModel(sys, [:i], [:y]) for _ in 1:8])
    end
    @test compile_counts() == (1, 7)
    @test all(c -> isempty(model_diff(fresh, c)), built)
end
