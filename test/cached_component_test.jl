using ModelingToolkitBase
using ModelingToolkitBase: t_nounits as t, D_nounits as Dt
using NetworkDynamics
using NetworkDynamics: get_metadata
using Symbolics: SymbolicUtils
using StableRNGs: StableRNG
using Test
using IOCapture: IOCapture
mtkext = Base.get_extension(NetworkDynamics, :NetworkDynamicsMTKExt)

component_counts() = (mtkext.COMPONENT_CACHE.misses, mtkext.COMPONENT_CACHE.hits)
freshhash(sys) = mtkext._structural_hash(sys, zero(UInt))
# symbols of `sys` as they appear in its own equations, e.g. `machine₊M`
sym(sys, name) = mtkext.getproperty_symbolic(sys, name; might_contain_toplevel_ns=false)
defaultof(sys, name) = SymbolicUtils.unwrap_const(initial_conditions(sys)[sym(sys, name)])

# A machine with an optional lag on the mechanical power, and a bus around it, much like the
# models in PowerDynamics. `lag` is structural without being listed anywhere.
@cached_component function TestMachine(; name, lag=true)
    @variables δ(t)=0 ω(t)=0 P(t) Pmech(t)
    @parameters M=1 D=0.1 Pm=1
    eqs = [Dt(δ) ~ ω, Dt(ω) ~ (Pmech - D*ω - P)/M]
    if lag
        @parameters T=1
        @variables x(t)=0
        append!(eqs, [Dt(x) ~ (Pm - x)/T, Pmech ~ x])
    else
        push!(eqs, Pmech ~ Pm)
    end
    System(eqs, t; name)
end

@cached_component function TestBus(; name, lag=true)
    @named machine = TestMachine(; lag)
    @variables P(t) θ(t)
    System([machine.P ~ P, θ ~ machine.δ], t; systems=[machine], name)
end

"""
A lag, as in the docstring of `@cached_component`.
"""
@cached_component function TestLag(; name, T=1.0)
    @structural iszero(T)
    @variables u(t) y(t)
    if iszero(T)
        System([y ~ u], t; name)
    else
        @parameters T = T
        System([Dt(y) ~ (u - y) / T], t; name)
    end
end

# a value which `@structural` doesn't mention is part of the key
@cached_component function TestGain(; name, k=1.0)
    @variables y(t)=0 u(t)
    @parameters k = k
    System([Dt(y) ~ k*u - y], t; name)
end

# the projection is too coarse, it hides whether the lag exists
@cached_component function TestBadLag(; name, T=1.0)
    @structural T >= 0
    @variables y(t)=0 u(t)
    if iszero(T)
        System([y ~ u], t; name)
    else
        @parameters K=1
        System([Dt(y) ~ K*(u - y)], t; name)
    end
end

# derives a default from a projected kwarg, which a hit cannot replay
@cached_component function TestDerived(; name, T=1.0)
    @structural iszero(T)
    @variables y(t)=0 u(t)
    Kval = 1/T
    @parameters T=T K=Kval
    System([Dt(y) ~ K*(u - y)], t; name)
end

@cached_component function TestTwo(; name, T=1.0, lag=true)
    @structural begin
        iszero(T)
        lag
    end
    @variables y(t)=0 u(t)
    if iszero(T) || !lag
        System([y ~ u], t; name)
    else
        @parameters T = T
        System([Dt(y) ~ (u - y)/T, ], t; name)
    end
end

# `T == 0` disables the lag, `nothing` keeps it but leaves `T` without a default
@cached_component function TestOptLag(; name, T=1.0)
    @structural T == 0
    @variables y(t)=0 u(t)
    if T == 0
        System([y ~ u], t; name)
    else
        @parameters T = T
        System([Dt(y) ~ (u - y)/T], t; name)
    end
end

@cached_component function TestDoubled(; name, T=1.0)
    @structural iszero(T)
    @variables y(t)=0 u(t)
    Tval = 2T
    @parameters T = Tval
    System([Dt(y) ~ (u - y)/T], t; name)
end

# passes a projected kwarg on to a subsystem, where a hit cannot replay it
@cached_component function TestForward(; name, T=1.0)
    @structural iszero(T)
    machine = iszero(T) ? TestMachine(; name=:machine, lag=false) : TestMachine(; name=:machine, lag=true, T)
    @variables P(t)
    System([machine.P ~ P], t; systems=[machine], name)
end

@testset "hits differ only in name and defaults" begin
    NetworkDynamics.with_mtk_model_cache() do
        NetworkDynamics.wipe_mtk_model_cache!()
        a = TestMachine(; name=:a, M=2.0)
        b = TestMachine(; name=:b, D=0.5)
        c = TestMachine(; name=:c, T=0.3)
        @test component_counts() == (1, 2)
        @test nameof.((a, b, c)) == (:a, :b, :c)
        # hits share the cached equations, which match an uncached build
        @test ModelingToolkitBase.get_eqs(a) === ModelingToolkitBase.get_eqs(c)
        @test freshhash(c) == freshhash(TestMachine(; name=:x, T=0.3, cache=:off))

        # defaults of one instance don't show up in another
        @test defaultof(a, :M) == 2.0
        @test defaultof(b, :M) == 1
        @test defaultof(b, :D) == 0.5
        @test defaultof(a, :D) == 0.1
        @test defaultof(c, :T) == 0.3
        @test defaultof(b, :T) == 1

        # the bus hits as well
        bus1 = TestBus(; name=:bus1)
        bus2 = TestBus(; name=:bus2, machine₊M=3.0)
        @test freshhash(bus2) == freshhash(TestBus(; name=:x, cache=:off))
        @test defaultof(bus2, :machine₊M) == 3.0
        @test defaultof(bus1, :machine₊M) == 1

        # `nothing` can't remove the default of a subsystem
        @test_throws r"Cannot remove the default of machine₊M" TestBus(; name=:bus3, machine₊M=nothing)
        @test getfield(bus2, :gui_metadata).type == GlobalRef(@__MODULE__, :TestBus)

        # different structural choices are different entries
        nolag = TestMachine(; name=:a, lag=false)
        @test freshhash(nolag) != freshhash(a)
        @test !mtkext._names_symbol(nolag, :T)
    end
end

@testset "explicit kwargs and @structural" begin
    @test occursin("A lag, as in the docstring", string(@doc TestLag))
    NetworkDynamics.with_mtk_model_cache(:check) do
        NetworkDynamics.wipe_mtk_model_cache!()
        # every value of an explicit kwarg is its own entry
        TestGain(; name=:a)
        TestGain(; name=:b)
        g = TestGain(; name=:c, k=2.0)
        @test component_counts() == (2, 1)
        @test defaultof(g, :k) == 2.0

        # a projected kwarg which names a symbol is applied on the hit, an unused one is quiet
        TestLag(; name=:a, T=0.5)
        @test defaultof(TestLag(; name=:b, T=0.2), :T) == 0.2
        @test isempty(ModelingToolkitBase.get_ps(@test_logs TestLag(; name=:c, T=0)))
        @test component_counts() == (4, 2)

        # values of different type are different choices
        @test defaultof(TestGain(; name=:d, k=1), :k) === 1
        @test component_counts() == (5, 2)

        # several expressions, also as a block
        TestTwo(; name=:a, T=0.5, lag=true)
        b = TestTwo(; name=:b, T=0.2, lag=true)
        @test defaultof(b, :T) == 0.2
        @test length(ModelingToolkitBase.get_eqs(TestTwo(; name=:c, T=0, lag=true))) == 1
        @test length(ModelingToolkitBase.get_eqs(TestTwo(; name=:d, T=0.2, lag=false))) == 1
        @test component_counts() == (8, 3)

        # `T=nothing` falls into the same structural class as `T=1` and removes the default
        TestOptLag(; name=:a, T=1.0)
        o = TestOptLag(; name=:b, T=nothing)
        @test !haskey(initial_conditions(o), sym(o, :T))
        @test defaultof(TestOptLag(; name=:c, T=2.0), :T) == 2.0

        # the argument replaces the body's own value, which the miss warns about
        o2 = @test_logs (:warn, r"Keyword argument T of cached component TestDoubled .* replaces the value 2.0") TestDoubled(; name=:a, T=1.0)
        @test defaultof(o2, :T) == 1.0
    end
end

@testset "a projection which hides a structural choice" begin
    NetworkDynamics.with_mtk_model_cache() do
        NetworkDynamics.wipe_mtk_model_cache!()
        TestBadLag(; name=:a, T=0)
        # the hit has no parameter `K`
        @test_throws ArgumentError TestBadLag(; name=:b, T=0.5, K=2.0)

        # the subsystem has a `T`, which the miss warns about
        TestForward(; name=:a, T=0.5)
        @test defaultof(TestForward(; name=:b, T=0.2), :machine₊T) == 0.5
    end
end

@testset "check mode" begin
    NetworkDynamics.with_mtk_model_cache(:check) do
        NetworkDynamics.wipe_mtk_model_cache!()
        TestMachine(; name=:a, M=2.0)
        b = TestMachine(; name=:b, D=0.5, T=0.3)
        TestBus(; name=:bus1)
        TestBus(; name=:bus2, machine₊M=3.0)
        # the bus hits the machine, and so does its fresh build in the check
        @test component_counts() == (2, 4)
        @test defaultof(b, :T) == 0.3

        TestBadLag(; name=:a, T=0.5)
        @test_throws r"TestBadLag.*differs from a fresh build in its structure" TestBadLag(; name=:b, T=0)

        TestDerived(; name=:a, T=1.0)
        TestDerived(; name=:b, T=1.0)
        @test_throws r"TestDerived.*differs from a fresh build in its initial conditions" TestDerived(; name=:c, T=2.0)

        TestForward(; name=:a, T=0.5)
        @test_throws r"TestForward.*differs from a fresh build in its initial conditions" TestForward(; name=:b, T=0.2)
    end

    # without the check, the derived default stays stale
    NetworkDynamics.with_mtk_model_cache(:on) do
        NetworkDynamics.wipe_mtk_model_cache!()
        TestDerived(; name=:a, T=1.0)
        @test defaultof(TestDerived(; name=:b, T=2.0), :K) == 1.0
        @test defaultof(TestDerived(; name=:b, T=2.0, cache=:off), :K) == 0.5
    end
end

# no slurp, a typed kwarg, a kwarg default depending on another one, and an early return
@cached_component function TestPlain(; name, a::Float64=1.0, b=2a)
    @variables y(t)=0
    if a > 2
        @parameters a=a
        return System([Dt(y) ~ a - y], t; name)
    end
    @parameters a=a b=b
    System([Dt(y) ~ a - b*y], t; name)
end

# the guess comes from a projected kwarg which is not a symbol
@cached_component function TestGuess(; name, g=1.0)
    @structural iszero(g)
    @variables y(t) [guess=g] u(t)
    System([0 ~ y^3 - u], t; name)
end

# the projection hides the recursion depth
@cached_component function TestChain(; name, n=2)
    @structural n > 0
    @variables y(t)=0
    systems = n > 0 ? [TestChain(; name=:inner, n=n-1)] : System[]
    System([Dt(y) ~ -y], t; systems, name)
end

@testset "constructor variants" begin
    NetworkDynamics.with_mtk_model_cache(:check) do
        NetworkDynamics.wipe_mtk_model_cache!()
        TestPlain(; name=:a)
        TestPlain(; name=:b)
        p = TestPlain(; name=:c, a=1.5)
        @test (defaultof(p, :a), defaultof(p, :b)) == (1.5, 3.0)
        @test length(ModelingToolkitBase.get_ps(TestPlain(; name=:d, a=3.0))) == 1
        @test_throws TypeError TestPlain(; name=:e, a=1)
        @test component_counts() == (3, 1)

        @test_logs TestGuess(; name=:a, g=1.0)
        @test_throws r"TestGuess.*differs from a fresh build in its.*guesses" TestGuess(; name=:b, g=2.0)

    end
    NetworkDynamics.with_mtk_model_cache() do
        @test_throws r"requested again while it was being built" TestChain(; name=:chain)
        # the failed builds leave no entries behind
        @test isempty(mtkext.COMPONENT_CACHE.entries)
    end
end

@testset "the cache kwarg and cache lifetime" begin
    NetworkDynamics.with_mtk_model_cache() do
        NetworkDynamics.wipe_mtk_model_cache!()
        TestMachine(; name=:a)
        TestMachine(; name=:b, cache=:off)
        @test component_counts() == (1, 0)
        @test length(mtkext.COMPONENT_CACHE.entries) == 1
        stats = IOCapture.capture(NetworkDynamics.mtk_cache_stats).output
        @test occursin(r"Cached Component Stats:\n  Entries:  1.*\n  Misses:   1\n  Hits:     0", stats)
        NetworkDynamics.wipe_mtk_model_cache!()
        @test isempty(mtkext.COMPONENT_CACHE.entries)
        TestMachine(; name=:a)
        NetworkDynamics.with_mtk_model_cache(false) do
            @test isempty(mtkext.COMPONENT_CACHE.entries)
            TestMachine(; name=:a)
        end
        @test component_counts() == (0, 0) # the nested block reset the counters
    end

    # the kwarg overrides the global mode for one call
    NetworkDynamics.wipe_mtk_model_cache!()
    TestMachine(; name=:a)
    @test isempty(mtkext.COMPONENT_CACHE.entries)
    TestMachine(; name=:a, cache=true)
    @test length(mtkext.COMPONENT_CACHE.entries) == 1
    TestDerived(; name=:a, T=1.0, cache=:on)
    @test_throws r"TestDerived.*differs from a fresh build" TestDerived(; name=:b, T=2.0, cache=:check)
    @test_throws r"Cache mode must be" TestMachine(; name=:a, cache=:maybe)
    NetworkDynamics.wipe_mtk_model_cache!()
end

# outputs and state derivatives at random but reproducible inputs
function fg_output(c)
    outs, du, u, ins, p, _t = NetworkDynamics.rand_inputs_fg(StableRNG(1), c)
    NetworkDynamics.compfg(c)(outs, du, u, ins, p, 0.0)
    (outs, du)
end

@testset "a cached component compiles like a plain one" begin
    build() = VertexModel(TestBus(; name=:bus, machine₊M=2.0, machine₊T=0.5), [:P], [:θ])
    plain = NetworkDynamics.with_mtk_model_cache(build, false)
    cached = NetworkDynamics.with_mtk_model_cache() do
        NetworkDynamics.wipe_mtk_model_cache!()
        TestBus(; name=:other)
        build()
    end
    @test component_counts() == (2, 1)
    for fn in (:name, :sym, :psym, :obssym, :insym, :outsym, :mass_matrix, :symmetadata)
        @test isequal(getfield(plain, fn), getfield(cached, fn))
    end
    @test get_metadata(cached, :machine₊M, :default) == 2.0
    @test isequal(fg_output(plain), fg_output(cached))
end

@testset "macro" begin
    expand(ex) = NetworkDynamics.cached_component_expr(Main, ex)
    @test_throws r"keyword arguments only" expand(:(function F(x; name) end))
    @test_throws r"needs a `name`" expand(:(function F(; a) end))
    @test_throws r"reserved" expand(:(function F(; name, cache=false) end))
    @test_throws r"return type" expand(:(function F(; name)::Int end))
    @test_throws r"first statement" expand(:(function F(; name, a)
        @structural a
        @structural a
    end))
    @test_throws r"first statement" expand(:(function F(; name, a)
        x = a
        @structural x
    end))
    @test_throws r"`name` can't be structural" expand(:(function F(; name)
        @structural length(string(name))
    end))
    @test_throws r"remove them from the signature" expand(:(function F(; name, defaults...) end))
end
