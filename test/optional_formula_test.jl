using NetworkDynamics
using NetworkDynamics: get_initformulas, set_default!, set_guess!
using ModelingToolkitBase
using ModelingToolkitBase: t_nounits as t, D_nounits as D, System, @named, @variables, @parameters
using SciCompDSL: @mtkmodel
using Test

# `optional` is the source-side escape hatch: a formula whose inputs never became known is
# skipped instead of failing the initialization. That is what lets a *library block* ship both
# directions of one relation and let the surrounding model decide which one runs.
#
# The block below is that case. A first-order lag whose steady state is `y = K*u` carries both
# directions, and the container anchors either end — or neither.

# Not an `@mtkmodel`: a variable option can only reference variables declared before it, and the
# two directions reference each other. The backward one goes on the system instead, which is the
# other spelling of the same thing.
"A first-order lag with gain. In steady state `x = u` and therefore `y = K*u`."
function Lag(; name)
    @parameters K=1.5 T=2.0
    @variables u(t)
    @variables y(t) [initf_optional = K * u]
    @variables x(t) [guess=0.0]
    eqs = [D(x) ~ (u - x) / T,
           y ~ K * x]
    sys = System(eqs, t; name)
    set_initf(sys, u => y / K; optional=true)
end

# The container wires the block between the component's input and output. The factors 2 and 3
# keep `ll₊u`/`ll₊y` from collapsing into `i`/`o` as plain aliases, so the block's own symbols
# stay visible during initialization.
@mtkmodel LagBus begin
    @components begin
        ll = Lag()
    end
    @variables begin
        i(t), [input = true]
        o(t), [output = true]
    end
    @equations begin
        ll.u ~ 2 * i
        o ~ 3 * ll.y
    end
end

function lagbus(; anchor...)
    @named sys = LagBus()
    v = VertexModel(sys, [:i], [:o]; verbose=false)
    set_guess!(v, :i, 0.5)
    set_guess!(v, :o, 0.5)
    for (s, val) in anchor
        set_default!(v, s, val)
    end
    v
end

# the operating point all three anchors describe: i = 1 → ll₊u = 2 = ll₊x → ll₊y = 3 → o = 9
FIX = (; i=1.0, o=9.0, x=2.0)
verbose_init(v) = (io = IOBuffer(); state = initialize_component(v; verbose=true, io);
                   (state, String(take!(io))))

@testset "both directions attached" begin
    fs = get_initformulas(lagbus())
    @test length(fs) == 2
    @test all(f -> f.optional, fs)
    @test all(f -> !f.weak, fs)
    # the flag round-trips into the printed recipe
    @test all(f -> occursin("@initformula optional=true", sprint(show, MIME"text/plain"(), f)), fs)
end

@testset "(a) no value to start from: neither direction fires" begin
    # the operating point is fixed by a constraint, not by a value, so the init pass has nothing
    # to run on — both formulas are skipped, which would be an error without the flag
    v = lagbus()
    add_initconstraint!(v, @initconstraint :i - 1.0)
    state, log = verbose_init(v)

    @test state[:i] ≈ FIX.i
    @test state[:o] ≈ FIX.o
    @test occursin("skipped, it is optional", log)
    @test !occursin("InitFormulas set:", log)
end

@testset "(b) anchored at the output: the backward direction fires" begin
    state, log = verbose_init(lagbus(; :o => FIX.o))

    @test state[:i] ≈ FIX.i            # travelled o → ll₊y → ll₊u → i
    @test state[Symbol("ll₊x")] ≈ FIX.x
    @test occursin("ll₊u = ll₊y / ll₊K", log)
    @test !occursin("skipped", log)
end

@testset "(c) anchored at the input: the forward direction fires" begin
    state, log = verbose_init(lagbus(; :i => FIX.i))

    @test state[:o] ≈ FIX.o            # travelled i → ll₊u → ll₊y → o
    @test state[Symbol("ll₊x")] ≈ FIX.x
    @test occursin("ll₊y = ll₊K*ll₊u", log)
    @test !occursin("skipped", log)
end

@testset "without the flag, an unanchored direction is an error" begin
    v = lagbus()
    add_initconstraint!(v, @initconstraint :i - 1.0)
    strict = [InitFormula(f.f, f.outsym, f.sym, f.prettyprint; label=f.label)
              for f in get_initformulas(v)]
    delete_initformulas!(v)
    err = try
        initialize_component(v; additional_initformula=strict, verbose=false)
    catch e; e end
    @test err isa ArgumentError
    @test occursin("could not be resolved", sprint(showerror, err))
end
