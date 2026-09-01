using NetworkDynamics
using NetworkDynamics: get_initformulas, set_default!, set_guess!, get_metadata, get_aliasmap, obssym
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

# The other way a value travels: across a residual equation. The bus below is a current source on
# a node, where the injector's current depends on the terminal voltage, so the reduction matches
# the KCL to the output and solves `j` from the physics instead:
#
#     0 ~ i + inj₊j          the KCL, a residual equation
#     inj₊j ~ K*x*cos(o)     the injector's physics, an observable of the states
#
# `i` is an input, so handing it in determines `inj₊j` through the KCL and the backward formula
# continues from there. This is what a power flow does with a bus interface current.
#
# `assume_io_coupling` is only the lever that produces this shape deterministically — it makes the
# KCL depend on the output opaquely, so matching the output to it is the cheap assignment. Without
# the flag both assignments cost the same and a model this small breaks the tie the other way.

"A current source whose output `j` lags the setpoint `pset`, modulated by the terminal `u`."
function injector(name)
    @parameters K=2.0 T=1.0 pset
    @variables u(t) j(t) x(t) [guess=0.5]
    sys = System([D(x) ~ (pset - x) / T,
                  j ~ K * x * cos(u)], t; name)
    # backward recipe: the terminal current and voltage fix the setpoint
    set_initf(sys, pset => j / (K * cos(u)); optional=true)
end

function csbus(; names=(:inj,), coeff=1)
    injs = [injector(nm) for nm in names]
    @parameters G=1.0
    @variables i(t) [input=true] o(t) [output=true, guess=0.1]
    k = coeff === :G ? G : coeff
    eqs = Equation[inj.u ~ o for inj in injs]
    push!(eqs, 0 ~ i + sum(k * inj.j for inj in injs))
    @named bus = System(eqs, t; systems=injs)
    v = VertexModel(bus, [:i], [:o]; verbose=false, assume_io_coupling=true)
    foreach(nm -> set_guess!(v, Symbol(nm, :₊pset), 0.5), names)
    v
end

# i = -0.9, o = 0.2 leave the injector at x = pset = 0.9/(K*cos(0.2))
XFIX = 0.9 / (2.0 * cos(0.2))
bus_init(v) = (io = IOBuffer();
               state = initialize_component(v; default_overrides=Dict(:i => -0.9, :o => 0.2),
                                            verbose=true, io);
               (state, String(take!(io))))

@testset "a provided input travels across a residual equation" begin
    v = csbus()
    # the voltage half of the interface is an alias class, the current half is not
    @test get_aliasmap(v) == Dict(:inj₊u => :o)
    @test :inj₊j ∈ obssym(v)
    kcl = only(filter(eq -> repr(eq.lhs) == "0", get_metadata(v, :equations)))
    @test occursin("i(t)", repr(kcl.rhs)) && occursin("inj₊j(t)", repr(kcl.rhs))

    # `inj₊j = -i` comes out of the residual, and the formula continues backwards from it
    state, log = bus_init(v)
    @test occursin("InitFormulas set:", log)
    @test !occursin("never became known", log)
    @test state[:inj₊pset] ≈ XFIX
end

@testset "a non-numeric coefficient is not inverted" begin
    # `0 ~ i + G*inj₊j` is just as invertible in principle, but the rule only accepts numeric
    # coefficients, so the formula stays skipped rather than resolving by accident
    v = csbus(; coeff=:G)
    _, log = bus_init(v)
    @test occursin("skipped, it is optional and [:inj₊j] never became known", log)
end

@testset "two injectors leave the KCL genuinely uninvertible" begin
    # `0 ~ i + inj1₊j + inj2₊j` has two unresolved symbols, so nothing can split the current
    # between them — this must stay skipped
    _, log = bus_init(csbus(; names=(:inj1, :inj2)))
    @test occursin("skipped, it is optional and [:inj1₊j] never became known", log)
    @test occursin("skipped, it is optional and [:inj2₊j] never became known", log)
end
