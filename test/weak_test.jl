using NetworkDynamics
using NetworkDynamics: psym, get_initformulas,
                      ResolutionRule, resolve_rules, AliasMap,
                      extend_knowns_by_formulas!, get_aliasmap, initialize_component!,
                      get_default_or_init, has_default, has_init, set_default!, set_guess!
using ModelingToolkitBase
using ModelingToolkitBase: t_nounits as t, D_nounits as D, System, @variables, @parameters, @named
using SciCompDSL: @mtkmodel
using Test

# shared per-unit models live in the ComponentLibrary (`Lib`): `pu_bus` carries an injector
# rating `Sn` that *weakly* defaults to the busbar base (`initf_weak = busbar.S_b`).
@__MODULE__() == Main ? includet(joinpath(pkgdir(NetworkDynamics), "test", "ComponentLibrary.jl")) : (const Lib = Main.Lib)

# A minimal solvable weak-init component: `p` weakly defaults to `p_src`. Both symbols are used
# in equations so simplification keeps them — an unused source would be pruned and the lowered
# formula silently dropped (the same gotcha `bound_to` has). Pass `p_default` to give the target
# a genuine user default (so the weak formula yields).
function _weakdev(; name=:wd, p_default=nothing)
    @parameters p_src = 2.0
    @variables x(t) y(t) i(t)=0.0 o(t)
    if p_default === nothing
        @parameters p [initf_weak = p_src]
    else
        @parameters p=p_default [initf_weak = p_src]
    end
    System([D(x) ~ p - x + i, D(y) ~ p_src - y, o ~ x], t; name)
end

# Same shape, but `p` carries no formula — the weak recipe is attached externally via set_initf.
function _plaindev(; name=:sub)
    @parameters p_src = 2.0 p
    @variables x(t) y(t) i(t)=0.0 o(t)
    System([D(x) ~ p - x + i, D(y) ~ p_src - y, o ~ x], t; name)
end

@testset "surface: struct field, macro, show round-trips" begin
    # both the block form and the single-line form (no begin/end) parse and carry `weak`
    weak   = @initformula weak=true :Sn = :S_b
    strong = @initformula :Sn = :S_b
    @test weak.weak == true
    @test strong.weak == false
    # the printed recipe renders `weak=true` into the macro header (not a trailing marker) so
    # it stays copy-pasteable, and re-parsing it reconstructs an identical weak formula
    wshow = sprint(show, MIME"text/plain"(), weak)
    @test occursin("@initformula weak=true", wshow)
    @test occursin(":Sn = :S_b", wshow)
    @test !occursin("[weak]", wshow)
    @test eval(Meta.parse(wshow)).weak == true
    # strong prints the plain header, no weak spelling anywhere
    sshow = sprint(show, MIME"text/plain"(), strong)
    @test !occursin("weak", sshow)
    # `weak=` is the only accepted option
    @test_throws Exception macroexpand(@__MODULE__, :(@initformula bogus=1 begin :Sn = :S_b end))

    # the no-prettyprint fallback (raw constructor) stays a valid constructor call, never a
    # ` [weak]` marker: weak=true keeps the full positional form, weak=false trims to the
    # compact 3-arg form (weak omitted at its default)
    rawweak   = InitFormula((o,u)->nothing, [:Sn], [:S_b]; weak=true)
    rawstrong = InitFormula((o,u)->nothing, [:Sn], [:S_b]; weak=false)
    rw = sprint(show, MIME"text/plain"(), rawweak)
    @test !occursin("[weak]", rw)
    @test occursin(", nothing, true)", rw)   # weak shown, positionally valid
    @test endswith(sprint(show, MIME"text/plain"(), rawstrong), "[:Sn], [:S_b])")  # trimmed
end

@testset "weak: fire vs. yield" begin
    # Weakness is no longer a pre-pass that drops the formula, it is two things the resolution
    # graph does anyway: a target that already carries a value outranks a weak writer, so the
    # rule prunes; with nothing on the target the rule is the only writer and fires.
    vm = VertexModel(_weakdev(), [:i], [:o]; verbose=false)
    wf = only(get_initformulas(vm))

    # no default on the target -> it fires
    d1 = Dict{Symbol,Float64}(:p_src => 100.0)
    extend_knowns_by_formulas!(d1, vm, [wf])
    @test d1[:p] == 100.0

    # target already carries a default -> the rule prunes and the default stands
    d2 = Dict{Symbol,Float64}(:p_src => 100.0, :p => 250.0)
    log = sprint() do io
        extend_knowns_by_formulas!(d2, vm, [wf]; verbose=true, io)
    end
    @test d2[:p] == 250.0
    @test occursin("yields", log)

    # a strong formula on the same target wins, and the weak one yields silently — no default
    # needed, `:strong_formula` simply outranks `:weak_formula` at the write
    d3 = Dict{Symbol,Float64}(:p_src => 100.0)
    extend_knowns_by_formulas!(d3, vm, [wf, @initformula :p = :p_src + 1])
    @test d3[:p] == 101.0
end

@testset "weak formulas are single-output (constructor rejects multi-output)" begin
    # weak defaulting is single-target: a multi-output weak formula is rejected at construction,
    # which is what makes the drop pin-safe (it can never strand a uniquely-pinned sibling output)
    @test_throws ArgumentError (@initformula weak=true begin
        :a = :src
        :b = :src
    end)
    # a strong multi-output formula is still fine
    @test length((@initformula begin
        :a = :src
        :b = :src
    end).outsym) == 2
end

@testset "weak yields to a strong co-writer on one target" begin
    strong = @initformula begin :Sn = :S_b end
    weak   = @initformula weak=true begin :Sn = :S_b end
    # both rules fire; the weak write lands on a `:strong_formula` value and yields silently —
    # an equal-rank landing would be a check, a lower-rank one is not even that
    res = resolve_rules(Dict(:S_b => 100.0), [ResolutionRule(weak, AliasMap()),
                                              ResolutionRule(strong, AliasMap())])
    @test res.vals[:Sn] == 100.0
    @test res.provenance[:Sn] == :strong_formula
    @test isempty(res.conflicts)

    # two *weak* writers on one target used to be rejected structurally, before any value
    # existed. Now it is an ordinary equal-rank collision: agreeing is fine, and only a
    # disagreement is reported — as a conflict naming both values.
    weak2 = @initformula weak=true begin :Sn = :S_b + 1 end
    conflicting = resolve_rules(Dict(:S_b => 100.0), [ResolutionRule(weak, AliasMap()),
                                                      ResolutionRule(weak2, AliasMap())])
    c = only(conflicting.conflicts)
    # which of the two became the assignment and which the check is deliberately unspecified
    @test c.sym == :Sn
    @test Set([c.held, c.offered]) == Set([100.0, 101.0])
end

@testset "E2E component: weak fires, reinit re-fires, never freezes into a default" begin
    # The load-bearing reinit guard. `p` weakly defaults to `p_src` (a plain parameter). Init
    # once; change the source; reinit. `p` must track the source *both* times — i.e. the first
    # init did not persist the weak output as a `default` (which would freeze it and block the
    # weak formula on the second run). A dedicated solvable component: the shared `pu_bus` is
    # deliberately not at a fixpoint, so it cannot stand in for a real init here.
    vm = VertexModel(_weakdev(), [:i], [:o]; verbose=false)
    for s in (:x, :y, :o); set_guess!(vm, s, 0.5); end

    initialize_component!(vm; verbose=false)
    @test get_default_or_init(vm, :p) == 2.0    # fired: p == p_src
    @test get_default_or_init(vm, :x) == 2.0
    @test !has_default(vm, :p)                  # fresh output stored as `init`, not `default`
    @test has_init(vm, :p)

    set_default!(vm, :p_src, 5.0)
    initialize_component!(vm; verbose=false)
    @test get_default_or_init(vm, :p) == 5.0    # re-fired after the source moved
    @test get_default_or_init(vm, :x) == 5.0

    # with a genuine *user* default on the target, weak yields on both runs
    vmd = VertexModel(_weakdev(; p_default=99.0), [:i], [:o]; verbose=false)
    for s in (:x, :y, :o); set_guess!(vmd, s, 0.5); end
    initialize_component!(vmd; verbose=false)
    @test get_default_or_init(vmd, :p) == 99.0
    set_default!(vmd, :p_src, 5.0)
    initialize_component!(vmd; verbose=false)
    @test get_default_or_init(vmd, :p) == 99.0  # still yielded
end

@testset "MTK spelling: initf_weak variable option and set_initf(; weak)" begin
    # variable-level `initf_weak` lowers to a weak InitFormula, and `Sn` stays a live parameter
    # (unlike `bound_to`, which eliminates it)
    vb = Lib.pu_bus_vertex()
    @test :Sn ∈ psym(vb)
    wforms = filter(f -> f.outsym == [:Sn], collect(get_initformulas(vb)))
    @test only(wforms).weak == true

    # system-level set_initf(...; weak=true) lowers to a weak InitFormula on a subsystem target
    @named sub = _plaindev()
    @parameters q = 1.0
    @variables z(t) = 0.0
    parent = System([D(z) ~ -z + q], t; name=:par, systems=[sub])
    parent = set_initf(parent, sub.p => sub.p_src; weak=true)
    vp = VertexModel(parent, [:sub₊i], [:sub₊o]; verbose=false)
    wf = only(filter(f -> f.outsym == [:sub₊p], collect(get_initformulas(vp))))
    @test wf.weak == true
end

@testset "@mtkmodel bare-symbol initf RHS is recovered by name" begin
    # `@mtkmodel` stores a lone-identifier metadata RHS as a plain `Symbol` (not the variable), so
    # `[initf_weak = S_b]` would otherwise lower to an input-less formula that fails at init. The
    # collector recovers it by name against the model's own variables, so a bare identifier works
    # exactly like the `1*S_b` compound spelling. An unresolvable name is a hard error.
    mtkext = Base.get_extension(NetworkDynamics, :NetworkDynamicsMTKExt)

    @mtkmodel BareRecover begin
        @parameters begin
            S_b = 5.0
            Sn = 1.0, [initf_weak = S_b]     # bare identifier — mangled to `:S_b` by @mtkmodel
        end
        @variables begin
            x(t) = 0.0
        end
        @equations begin
            D(x) ~ S_b - Sn - x              # keep S_b and Sn from being pruned
        end
    end
    sys = BareRecover(; name=:c)
    e = only(filter(en -> Symbol(en.target) == :Sn, mtkext.collect_initf(sys)))
    @test e.weak == true
    @test !(e.expr isa Symbol)     # recovered to a real symbolic (not the mangled Symbol) …
    @test Symbol(e.expr) == :S_b   # … which is the S_b variable

    # a bare identifier that names no variable of the model → ArgumentError (naming the candidates)
    @mtkmodel BareUnresolvable begin
        @parameters begin
            Sn = 1.0, [initf_weak = NOPE]
        end
        @variables begin
            x(t) = 0.0
        end
        @equations begin
            D(x) ~ Sn - x
        end
    end
    @test_throws ArgumentError mtkext.collect_initf(BareUnresolvable(; name=:c))
end

@testset "metadata initf + initf_weak on one target: weak yields to strong" begin
    # a target carrying both a strong `initf` and a weak `initf_weak` keeps only the strong one
    # (weak yields to a strong writer, no conflict error) — whether the two rhs match or differ.
    strong_wins(strong_rhs) = begin
        @named sub = _weakdev(name=:sub)   # child already carries `initf_weak = p_src` on p
        @variables z(t) = 0.0
        parent = System([D(z) ~ -z], t; name=:par, systems=[sub])   # steady state z=0 (= default)
        parent = set_initf(parent, sub.p => strong_rhs(sub))
        vp = VertexModel(parent, [:sub₊i], [:sub₊o]; verbose=false)
        only(filter(f -> f.outsym == [:sub₊p], collect(get_initformulas(vp))))
    end
    @test strong_wins(sub -> sub.p_src).weak == false        # identical rhs
    @test strong_wins(sub -> sub.p_src + 10).weak == false   # differing rhs

    # and end to end: the strong value wins at init, the weak default stands down
    @named sub = _weakdev(name=:sub)
    @variables z(t) = 0.0
    parent = System([D(z) ~ -z], t; name=:par, systems=[sub])
    parent = set_initf(parent, sub.p => sub.p_src + 10)   # strong, p_src=2 -> p should be 12
    vp = VertexModel(parent, [:sub₊i], [:sub₊o]; verbose=false)
    for s in (:sub₊x, :sub₊y, :sub₊o); set_guess!(vp, s, 0.5); end
    initialize_component!(vp; verbose=false)
    @test get_default_or_init(vp, :sub₊p) == 12.0
end
