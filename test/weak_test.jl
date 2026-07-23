using NetworkDynamics
using NetworkDynamics: psym, get_initformulas,
                      drop_weak_formulas, apply_init_formulas!, topological_sort_formulas,
                      normalize, get_aliasmap, initialize_component!,
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
    @test occursin(", nothing, nothing, true)", rw)   # weak shown, positionally valid
    @test endswith(sprint(show, MIME"text/plain"(), rawstrong), "[:Sn], [:S_b])")  # trimmed
end

@testset "drop_weak_formulas: fire vs. yield" begin
    wf = @initformula weak=true begin :Sn = :S_b end

    # no default on the target -> kept, and it fires
    d1 = Dict{Symbol,Float64}(:S_b => 100.0)
    kept = drop_weak_formulas([wf], d1)
    @test length(kept) == 1
    apply_init_formulas!(d1, kept)
    @test d1[:Sn] == 100.0

    # target already carries a default -> dropped, yields (verbose still accounts for it)
    d2 = Dict{Symbol,Float64}(:S_b => 100.0, :Sn => 250.0)
    log = sprint() do io
        @test isempty(drop_weak_formulas([wf], d2; verbose=true, io))
    end
    @test occursin("dropping weak formula", log)
    apply_init_formulas!(d2, drop_weak_formulas([wf], d2))
    @test d2[:Sn] == 250.0

    # a strong formula is never touched by the drop, default or not
    sf = @initformula begin :Sn = :S_b end
    @test length(drop_weak_formulas([sf], d2)) == 1
end

@testset "multi-output weak with a mixed default: whole formula dropped + warn" begin
    mf = @initformula weak=true begin
        :a = :src
        :b = :src
    end
    d = Dict{Symbol,Float64}(:src => 5.0, :a => 1.0)  # only :a has a default
    kept = @test_logs (:warn, r"apply only partially") drop_weak_formulas([mf], d)
    @test isempty(kept)
    # no default on either output -> the multi-output formula participates normally
    @test length(drop_weak_formulas([mf], Dict{Symbol,Float64}(:src => 5.0))) == 1
end

@testset "weak + strong on one target: duplicate-writer error is weak-aware" begin
    strong = @initformula begin :Sn = :S_b end
    weak   = @initformula weak=true begin :Sn = :S_b end
    # both survive the pre-filter only when no default exists -> duplicate-writer error, whose
    # message must explain that weak yields to a default, not to another formula
    err = try
        topological_sort_formulas([strong, weak]); nothing
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("weak", sprint(showerror, err))
end

@testset "normalize preserves weak" begin
    vm = VertexModel(_weakdev(), [:i], [:o]; verbose=false)
    wf = @initformula weak=true begin :p = :p_src end
    nf = normalize(wf, get_aliasmap(vm), vm)
    @test nf.weak == true
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

@testset "dedupe: identical initf + initf_weak on one target keeps strong" begin
    # a target carrying both a strong `initf` and an identical weak `initf_weak` is one recipe
    # with ANDed weak flags -> strong wins (it overwrites, does not yield)
    @named sub = _weakdev(name=:sub)   # child already carries `initf_weak = p_src` on p
    @parameters q = 1.0
    @variables z(t) = 0.0
    parent = System([D(z) ~ -z + q], t; name=:par, systems=[sub])
    parent = set_initf(parent, sub.p => sub.p_src)   # strong, identical rhs
    vp = VertexModel(parent, [:sub₊i], [:sub₊o]; verbose=false)
    wf = only(filter(f -> f.outsym == [:sub₊p], collect(get_initformulas(vp))))
    @test wf.weak == false
end
