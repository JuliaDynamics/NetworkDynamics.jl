using NetworkDynamics
using NetworkDynamics: AliasMap, canonicalize, normalize_valuedict, normalize_bounds
using ForwardDiff
using Test

# This file has two halves. Up to the "formula normalization" section the primitives under
# test are pure functions of an AliasMap and a Dict: no component, no Symbolics and
# (deliberately, cf. I10) no ModelingToolkit is involved. Everything from there on needs a
# compiled component and hence MTK.

@testset "canonicalize" begin
    am = AliasMap(:θ => (-1.0, :u_r), :s => (2.5, :x))

    @test canonicalize(am, :θ) == (-1.0, :u_r)
    @test canonicalize(am, :s) == (2.5, :x)

    # everything which is not an alias key passes through untouched (I4): canonical
    # symbols, non-alias observables and symbols unknown to any component alike
    @test canonicalize(am, :u_r) == (1.0, :u_r)
    @test canonicalize(am, :nonalias_obs) == (1.0, :nonalias_obs)
    @test canonicalize(am, :totally_unknown) == (1.0, :totally_unknown)
    @test canonicalize(AliasMap(), :θ) == (1.0, :θ)
end

@testset "normalize_valuedict" begin
    @testset "moves values onto the canonical symbol" begin
        am = AliasMap(:θ => (-1.0, :u_r), :s => (2.5, :x))
        d = Dict(:θ => 0.3, :s => 5.0, :u_i => 1.0)
        @test normalize_valuedict(am, d) == Dict(:u_r => -0.3, :x => 2.0, :u_i => 1.0)
        @test d == Dict(:θ => 0.3, :s => 5.0, :u_i => 1.0) # I3: input untouched
    end

    @testset "pass-through" begin
        d = Dict(:θ => 0.3, :nonalias_obs => 1.0)
        # empty map is the zero-work path (I1) and hands back the very same dict
        @test normalize_valuedict(AliasMap(), d) === d
        # non-alias keys survive a non-empty map, so `broken_observable_defaults` keeps
        # seeing defaults placed on genuinely algebraic observables
        @test normalize_valuedict(AliasMap(:other => (1.0, :x)), d) == d
    end

    @testset "collisions are sign-aware (I5)" begin
        am = AliasMap(:a => (-1.0, :b))

        # a = 1.0 and b = -1.0 under `a ~ -b` describe the same state, so they merge
        @test normalize_valuedict(am, Dict(:a => 1.0, :b => -1.0)) == Dict(:b => -1.0)
        # ... and would be a contradiction without the sign
        @test_throws ArgumentError normalize_valuedict(am, Dict(:a => 1.0, :b => 1.0))

        # tolerances: agreement is approximate, roundoff between the two must not throw
        @test normalize_valuedict(am, Dict(:a => 1.0, :b => -1.0 - 1e-14))[:b] ≈ -1.0

        err = try
            normalize_valuedict(am, Dict(:a => 1.0, :b => 5.0); what="default")
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        msg = sprint(showerror, err)
        @test occursin("default", msg)
        for s in [":a", ":b", "1", "5", "-1"] # both symbols, both raw and both moved values
            @test occursin(s, msg)
        end
    end

    @testset "collision resolution is deterministic" begin
        # values agree only approximately, so which one survives must not depend on hash
        # order: the canonical symbol wins over any alias ...
        am = AliasMap(:a => (-1.0, :b))
        for _ in 1:20
            @test normalize_valuedict(am, Dict(:a => -1.0, :b => 1.0 + 1e-14))[:b] == 1.0 + 1e-14
        end
        # ... and between two aliases the sorted-first one does
        am2 = AliasMap(:a => (1.0, :x), :z => (1.0, :x))
        for _ in 1:20
            @test normalize_valuedict(am2, Dict(:a => 1.0, :z => 1.0 + 1e-14))[:x] == 1.0
        end
    end

    @testset "nothing markers follow their key" begin
        am = AliasMap(:a => (-1.0, :b))
        d = Dict{Symbol,Union{Float64,Nothing}}(:a => nothing, :c => 1.0)
        @test normalize_valuedict(am, d) == Dict(:b => nothing, :c => 1.0)
        # two removal markers in one class agree, remove-vs-set does not
        @test_throws ArgumentError normalize_valuedict(am, Dict{Symbol,Any}(:a => nothing, :b => 1.0))
    end

    @testset "verbose reports the moves" begin
        am = AliasMap(:θ => (-1.0, :u_r))
        out = sprint() do io
            normalize_valuedict(am, Dict(:θ => 0.3); what="default", verbose=true, io)
        end
        @test occursin("θ", out) && occursin("u_r", out) && occursin("0.3", out)
        # silent by default
        @test sprint(io -> normalize_valuedict(am, Dict(:θ => 0.3); io)) == ""
    end
end

@testset "normalize_bounds" begin
    @testset "scaling and endpoint swap" begin
        am = AliasMap(:θ => (-1.0, :u_r), :s => (2.5, :x))
        d = Dict(:θ => (-1.0, 2.0), :s => (5.0, 10.0), :u_i => (0.0, 1.0))
        # a negative factor flips the interval, so the endpoints swap
        @test normalize_bounds(am, d) == Dict(:u_r => (-2.0, 1.0), :x => (2.0, 4.0), :u_i => (0.0, 1.0))
        @test d[:θ] == (-1.0, 2.0) # I3
    end

    @testset "collisions" begin
        am = AliasMap(:a => (-1.0, :b))
        # (-1, 2) on :a is the same interval as (-2, 1) on :b
        @test normalize_bounds(am, Dict(:a => (-1.0, 2.0), :b => (-2.0, 1.0))) == Dict(:b => (-2.0, 1.0))
        # one endpoint off is enough to conflict
        @test_throws ArgumentError normalize_bounds(am, Dict(:a => (-1.0, 2.0), :b => (-2.0, 9.0)))
    end

    @testset "pass-through" begin
        d = Dict(:θ => (0.0, 1.0))
        @test normalize_bounds(AliasMap(), d) === d
    end
end

####
#### Formula normalization — needs a compiled component, hence MTK.
####
using ModelingToolkitBase
using ModelingToolkitBase: t_nounits as t, D_nounits as Dt
using SciCompDSL
using NetworkDynamics: generate_obs_expansion, settable_symbols, obssym, normalize,
                       get_aliasmap, delete_aliasmap!, apply_init_formulas!,
                       apply_guess_formulas!, topological_sort_formulas, dim,
                       delete_metadata!
using Graphs: path_graph

# One bus carrying every shape normalization has to tell apart: a sign-flipped alias, an
# alias chain, an alias of a parameter, a genuinely algebraic (non-alias) observable, a
# multi-root one, and one that depends explicitly on time.
@mtkmodel AliasNormBus begin
    @variables begin
        u_r(t) = 1.0
        u_i(t) = 0.0
        θ(t); scaled(t); nl(t); Vmeas(t); summed(t); timed(t)
        i_r(t), [input=true]
        P(t), [output=true]
    end
    @parameters begin
        Vset = 1.0
    end
    @equations begin
        Dt(u_r) ~ -u_r + i_r
        Dt(u_i) ~ -u_i
        θ ~ -u_r                    # sign flipped alias of a state
        scaled ~ 2*θ                # chain: scaled = 2θ = -2*u_r
        nl ~ u_r^2                  # not an alias
        Vmeas ~ Vset                # alias of a parameter
        summed ~ u_r + u_i + Vset   # multi-root, not an alias
        timed ~ u_r * t             # explicitly time dependent
        P ~ u_r + nl
    end
end
@named _anb = AliasNormBus()
const VM = VertexModel(_anb, [:i_r], [:P])
const AM = get_aliasmap(VM)

@testset "generate_obs_expansion" begin
    @testset "expands to settable roots" begin
        # chain through the sign-flipped alias: scaled = 2θ = -2*u_r
        roots, f = generate_obs_expansion(VM, [:scaled])
        @test roots == [:u_r]
        @test f([3.0], NaN) ≈ [-6.0]

        roots, f = generate_obs_expansion(VM, [:summed])
        @test Set(roots) == Set([:u_r, :u_i, :Vset])
        vals = [Dict(:u_r => 2.0, :u_i => 0.5, :Vset => 1.5)[r] for r in roots]
        @test f(vals, NaN) ≈ [4.0]
    end

    @testset "symbols without an observed equation are their own root" begin
        roots, f = generate_obs_expansion(VM, [:u_r, :Vset])
        @test roots == [:u_r, :Vset]
        @test f([7.0, 9.0], NaN) ≈ [7.0, 9.0]
    end

    # this is what lets `normalize` hand over a formula's whole input list in one sweep
    @testset "one sweep over mixed symbols, roots deduplicated" begin
        roots, f = generate_obs_expansion(VM, [:θ, :u_r, :nl])
        @test roots == [:u_r]
        @test f([3.0], NaN) ≈ [-3.0, 3.0, 9.0]
    end

    @testset "time is an argument, never a root" begin
        roots, f = generate_obs_expansion(VM, [:timed])
        @test roots == [:u_r]
        @test f([2.0], 5.0) ≈ [10.0]
    end

    @testset "components without observed equations" begin
        cf = VertexModel(f=(du, u, in, p, t) -> du .= u, g=1:1, sym=[:x], outsym=[:o])
        @test_throws ArgumentError generate_obs_expansion(cf, [:x])
    end
end

@testset "normalize(::InitFormula)" begin
    @testset "identity fast path returns ===" begin
        f = @initformula :Vset = :u_r + :u_i
        @test normalize(f, AM, VM) === f
        @test normalize(f, AliasMap(), VM) === f
    end

    @testset "sign-flipped alias output" begin
        f = @initformula :θ = 0.25          # θ = -u_r, so u_r must end up at -0.25
        n = normalize(f, AM, VM)
        @test n.outsym == [:u_r]
        @test n.derived_from === f
        d = Dict{Symbol,Float64}()
        apply_init_formulas!(d, [n])
        @test d == Dict(:u_r => -0.25)
    end

    @testset "sign-flipped alias input" begin
        f = @initformula :Vset = :θ
        n = normalize(f, AM, VM)
        @test n.sym == [:u_r]
        d = Dict{Symbol,Float64}(:u_r => 3.0)
        apply_init_formulas!(d, [n])
        @test d[:Vset] ≈ -3.0
    end

    # the formula reads an observable nobody pinned a default on; it runs anyway, from the
    # defaults of the roots
    @testset "expansion of a genuinely algebraic obs input" begin
        f = @initformula :P = :summed * 2
        n = normalize(f, AM, VM)
        @test Set(n.sym) == Set([:u_r, :u_i, :Vset])
        d = Dict{Symbol,Float64}(:u_r => 1.0, :u_i => 2.0, :Vset => 4.0)
        apply_init_formulas!(d, [n])
        @test d[:P] ≈ 14.0
    end

    @testset "time-dependent obs input" begin
        f = @initformula :P = :timed        # timed = u_r * t
        d = Dict{Symbol,Float64}(:u_r => 3.0)
        apply_init_formulas!(d, [normalize(f, AM, VM; t=2.0)])
        @test d[:P] ≈ 6.0

        # roots resolved but t did not. This is an unresolvable *input*, just noticed one
        # layer deeper, so the caller decides: initialize_component errors ...
        err = try
            apply_init_formulas!(Dict{Symbol,Float64}(:u_r => 3.0), [normalize(f, AM, VM)])
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("depend explicitly on time", sprint(showerror, err))

        # ... and NWState skips (see the NWState testset)
        d = Dict{Symbol,Float64}(:u_r => 3.0)
        @test_nowarn apply_init_formulas!(d, [normalize(f, AM, VM)]; error_unresolvable=false)
        @test !haskey(d, :P)
    end
end

@testset "normalize: structural errors" begin
    @testset "output which is not translatable" begin
        f = InitFormula([:nl], [:u_i]) do out, u   # :nl = u_r^2, no slot to write to
            out[:nl] = u[:u_i]
        end
        @test_throws ArgumentError normalize(f, AM, VM)
    end

    @testset "two outputs collapsing onto one canonical symbol" begin
        f = InitFormula([:θ, :u_r], [:u_i]) do out, u
            out[:θ] = u[:u_i]
            out[:u_r] = u[:u_i]
        end
        err = try; normalize(f, AM, VM); catch e; e; end
        @test err isa ArgumentError
        @test occursin("u_r", sprint(showerror, err))
    end

    # neither of these overlaps on the raw symbol lists — normalization is what uncovers them
    @testset "hidden self-dependency" begin
        alias_induced = @initformula :θ = :u_r * 2        # writes -u_r, reads u_r
        err = try; normalize(alias_induced, AM, VM); catch e; e; end
        @test err isa ArgumentError
        @test occursin("[:u_r] → [:θ]", sprint(showerror, err)) # names raw and canonical

        expansion_induced = @initformula :u_r = :summed   # summed = u_r + u_i + Vset
        @test_throws ArgumentError normalize(expansion_induced, AM, VM)
    end
end

# A formula writing a settable symbol and one reading an observable of it share no raw
# symbol, so before normalization the DAG had no edge between them and the order was
# arbitrary. This is the central point of the whole exercise.
@testset "normalize: dependency through an observable" begin
    A = @initformula :u_i = 5.0
    B = @initformula :P = :summed                        # summed = u_r + u_i + Vset
    nA, nB = normalize(A, AM, VM), normalize(B, AM, VM)
    @test isdisjoint(A.sym, B.outsym) && isdisjoint(B.sym, A.outsym) # no raw edge
    @test topological_sort_formulas([nB, nA]) == [nA, nB]

    d = Dict{Symbol,Float64}(:u_r => 1.0, :Vset => 2.0)
    apply_init_formulas!(d, [nB, nA])
    @test d[:u_i] == 5.0
    @test d[:P] ≈ 8.0    # 1 + 5 + 2: B saw A's freshly written value
end

@testset "normalize(::GuessFormula)" begin
    f = @guessformula :P = :summed
    n = normalize(f, AM, VM)
    @test n.outsym == [:P]
    @test Set(n.sym) == Set([:u_r, :u_i, :Vset])

    # root lookup keeps the documented defaults-before-guesses priority
    guesses = Dict{Symbol,Float64}(:u_i => 100.0)
    apply_guess_formulas!(guesses, Dict{Symbol,Float64}(:u_r => 1.0, :u_i => 2.0, :Vset => 4.0), [n])
    @test guesses[:P] ≈ 7.0
    guesses = Dict{Symbol,Float64}(:u_i => 10.0)
    apply_guess_formulas!(guesses, Dict{Symbol,Float64}(:u_r => 1.0, :Vset => 4.0), [n])
    @test guesses[:P] ≈ 15.0

    # aliased output lands on the canonical symbol, and in `guesses` only
    g = normalize((@guessformula :θ = :u_i), AM, VM)
    @test g.outsym == [:u_r]
    guesses = Dict{Symbol,Float64}()
    defaults = Dict{Symbol,Float64}(:u_i => 4.0)
    apply_guess_formulas!(guesses, defaults, [g])
    @test guesses[:u_r] ≈ -4.0
    @test !haskey(defaults, :u_r)
end

# InitConstraints are deliberately left alone: they are evaluated against a full candidate
# state, so the observable mapping already reads an alias correctly, and `c.sym` feeds
# nothing structural (only `dim` does). This pins that the untouched path is correct.
@testset "InitConstraint needs no normalization" begin
    c = @initconstraint :θ + :u_i     # θ = -u_r
    @test :θ ∈ obssym(VM)             # ... reached through the observable mapping
    res = zeros(1)
    c(res, [-3.0, 10.0])              # the mapping hands it θ, not u_r
    @test res ≈ [7.0]

    # and it stays AD-transparent, which is what the nonlinear solve needs
    duals = [ForwardDiff.Dual(-3.0, 1.0), ForwardDiff.Dual(10.0, 0.0)]
    dres = zeros(eltype(duals), 1)
    c(dres, duals)
    @test ForwardDiff.value(only(dres)) ≈ 7.0
    @test ForwardDiff.partials(only(dres), 1) ≈ 1.0
end

@testset "attach-time validation (D2)" begin
    @test add_initformula!(copy(VM), @initformula :θ = :u_i) !== nothing      # alias target
    @test add_initformula!(copy(VM), @initformula :u_i = :summed) !== nothing # obs input
    @test_throws ArgumentError add_initformula!(copy(VM), @initformula :nl = :u_i)
    @test_throws ArgumentError add_initformula!(copy(VM), @initformula :nonexistent = :u_i)

    # observable inputs are expanded to roots, for GuessFormulas just as for InitFormulas
    @test add_guessformula!(copy(VM), @guessformula :P = :summed) !== nothing
    @test_throws ArgumentError add_guessformula!(copy(VM), @guessformula :nl = :u_i)

    # with no aliasmap this degrades to the plain "outputs must be settable" rule
    bare = copy(VM)
    delete_aliasmap!(bare)
    @test_throws ArgumentError add_initformula!(bare, @initformula :θ = :u_i)
end

@testset "normalize: provenance" begin
    f = @initformula :θ = :Vmeas + :u_i   # writes -u_r; reads Vset, u_i
    before = (copy(f.sym), copy(f.outsym), f.prettyprint)
    n = normalize(f, AM, VM)

    @testset "I3: the original is never mutated" begin
        @test (f.sym, f.outsym, f.prettyprint) == before
        @test f.derived_from === nothing
        @test n.derived_from === f
        @test n.sym !== f.sym    # nor does the copy share its vectors
    end

    @testset "show reports the original recipe plus what is actually in play" begin
        s = sprint(show, MIME"text/plain"(), n)
        @test occursin(f.prettyprint, s)
        @test occursin("(normalized: [:Vset, :u_i] → [:u_r], derived from [:Vmeas, :u_i] → [:θ])", s)
        # constructor-built objects (no prettyprint) survive too
        plain = InitFormula([:θ], [:u_i]) do out, u
            out[:θ] = u[:u_i]
        end
        @test occursin("normalized:", sprint(show, MIME"text/plain"(), normalize(plain, AM, VM)))
    end

    @testset "unresolvable roots point back at what was asked for" begin
        nb = normalize((@initformula :P = :summed), AM, VM)
        d = Dict{Symbol,Float64}(:u_r => 1.0, :Vset => 2.0) # :u_i missing
        err = try; apply_init_formulas!(d, [nb]); catch e; e; end
        @test err isa ArgumentError
        msg = sprint(showerror, err)
        @test occursin("u_i", msg)      # the root that is actually missing
        @test occursin("summed", msg)   # ... and the symbol the user wrote
        @test occursin("Defaults on observables are not consumed", msg)
    end
end

# `NWState` fills states/parameters from defaults and formulas, so it normalizes for the
# same reasons `initialize_component` does — see `_get_appropriate_dict`.
@testset "NWState" begin
    line = EdgeModel(g=AntiSymmetric((y, u_s, u_d, p, t) -> y .= 0.0), outsym=[:P], insym=[:P])
    twonode(v) = Network(path_graph(2), [v, copy(VM)], line)
    # strip the MTK-declared defaults so each testset controls its own values
    freshvm() = (v = copy(VM); foreach(s -> delete_metadata!(v, s, :default), [:u_r, :u_i, :Vset]); v)

    @testset "default on an alias reaches the canonical state" begin
        v = freshvm()
        set_default!(v, :θ, 0.3)   # θ = -u_r, and :θ is an observable
        # the dict is filtered down to sym ∪ psym, which :θ is not a member of, so this
        # only survives by landing on :u_r first
        @test NWState(twonode(v))[VIndex(1, :u_r)] ≈ -0.3
    end

    @testset "contradicting members of one alias class are caught" begin
        v = freshvm()
        set_default!(v, :θ, 0.3)
        set_default!(v, :u_r, 1.0)  # -0.3 vs 1.0, one variable, two values
        @test_throws ArgumentError NWState(twonode(v))

        v2 = freshvm()              # ... while agreeing ones merge silently
        set_default!(v2, :θ, 0.3)
        set_default!(v2, :u_r, -0.3)
        @test NWState(twonode(v2))[VIndex(1, :u_r)] ≈ -0.3
    end

    @testset "formula targeting an alias" begin
        v = freshvm()
        set_default!(v, :u_i, 2.0)
        add_initformula!(v, @initformula :θ = :u_i * 3)
        @test NWState(twonode(v))[VIndex(1, :u_r)] ≈ -6.0
    end

    @testset "formula reading an observable" begin
        v = freshvm()
        set_default!(v, :u_r, 3.0)
        add_initformula!(v, @initformula :u_i = :nl)   # nl = u_r^2
        @test NWState(twonode(v))[VIndex(1, :u_i)] ≈ 9.0
    end

    # NWState has no time to offer, so this is unresolvable like a missing default — and
    # NWState skips those rather than throwing
    @testset "unresolvable expansion is skipped, not thrown" begin
        v = freshvm()
        set_default!(v, :u_r, 3.0)
        add_initformula!(v, @initformula :u_i = :timed)  # timed = u_r * t
        @test isnan(NWState(twonode(v))[VIndex(1, :u_i)])
    end
end
