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
            normalize_valuedict(am, Dict(:a => 1.0, :b => 5.0); what=:default)
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

    @testset "on_conflict=:keepfirst drops instead of erroring (guesses)" begin
        am = AliasMap(:a => (-1.0, :b))
        # inconsistent members would throw under :error, but a soft seed just keeps one
        @test normalize_valuedict(am, Dict(:a => 1.0, :b => 5.0); on_conflict=:keepfirst) ==
              Dict(:b => 5.0)  # value on the canonical symbol wins
        # deterministic winner between two aliases: sorted-first
        am2 = AliasMap(:a => (1.0, :x), :z => (1.0, :x))
        for _ in 1:20
            @test normalize_valuedict(am2, Dict(:a => 1.0, :z => 7.0); on_conflict=:keepfirst) ==
                  Dict(:x => 1.0)
        end
        # the drop is noted under verbose, silent otherwise
        d = Dict(:a => 1.0, :b => 5.0)
        out = sprint(io -> normalize_valuedict(am, d; what=:guess, on_conflict=:keepfirst, verbose=true, io))
        @test occursin("conflicting dropped", out) && occursin("a", out) && occursin("b", out)
        @test sprint(io -> normalize_valuedict(am, d; on_conflict=:keepfirst, io)) == ""
        # agreeing members still merge silently (no spurious drop line)
        out2 = sprint(io -> normalize_valuedict(am, Dict(:a => 1.0, :b => -1.0); on_conflict=:keepfirst, verbose=true, io))
        @test !occursin("dropped", out2)
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
            normalize_valuedict(am, Dict(:θ => 0.3); what=:default, verbose=true, io)
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
                       delete_metadata!, pinned_obssyms
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

    @testset "stop_at halts expansion at pinned symbols" begin
        # scaled = 2θ, but θ is a frontier symbol now: it must survive as the root even
        # though it is neither settable nor equation-free
        roots, f = @test_nowarn generate_obs_expansion(VM, [:scaled]; stop_at=Set([:θ]))
        @test roots == [:θ]
        @test f([3.0], NaN) ≈ [6.0]

        # a stopped symbol requested directly is its own root, identity pass-through
        roots, f = generate_obs_expansion(VM, [:θ]; stop_at=Set([:θ]))
        @test roots == [:θ]
        @test f([0.5], NaN) ≈ [0.5]

        # mixed sweep: one sym stops at the frontier, the other expands past it as usual
        roots, f = generate_obs_expansion(VM, [:scaled, :nl]; stop_at=Set([:θ]))
        @test Set(roots) == Set([:θ, :u_r])
        vals = [Dict(:θ => -2.0, :u_r => 3.0)[r] for r in roots]
        @test f(vals, NaN) ≈ [-4.0, 9.0]

        # empty frontier is exactly the default behavior
        roots, f = generate_obs_expansion(VM, [:scaled]; stop_at=Set{Symbol}())
        @test roots == [:u_r]
        @test f([3.0], NaN) ≈ [-6.0]
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

# An observable written by an InitFormula is "pinned": it becomes an init-time dataflow
# node. Readers stop expansion at it and consume the written value instead of the defining
# equation — this is what makes backward flow through a model possible (a parent formula
# states what a child's output must be, the child's formula inverts its own equation).
@testset "pinned observables" begin
    @testset "pin classification and the write path" begin
        w = @initformula :nl = :u_r + 1      # nl = u_r^2 is a non-alias obs
        r = @initformula :Vset = :nl
        a = @initformula :θ = 0.5            # alias output — canonicalizes, no pin

        @test pinned_obssyms([w, r, a], VM) == Set([:nl])
        @test pinned_obssyms(nothing, VM) == Set{Symbol}()

        # without the pin declared, writing the obs still throws...
        @test_throws ArgumentError normalize(w, AM, VM)
        # ... with it, there is nothing left to rewrite on either formula
        @test normalize(w, AM, VM; pinned=Set([:nl])) === w
        @test normalize(r, AM, VM; pinned=Set([:nl])) === r
    end

    @testset "expansion stops at a pin along the way" begin
        h = @initformula :Vset = :scaled     # scaled = 2θ
        n = normalize(h, AM, VM; pinned=Set([:θ]))
        @test n.sym == [:θ]
    end

    @testset "writer → reader dataflow through the pin" begin
        w = @initformula :nl = :u_r + 1
        r = @initformula :Vset = :nl
        pins = pinned_obssyms([w, r], VM)
        nfs = [normalize(f, AM, VM; pinned=pins) for f in [r, w]]  # wrong order on purpose
        @test only(topological_sort_formulas(nfs)[1].outsym) == :nl
        d = Dict{Symbol,Float64}(:u_r => 2.0)
        apply_init_formulas!(d, nfs)
        @test d[:nl] ≈ 3.0
        @test d[:Vset] ≈ 3.0
    end

    @testset "two writers of one pin are duplicate writers" begin
        w1 = @initformula :nl = :u_r + 1
        w2 = @initformula :nl = 5.0
        pins = pinned_obssyms([w1, w2], VM)
        nfs = [normalize(f, AM, VM; pinned=pins) for f in [w1, w2]]
        @test_throws ArgumentError topological_sort_formulas(nfs)
    end

    @testset "GuessFormulas pin too, as hints" begin
        # a guess formula may write a pin (it lands in the guesses dict) and read one; the
        # frontier for guess formulas is the union of init and guess pins
        gw = @guessformula :nl = :u_r + 1
        gr = @guessformula :Vset = :nl
        @test pinned_obssyms([gw, gr], VM) == Set([:nl])
        @test normalize(gw, AM, VM; pinned=Set([:nl])) === gw
        @test normalize(gr, AM, VM; pinned=Set([:nl])) === gr
        @test add_guessformula!(copy(VM), gw) !== nothing

        guesses = Dict{Symbol,Float64}()
        apply_guess_formulas!(guesses, Dict(:u_r => 2.0), [gw, gr])
        @test guesses[:nl] ≈ 3.0
        @test guesses[:Vset] ≈ 3.0

        # an init pin on the same symbol shadows the guess pin: defaults take precedence
        guesses = Dict{Symbol,Float64}()
        apply_guess_formulas!(guesses, Dict(:u_r => 2.0, :nl => 10.0), [gw, gr])
        @test guesses[:Vset] ≈ 10.0
    end
end

# Reading through one's own pin is a genuine cycle — the pin stops the expansion, so the
# formula's effective inputs contain its own output. Needs an obs-of-obs, which the main
# fixture doesn't have.
@mtkmodel PinSelfDepBus begin
    @variables begin
        x(t) = 1.0
        y(t); yplus(t)
        i(t), [input=true]
        o(t), [output=true]
    end
    @parameters begin
        K = 2.0
    end
    @equations begin
        Dt(x) ~ -x + i
        y ~ K * x          # parameter factor: not an alias
        yplus ~ y + 1
        o ~ x
    end
end
@testset "self-dependency through one's own pin" begin
    @named _psd = PinSelfDepBus()
    VMsd = VertexModel(_psd, [:i], [:o])
    f = @initformula :y = :yplus - 1   # yplus stops at pinned :y → reads what it writes
    err = try; normalize(f, get_aliasmap(VMsd), VMsd; pinned=Set([:y])); catch e; e; end
    @test err isa ArgumentError
    @test occursin("depends on its own", sprint(showerror, err))
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
    # a non-alias obs output is a *pin* — attachable for InitFormulas, effective once the
    # init pipeline computes the pin-set (see the "pinned observables" testset)
    @test add_initformula!(copy(VM), @initformula :nl = :u_i) !== nothing
    @test_throws ArgumentError add_initformula!(copy(VM), @initformula :nonexistent = :u_i)

    # observable inputs are expanded to roots, for GuessFormulas just as for InitFormulas,
    # and obs outputs (pins) are attachable for both kinds
    @test add_guessformula!(copy(VM), @guessformula :P = :summed) !== nothing
    @test add_guessformula!(copy(VM), @guessformula :nl = :u_i) !== nothing

    # with no aliasmap, an obs-targeting formula is still attachable — it degrades from
    # "transported to the canonical symbol" to "pin on the obs"
    bare = copy(VM)
    delete_aliasmap!(bare)
    @test add_initformula!(bare, @initformula :θ = :u_i) !== nothing
    @test add_guessformula!(bare, @guessformula :θ = :u_i) !== nothing
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

# `initialize_component` normalizes values and formulas onto canonical settable symbols, so
# the init result does not depend on which member of an alias class a default/guess/formula
# was written against. This is the component-init counterpart of the `NWState` testset above.
@testset "initialize_component normalization" begin
    # strip the MTK-declared defaults so each case controls its own values
    freshvm() = (v = copy(VM); foreach(s -> delete_metadata!(v, s, :default), [:u_r, :u_i, :Vset]); v)
    # a steady, fully-pinned surrounding: i_r = u_r keeps Dt(u_r)=0, u_i=0 keeps Dt(u_i)=0,
    # :P is the free output solved from P ~ u_r + u_r^2
    base!(v) = (set_default!(v, :i_r, 0.5); set_default!(v, :u_i, 0.0);
                set_default!(v, :Vset, 1.0); set_guess!(v, :P, 0.0); v)
    eq(a, b) = keys(a) == keys(b) && all(isapprox(a[k], b[k]; atol=1e-7) for k in keys(a))

    # For the sign-flipped alias pair θ = -u_r, the same scenario carried once on the
    # canonical :u_r and once on the alias :θ must init identically — one case per carrier.
    @testset "placement invariance: a default" begin
        canon = let v = base!(freshvm()); set_default!(v, :u_r, 0.5); initialize_component(v; verbose=false) end
        alias = let v = base!(freshvm()); set_default!(v, :θ, -0.5); initialize_component(v; verbose=false) end
        @test eq(canon, alias)
        @test canon[:u_r] ≈ 0.5              # landed on the canonical symbol, not :θ
        @test !haskey(alias, :θ)
    end

    @testset "placement invariance: a guess" begin  # u_r free, solved to i_r=0.5
        canon = let v = base!(freshvm()); set_guess!(v, :u_r, 0.1);  initialize_component(v; verbose=false) end
        alias = let v = base!(freshvm()); set_guess!(v, :θ, -0.1);   initialize_component(v; verbose=false) end
        @test eq(canon, alias)
    end

    @testset "placement invariance: a formula output" begin
        canon = let v = base!(freshvm()); add_initformula!(v, @initformula :u_r = 0.5);  initialize_component(v; verbose=false) end
        alias = let v = base!(freshvm()); add_initformula!(v, @initformula :θ = -0.5);   initialize_component(v; verbose=false) end
        @test eq(canon, alias)
    end

    @testset "placement invariance: a formula input" begin
        # the formula reads the voltage through either member (adjusting for the factor) and
        # computes the same u_i; the read must resolve to the one canonical value
        mk(read) = let v = base!(freshvm())
            set_default!(v, :u_r, 0.5); delete_metadata!(v, :u_i, :default)
            add_initformula!(v, read); initialize_component(v; verbose=false)
        end
        @test eq(mk(@initformula :u_i = :u_r - 0.5), mk(@initformula :u_i = -:θ - 0.5))
    end

    @testset "placement invariance: a constraint input" begin
        # InitConstraints are not normalized, but reach an alias through the observable
        # mapping, so placement invariance holds there too (u_r free + guess, constraint pins)
        mk(con) = let v = base!(freshvm())
            set_guess!(v, :u_r, 0.1); add_initconstraint!(v, con); initialize_component(v; verbose=false)
        end
        @test eq(mk(@initconstraint :u_r - 0.5), mk(@initconstraint :θ + 0.5))
    end

    # An override may name any member of the class and must beat the value the class holds,
    # whichever member the metadata used — one case per carrier, as above.
    @testset "placement invariance: a default override" begin
        # the metadata value is the wrong one (u_r must equal i_r = 0.5 for steady state),
        # so only an override that lands on the class gets this to init at all
        mk(ov) = let v = base!(freshvm()); set_default!(v, :u_r, 0.1)
            initialize_component(v; default_overrides=ov, verbose=false)
        end
        canon = mk(Dict(:u_r => 0.5))
        alias = mk(Dict(:θ => -0.5))
        @test eq(canon, alias)
        @test alias[:u_r] ≈ 0.5    # the override, not the 0.1 from metadata
    end

    @testset "guess override on an alias beats the metadata guess" begin
        # guesses only seed the solve, so the working dict is the only place the override
        # is observable — the solution is the same either way
        g = Ref{Any}()
        v = base!(freshvm()); set_guess!(v, :u_r, 0.1)
        initialize_component(v; guess_overrides=Dict(:θ => -0.7), verbose=false, _final_guesses=g)
        @test g[][:u_r] ≈ 0.7
        @test !haskey(g[], :θ)
    end

    @testset "a removal marker reaches the class through any member" begin
        d = Ref{Any}()
        v = base!(freshvm()); set_default!(v, :u_r, 0.5); set_guess!(v, :u_r, 0.1)
        initialize_component(v; default_overrides=Dict(:θ => nothing), verbose=false, _final_defaults=d)
        @test !haskey(d[], :u_r)   # the default is gone, u_r is free and solved for again
        @test !haskey(d[], :θ)
    end

    @testset "write-back: alias readable as factor * canonical" begin
        v = base!(freshvm()); set_default!(v, :u_r, 0.5)
        initialize_component!(v; verbose=false)
        @test get_initial_state(v, :u_r) ≈ 0.5
        @test get_initial_state(v, :θ) ≈ -1 * get_initial_state(v, :u_r)      # θ = -u_r
        @test get_initial_state(v, :scaled) ≈ -2 * get_initial_state(v, :u_r) # scaled = -2u_r
    end

    @testset "t is threaded into the expansion" begin
        # the formula targets the parameter :Vset (absent from the dynamics) so it cannot
        # disturb the residual; :timed = u_r * t
        tbase!(v) = (set_default!(v, :i_r, 3.0); set_default!(v, :u_i, 0.0);
                     set_default!(v, :u_r, 3.0); set_guess!(v, :P, 0.0);
                     delete_metadata!(v, :Vset, :default); v)

        v = tbase!(freshvm()); add_initformula!(v, @initformula :Vset = :timed)
        @test initialize_component(v; verbose=false, t=2.0)[:Vset] ≈ 6.0   # 3 * 2

        # no time given (t=NaN): the expanded input is unresolvable, so it errors here (unlike
        # the NWState path which skips)
        v2 = tbase!(freshvm()); add_initformula!(v2, @initformula :Vset = :timed)
        err = try; initialize_component(v2; verbose=false); catch e; e end
        @test err isa Exception
        @test occursin("depend explicitly on time", sprint(showerror, err))
    end

    @testset "default on a non-alias observable is not consumed" begin
        # :nl = u_r^2 is a genuine observable, not an alias; a default pinned on it is not a
        # formula input, and with the root :u_r missing the formula errors with the migration
        # instruction naming the roots
        v = freshvm()
        set_default!(v, :i_r, 0.5); set_default!(v, :u_i, 0.0)
        set_default!(v, :Vset, 1.0); set_guess!(v, :P, 0.0)
        set_default!(v, :nl, 4.0)   # pinned on the observable itself
        add_initformula!(v, @initformula :Vset = :nl)
        err = try; initialize_component(v; verbose=false); catch e; e end
        @test err isa ArgumentError
        msg = sprint(showerror, err)
        @test occursin("roots of the originally requested [:nl]", msg)
        @test occursin("provide defaults for the roots", msg)
    end

    # I1: a hand-built, non-MTK component has no aliasmap and no observable-input formulas, so
    # normalization is a no-op and init is unaffected.
    @testset "empty-map component is unaffected (I1)" begin
        v = VertexModel(f=(du, u, in, p, t) -> du .= .-u .+ in, g=1:1,
                        sym=[:x], insym=[:i], outsym=[:o], name=:bare)
        @test isempty(get_aliasmap(v))
        set_default!(v, :i, 0.7); set_guess!(v, :x, 0.0); set_guess!(v, :o, 0.0)
        res = initialize_component(v; verbose=false)
        @test res[:x] ≈ 0.7 && res[:o] ≈ 0.7
    end
end

####
#### Pinned observables in the full init pipeline: a PI-flavored vertex where the
#### backward-init chain runs *through* an observable. The parent-side knowledge "in steady
#### state the PI output must hold the plant" pins `y`; the child-side knowledge "my
#### integrator state follows from my output" inverts the PI equation. Neither formula
#### works without the other, and together they determine everything — no nonlinear solve.
####
@mtkmodel PinPipelineBus begin
    @variables begin
        x(t), [guess=0]      # integrator state
        v(t) = 1.0           # plant state
        err(t); y(t)
        i(t), [input=true]
        o(t), [output=true]
    end
    @parameters begin
        Kp = 20.0
        Ki = 5.0
        K = 2.0
        vref = 1.0
    end
    @equations begin
        Dt(x) ~ err
        err ~ vref - v
        y ~ Kp*err + Ki*x
        Dt(v) ~ y - K*v + i
        o ~ v
    end
end
@testset "pinned observables: init pipeline" begin
    @named _ppb = PinPipelineBus()
    PVM = VertexModel(_ppb, [:i], [:o])
    seeds = Dict(:i => 0.5, :o => 1.0)

    pin_y  = @initformula :y = :K * :v - :i          # steady state of the plant equation
    calc_x = @initformula :x = (:y - :Kp * :err) / :Ki   # PI equation solved for x

    @testset "backward flow through the pin, zero free variables" begin
        io = IOBuffer()
        state = initialize_component(PVM;
            default_overrides=seeds,
            additional_initformula=[pin_y, calc_x],
            verbose=true, io)
        out = String(take!(io))
        @test state[:x] ≈ 0.3     # (1.5 - 20*0) / 5
        @test state[:v] ≈ 1.0
        @test !haskey(state, :y)  # the pin is init-time scratch, never part of the state
        @test occursin("(pinned observable)", out)
        @test occursin("No free variables!", out)
    end

    @testset "a wrong pin with a reader fails the residual check" begin
        bad_pin = @initformula :y = :K * :v - :i + 0.1
        @test_throws NetworkDynamics.ComponentInitError initialize_component(PVM;
            default_overrides=seeds,
            additional_initformula=[bad_pin, calc_x], verbose=false)
    end

    @testset "a wrong unread pin trips the consistency warning" begin
        # x is seeded manually, so nothing consumes the pin; the recomputed observable
        # disagrees with the asserted value and the post-solve check names the origin
        bad_pin = @initformula :y = :K * :v - :i + 0.1
        io = IOBuffer()
        initialize_component(PVM;
            default_overrides=merge(seeds, Dict(:x => 0.3)),
            additional_initformula=[bad_pin],
            verbose=false, tol=1.0, io)   # high tol: the warning must fire on its own
        out = String(take!(io))
        @test occursin("pinned by InitFormula", out)
        @test occursin("differ from their specified values", out)
    end

    @testset "NWState path applies the chain and drops the scratch value" begin
        cfm = copy(PVM)
        add_initformula!(cfm, pin_y); add_initformula!(cfm, calc_x)
        set_default!(cfm, :i, 0.5); set_default!(cfm, :o, 1.0)
        d = NetworkDynamics._get_appropriate_dict(nothing, cfm; guess=true,
                                                  apply_formulas=true, verbose=false)
        @test d[:x] ≈ 0.3
        @test !haskey(d, :y)
    end

    @testset "the same chain as guess formulas seeds the solve without committing" begin
        # spelled as guesses, the backward chain does not eliminate the free variables —
        # it starts the nonlinear solve at (what happens to be) the exact solution. The
        # library-author pattern: back-computing guess formulas are a safe default even
        # when they might be only approximately right.
        pin_y_g  = @guessformula :y = :K * :v - :i
        calc_x_g = @guessformula :x = (:y - :Kp * :err) / :Ki
        io = IOBuffer()
        state = initialize_component(PVM;
            default_overrides=seeds,
            additional_guessformula=[pin_y_g, calc_x_g],
            verbose=true, io)
        out = String(take!(io))
        @test occursin("(pinned observable)", out)   # marked in the guess rows
        @test occursin("NonlinearLeastSquaresProblem", out)  # x stays free, unlike initformula
        @test state[:x] ≈ 0.3 atol=1e-8
        @test !haskey(state, :y)
    end

    # :timed = u_r * t on the AliasNormBus fixture. The post-solve consistency check
    # recomputes the observable at the init-time `t`, so a time-dependent pin is verifiable
    # exactly when init ran at a concrete `t`.
    freshtimed() = begin
        v = copy(VM)
        set_default!(v, :i_r, 1.0); set_default!(v, :u_r, 1.0); set_default!(v, :u_i, 0.0)
        set_default!(v, :Vset, 1.0); set_default!(v, :P, 2.0)
        v
    end

    @testset "at t=NaN a time-dependent pin cannot be checked (no warning)" begin
        # recompute is u_r * NaN = NaN, which must read as "unverifiable", not a contradiction
        io = IOBuffer()
        initialize_component(freshtimed();
            additional_initformula=[@initformula :timed = 99.0],  # wrong, but unverifiable at t=NaN
            verbose=false, t=NaN, io)
        @test !occursin("WARNING", String(take!(io)))
    end

    @testset "at a concrete t a consistent time-dependent pin does not warn" begin
        # u_r*t at u_r=1, t=2 is 2.0, which matches the pinned 2*u_r
        io = IOBuffer()
        initialize_component(freshtimed();
            additional_initformula=[@initformula :timed = 2 * :u_r],
            verbose=false, t=2.0, io)
        @test !occursin("WARNING", String(take!(io)))
    end

    @testset "at a concrete t an inconsistent time-dependent pin warns" begin
        # pinned 2*u_r+0.5 = 2.5 disagrees with the recomputed u_r*t = 2.0 at t=2
        io = IOBuffer()
        initialize_component(freshtimed();
            additional_initformula=[@initformula :timed = 2 * :u_r + 0.5],
            verbose=false, t=2.0, io)
        out = String(take!(io))
        @test occursin("pinned by InitFormula", out)
        @test occursin("differ from their specified values", out)
    end
end

####
#### Integration test on a realistic composite: a synchronous machine with an AVR and a
#### governor, wired through a busbar. The field voltage `Efd` and mechanical power `Pm`
#### cross the machine↔AVR and machine↔governor boundaries, so each becomes an alias class
#### whose settable survivor is the control state (`avr₊Efd`, `gov₊Pm`). A machine back-init
#### formula produces them from the terminal condition; the AVR/governor back-init formulas
#### consume them. The whole point: those formulas share no *raw* symbol with the machine
#### one, yet must sort after it once normalization collapses the alias classes.
####
@mtkmodel GenBusBase begin
    @variables begin
        u_r(t)=1.0, [output=true]
        u_i(t)=0.0, [output=true]
        i_r(t), [input=true, guess=1.0]
        i_i(t), [input=true, guess=0.0]
        u_mag(t)
    end
    @equations begin
        u_mag ~ sqrt(u_r^2 + u_i^2)
    end
end
@mtkmodel GenMachine begin
    @variables begin
        u_r(t)=1.0, [output=true]
        u_i(t)=0.0, [output=true]
        i_r(t), [input=true, guess=1.0]
        i_i(t), [input=true, guess=0.0]
        θ(t), [guess=0.0]
        ω(t), [guess=0.0]
        Pel(t)
        Efd(t), [input=true, guess=1.0]   # field voltage, from the AVR
        Pm(t), [input=true, guess=1.0]    # mechanical power, from the governor
    end
    @parameters begin
        M = 1.0
        D = 0.1
    end
    @equations begin
        Pel ~ u_r*i_r + u_i*i_i
        Dt(θ) ~ ω
        Dt(ω) ~ 1/M*(Pm - D*ω + Pel)
        u_r ~ Efd*cos(θ)                  # stylized Park relation (algebraic, not an alias)
        u_i ~ Efd*sin(θ)
    end
end
@mtkmodel GenAVR begin
    @variables begin
        Efd(t), [output=true, guess=1.0]
        Vmeas(t), [input=true]
    end
    @parameters begin
        Ka = 10.0
        Vref, [guess=1.0]
    end
    @equations begin
        Dt(Efd) ~ Ka*(Vref - Vmeas) - Efd
    end
end
@mtkmodel GenGov begin
    @variables begin
        Pm(t), [output=true, guess=1.0]
        ω_meas(t), [input=true]
    end
    @parameters begin
        Kg = 5.0
        Pref, [guess=1.0]
    end
    @equations begin
        Dt(Pm) ~ Pref - Kg*ω_meas - Pm
    end
end
@mtkmodel GenBusModel begin
    @extend GenBusBase()
    @components begin
        machine = GenMachine()
        avr = GenAVR()
        gov = GenGov()
    end
    @equations begin
        machine.u_r ~ u_r
        machine.u_i ~ u_i
        0 ~ i_r - machine.i_r
        0 ~ i_i - machine.i_i
        machine.Efd ~ avr.Efd     # alias class {machine₊Efd, avr₊Efd}, survivor avr₊Efd
        machine.Pm ~ gov.Pm       # alias class {machine₊Pm, gov₊Pm}, survivor gov₊Pm
        avr.Vmeas ~ u_mag
        gov.ω_meas ~ machine.ω
    end
end

@testset "Step 5a: machine + AVR + governor stacking" begin
    @named genbus = GenBusModel()
    GB = VertexModel(genbus, [:i_r, :i_i], [:u_r, :u_i]; verbose=false)
    am = get_aliasmap(GB)

    # Efd and Pm survive as the control states, the machine-side names are their aliases
    @test am[:machine₊Efd] == (1.0, :avr₊Efd)
    @test am[:machine₊Pm]  == (1.0, :gov₊Pm)

    # the machine back-init produces Efd/Pm/θ from the terminal condition; the controls
    # back-compute their setpoints from Efd/Pm
    fM = @initformula begin
        :machine₊Efd = sqrt(:u_r^2 + :u_i^2)
        :machine₊θ   = atan(:u_i, :u_r)
        :machine₊Pm  = -(:u_r*:i_r + :u_i*:i_i)
    end
    fA = @initformula :avr₊Vref = :avr₊Efd/:avr₊Ka + :u_mag
    fG = @initformula :gov₊Pref = :gov₊Pm + :gov₊Kg*:machine₊ω

    # clear the machine-side voltage aliases so a terminal seed is unambiguous
    seed!(v; ur=1.0, ui=0.2, ir=0.5, ii=-0.1) = begin
        foreach(s -> delete_metadata!(v, s, :default), [:u_r, :machine₊u_r, :u_i, :machine₊u_i])
        set_default!(v, :u_r, ur); set_default!(v, :u_i, ui)
        set_default!(v, :i_r, ir); set_default!(v, :i_i, ii)
        set_default!(v, :machine₊ω, 0.0)
        v
    end

    @testset "back-init sorts before the controls after normalization" begin
        nM, nA, nG = normalize(fM, am, GB), normalize(fA, am, GB), normalize(fG, am, GB)
        # no shared *raw* symbol → before normalization the DAG has no edge and order is free
        @test isdisjoint(fA.sym, fM.outsym) && isdisjoint(fG.sym, fM.outsym)
        # normalization collapses the alias classes, so both controls now depend on M
        @test :avr₊Efd in nA.sym && :gov₊Pm in nG.sym
        @test first(topological_sort_formulas([nA, nG, nM])).outsym == nM.outsym
    end

    @testset "end-to-end init reproduces the hand-computed reference" begin
        v = seed!(copy(GB))
        add_initformula!(v, fM); add_initformula!(v, fA); add_initformula!(v, fG)
        s = initialize_component(v; verbose=false)
        ur, ui, ir, ii = 1.0, 0.2, 0.5, -0.1
        Efd = sqrt(ur^2 + ui^2); θ = atan(ui, ur); Pm = -(ur*ir + ui*ii)
        @test s[:avr₊Efd]  ≈ Efd
        @test s[:machine₊θ] ≈ θ
        @test s[:gov₊Pm]   ≈ Pm
        @test s[:avr₊Vref] ≈ Efd/10.0 + Efd        # Ka=10, Vmeas=u_mag=Efd
        @test s[:gov₊Pref] ≈ Pm                     # Kg=5, ω=0
        @test isapprox(s[:machine₊ω], 0.0; atol=1e-8)
    end

    # The Park relation u_r ~ Efd·cos(θ) is genuinely algebraic, not an alias — aliasing
    # handles the renamings, but this relation is enforced as a residual, which then acts as a
    # consistency check on the back-init: an inconsistent angle is caught rather than silently
    # accepted.
    @testset "wrong back-init is caught by the algebraic Park residual" begin
        fMbad = @initformula begin
            :machine₊Efd = sqrt(:u_r^2 + :u_i^2)
            :machine₊θ   = atan(:u_i, :u_r) + 0.5   # inconsistent with the seeded voltage
            :machine₊Pm  = -(:u_r*:i_r + :u_i*:i_i)
        end
        v = seed!(copy(GB))
        add_initformula!(v, fMbad); add_initformula!(v, fA); add_initformula!(v, fG)
        @test_throws NetworkDynamics.ComponentInitError initialize_component(v; verbose=false, warn=false)
    end

    @testset "duplicate writer into one alias class via different members" begin
        v = seed!(copy(GB))
        add_initformula!(v, fM)                              # writes :machine₊Efd
        add_initformula!(v, @initformula :avr₊Efd = 1.0)     # writes the survivor directly
        add_initformula!(v, fA); add_initformula!(v, fG)
        err = try; initialize_component(v; verbose=false); catch e; e end
        @test err isa ArgumentError
        @test occursin("avr₊Efd", sprint(showerror, err))
    end

    @testset "placement invariance across subcomponent levels" begin
        # the voltage seed carried on the bus symbol vs. the machine-level alias member
        initwith(seedsym) = begin
            v = copy(GB)
            foreach(s -> delete_metadata!(v, s, :default),
                    [:u_r, :machine₊u_r, :u_i, :machine₊u_i])
            set_default!(v, seedsym, 0.9)
            set_default!(v, :u_i, 0.2); set_default!(v, :i_r, 0.5); set_default!(v, :i_i, -0.1)
            set_default!(v, :machine₊ω, 0.0)
            add_initformula!(v, fM); add_initformula!(v, fA); add_initformula!(v, fG)
            initialize_component(v; verbose=false)
        end
        a = initwith(:u_r); b = initwith(:machine₊u_r)
        @test keys(a) == keys(b)
        @test all(isapprox(a[k], b[k]; atol=1e-9) for k in keys(a))
    end
end

####
#### Alias classes *without* a settable member: two names for one observable, created by a
#### connection between two subcomponents. `producer₊out_u` is defined by a many-to-one sum,
#### so neither it nor `consumer₊inp_u` has a storage slot — the class has no settable
#### survivor and is invisible to the AliasMap today. Yet the two names are provably one
#### value (MTK emits the identity `consumer₊inp_u ~ producer₊out_u` itself), so a pin on
#### either of them must be visible to a reader of the other: which of two interchangeable
#### names a formula happens to be written against must not decide whether the model
#### initializes.
####
@mtkmodel ObsAliasProducer begin
    @variables begin
        x1(t), [guess=1]
        x2(t), [guess=1]
        out_u(t)
    end
    @parameters begin
        K1 = 0.5
        K2 = 0.5
        r1, [guess=1]
        r2, [guess=1]
    end
    @equations begin
        Dt(x1) ~ r1 - x1
        Dt(x2) ~ r2 - x2
        out_u ~ K1*x1 + K2*x2   # many-to-one: a pure alias of no single state
    end
end
@mtkmodel ObsAliasConsumer begin
    @variables begin
        y(t) = 1.0
        inp_u(t)
    end
    @equations begin
        Dt(y) ~ inp_u - 2*y     # steady state ⇔ the demand `inp_u = 2y`
    end
end
@mtkmodel ObsAliasBus begin
    @components begin
        producer = ObsAliasProducer()
        consumer = ObsAliasConsumer()
    end
    @variables begin
        o(t), [output=true, guess=1]
    end
    @equations begin
        producer.out_u ~ consumer.inp_u   # the connection: an observable-only alias class
        o ~ consumer.y
    end
end

@testset "observable-only alias class" begin
    @named _oab = ObsAliasBus()
    OAB = VertexModel(_oab, [], [:o]; verbose=false)

    # the producer's own backward knowledge: "at rest my states follow from my output"
    invert = @initformula begin
        :producer₊x1 = :producer₊out_u / (:producer₊K1 + :producer₊K2)
        :producer₊x2 = :producer₊out_u / (:producer₊K1 + :producer₊K2)
    end
    # the same demand, once per member of the class
    demand_producer_name = @initformula :producer₊out_u = 2.0 * :consumer₊y
    demand_consumer_name = @initformula :consumer₊inp_u = 2.0 * :consumer₊y

    @testset "the class is detected, with the terminal observable as canonical" begin
        # MTK normalizes every reference onto `producer₊out_u` (it carries the defining
        # equation) and emits the leaf as a plain identity — that terminal is the only
        # possible canonical, since expansion bottoms out there by construction
        @test obssym(OAB) == [:producer₊out_u, :consumer₊inp_u]
        @test get_aliasmap(OAB)[:consumer₊inp_u] == (1.0, :producer₊out_u)
    end

    @testset "the pin is one node under both names" begin
        # whichever member is written, the frontier is expressed in canonical names, so the
        # producer's reader of `producer₊out_u` finds it
        @test pinned_obssyms([demand_producer_name, invert], OAB) == Set([:producer₊out_u])
        @test pinned_obssyms([demand_consumer_name, invert], OAB) == Set([:producer₊out_u])
    end

    # y = 1 ⇒ demanded out_u = 2 ⇒ x1 = x2 = 2/(K1+K2) = 2 ⇒ r1 = r2 = 2, residual 0.
    reference = Dict(:producer₊x1 => 2.0, :producer₊x2 => 2.0,
                     :producer₊r1 => 2.0, :producer₊r2 => 2.0,
                     :consumer₊y  => 1.0, :o => 1.0)

    @testset "both spellings of the demand initialize identically" begin
        for demand in (demand_producer_name, demand_consumer_name)
            state = initialize_component(OAB;
                additional_initformula=[demand, invert], verbose=false)
            for (s, v) in reference
                @test state[s] ≈ v
            end
            @test !haskey(state, :producer₊out_u)   # the pin is init-time scratch
            @test !haskey(state, :consumer₊inp_u)
        end
    end

    @testset "the boundary case still works (producer output settable)" begin
        # compiled standalone, `out_u` is a real output: settable, hence a root, which is
        # why the identical formula has always resolved here
        @named _oap = ObsAliasProducer()
        P = VertexModel(_oap, [], [:out_u]; verbose=false)
        @test :out_u ∉ obssym(P)   # it is a real output here, not an observable
        state = initialize_component(P;
            default_overrides=Dict(:out_u => 2.0),
            additional_initformula=[@initformula begin
                :x1 = :out_u / (:K1 + :K2)
                :x2 = :out_u / (:K1 + :K2)
            end], verbose=false)
        @test state[:x1] ≈ 2.0 && state[:x2] ≈ 2.0
    end

    @testset "two members of the class are two writers of one node" begin
        err = try
            initialize_component(OAB;
                additional_initformula=[demand_producer_name, demand_consumer_name, invert],
                verbose=false)
        catch e; e end
        @test err isa ArgumentError
        @test occursin("producer₊out_u", sprint(showerror, err))
    end
end
