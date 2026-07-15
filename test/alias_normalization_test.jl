using NetworkDynamics
using NetworkDynamics: AliasMap, canonicalize, normalize_valuedict, normalize_bounds
using Test

# The primitives under test are pure functions of an AliasMap and a Dict: no component, no
# Symbolics and (deliberately, cf. I10) no ModelingToolkit is involved.

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
