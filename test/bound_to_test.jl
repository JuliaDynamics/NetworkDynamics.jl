using NetworkDynamics
using NetworkDynamics: psym, obssym, get_metadata
using ModelingToolkitBase
using ModelingToolkitBase: t_nounits as t, D_nounits as D, System, @variables, @parameters, @named
using SciCompDSL: @mtkmodel
using Test

# shared per-unit models live in the ComponentLibrary (`Lib`): a bus with a `busbar` base + an
# injector `bound_to` it, and a line.
@__MODULE__() == Main ? includet(joinpath(pkgdir(NetworkDynamics), "test", "ComponentLibrary.jl")) : (const Lib = Main.Lib)

@testset "bound_to structural elimination" begin
    # On the shared bus, the injector's `S_b` is `bound_to = :busbar₊S_b`: it must leave `psym`
    # (it is no longer a free parameter) and reappear as an observable equal to the busbar base.
    vb = Lib.pu_bus_vertex()
    @test :S_b ∉ psym(vb)
    @test :busbar₊S_b ∈ psym(vb)
    @test :S_b ∈ obssym(vb)
    @test string(only(get_metadata(vb, :observed))) == "S_b ~ busbar₊S_b"
end

@testset "bound_to on both simplification backends" begin
    # The elimination must happen on ND's own reduction (`mtkcompile=false`, the default) and on
    # MTK's `mtkcompile` alike — both push the binding into the observed equations.
    for mtkcompile in (false, true)
        vm = Lib.pu_bus_vertex(; mtkcompile)
        @test :S_b ∉ psym(vm)
        @test string(only(get_metadata(vm, :observed))) == "S_b ~ busbar₊S_b"
    end
end

@testset "bound_to on the network: base propagates as an observable" begin
    # The real use case, end to end at the structural/observable level: the per-unit base lives
    # once (as `busbar₊S_b`), every bound device reads it as an observable, and changing the base
    # moves all of them — while a second, independent bus is untouched.
    nw = Lib.pu_network()
    s = NWState(nw)
    @test s[VIndex(1, :S_b)] == s[VIndex(1, :busbar₊S_b)] == 100.0

    s[VIndex(1, :busbar₊S_b)] = 250.0
    @test s[VIndex(1, :S_b)] == 250.0            # bound observable follows the base
    @test s[VIndex(2, :S_b)] == 100.0            # the other bus is independent
end

@mtkmodel MTKPuBus begin
    @components begin
        busbar = Lib.busbar_sys()
    end
    @variables begin
        u(t) = 1.0
        i(t)
        o(t)
    end
    @parameters begin
        P = 1.0
        S_b, [bound_to = :busbar₊S_b]
    end
    @equations begin
        D(u) ~ -u + i + P / S_b
        o ~ u
    end
end

@testset "bound_to via @mtkmodel (QuoteNode metadata)" begin
    # Under `@mtkmodel` the `[bound_to = :busbar₊S_b]` value arrives as a `QuoteNode` rather than
    # a bare `Symbol`; the result must be identical to the plain-`@parameters` authoring style.
    @named mb = MTKPuBus()
    vm = VertexModel(mb, [:i], [:o])
    @test :S_b ∉ psym(vm)
    @test :busbar₊S_b ∈ psym(vm)
    @test string(only(get_metadata(vm, :observed))) == "S_b ~ busbar₊S_b"
end

@testset "bound_to edge component" begin
    # `bound_to` works identically on edges: an internal line quantity aliased to a line parameter.
    @variables θ(t)=0.0 src_u(t) dst_u(t) flow(t)
    @parameters Y=0.5
    @parameters Y_pu [bound_to = :Y]
    @named line = System([D(θ) ~ Y - θ, flow ~ Y_pu * (src_u - dst_u)], t)
    em = EdgeModel(line, [:src_u], [:dst_u], AntiSymmetric([:flow]))
    @test :Y_pu ∉ psym(em)
    @test :Y ∈ psym(em)
    @test :Y_pu ∈ obssym(em)
end

@testset "bound_to error paths" begin
    # explicit default AND bound_to: contradictory intent → error naming both
    let
        @variables x(t)=0.0 i(t) o(t)
        @parameters base=100.0
        @parameters S_b=50.0 [bound_to = :base]
        @named sys = System([D(x) ~ -x * S_b + i * base, o ~ x], t)
        @test_throws ArgumentError VertexModel(sys, [:i], [:o])
    end

    # unresolvable target → error
    let
        @variables x(t)=0.0 i(t) o(t)
        @parameters S_b [bound_to = :nonexistent]
        @named sys = System([D(x) ~ -x * S_b + i, o ~ x], t)
        @test_throws ArgumentError VertexModel(sys, [:i], [:o])
    end

    # bound_to on a state/variable (not a parameter) → error
    let
        @variables x(t)=0.0 [bound_to = :base] i(t) o(t)
        @parameters base=100.0
        @named sys = System([D(x) ~ -x + i * base, o ~ x], t)
        @test_throws ArgumentError VertexModel(sys, [:i], [:o])
    end
end
