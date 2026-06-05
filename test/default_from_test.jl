using NetworkDynamics
using ModelingToolkitBase
using ModelingToolkitBase: t_nounits as t, D_nounits as Dt
using SciCompDSL
using Graphs
using Test

using NetworkDynamics: set_default_from!, get_default_from, has_default_from

# helper to build a simple vertex carrying namespaced parameters.
# `inherits` are `source => target` pairs (reads as "source feeds target").
function ivertex(name, psympairs; inherits=(), nodefault=())
    v = VertexModel(f=(dv, v, esum, p, t) -> (dv[1] = 0.0), g=1,
                    sym=[:x => 0.0], psym=first.(psympairs), name=name)
    for (s, d) in psympairs
        s in nodefault || set_default!(v, s, d)
    end
    for (source, target) in inherits
        set_default_from!(v, target, source)
    end
    v
end

@testset "same-component inherit (manual)" begin
    # device picks up busbar default
    v = ivertex(:bus, [:busbar₊Vbase => 1.5, :dev₊Vbase => 0.0];
                inherits=(:busbar₊Vbase => :dev₊Vbase,), nodefault=(:dev₊Vbase,))
    @test resolve_default_from!(v) == 1
    @test get_default(v, :dev₊Vbase) == 1.5

    # explicit-equal default is a no-op (no change, no warn)
    v = ivertex(:bus, [:busbar₊Vbase => 1.5, :dev₊Vbase => 1.5];
                inherits=(:busbar₊Vbase => :dev₊Vbase,))
    @test_logs resolve_default_from!(v)
    @test get_default(v, :dev₊Vbase) == 1.5

    # explicit-different default is preserved (+ warn when verbose)
    v = ivertex(:bus, [:busbar₊Vbase => 1.5, :dev₊Vbase => 9.0];
                inherits=(:busbar₊Vbase => :dev₊Vbase,))
    @test_logs (:warn,) resolve_default_from!(v; verbose=true)
    @test get_default(v, :dev₊Vbase) == 9.0
    # quiet by default
    v = ivertex(:bus, [:busbar₊Vbase => 1.5, :dev₊Vbase => 9.0];
                inherits=(:busbar₊Vbase => :dev₊Vbase,))
    @test_logs resolve_default_from!(v)
    @test get_default(v, :dev₊Vbase) == 9.0
end

@testset "inheritance is idempotent and tracks provenance" begin
    # dev inherits from busbar; re-resolution overwrites the *inherited* default
    # when the source changes, but never an explicit default.
    v = ivertex(:bus, [:busbar₊Vbase => 1.0, :dev₊Vbase => 0.0];
                inherits=(:busbar₊Vbase => :dev₊Vbase,), nodefault=(:dev₊Vbase,))
    resolve_default_from!(v)
    @test get_default(v, :dev₊Vbase) == 1.0

    # change the source and re-resolve -> inherited default follows
    set_default!(v, :busbar₊Vbase, 2.0)
    resolve_default_from!(v)
    @test get_default(v, :dev₊Vbase) == 2.0

    # explicit set_default! pins the value and clears provenance -> not overwritten
    set_default!(v, :dev₊Vbase, 99.0)
    set_default!(v, :busbar₊Vbase, 3.0)
    resolve_default_from!(v)
    @test get_default(v, :dev₊Vbase) == 99.0
end

@testset "raw default change invalidates inherited tag" begin
    v = ivertex(:bus, [:busbar₊Vbase => 1.0, :dev₊Vbase => 0.0];
                inherits=(:busbar₊Vbase => :dev₊Vbase,), nodefault=(:dev₊Vbase,))
    resolve_default_from!(v)
    @test get_default(v, :dev₊Vbase) == 1.0
    # bypass set_default!: change the default directly via metadata
    NetworkDynamics.set_metadata!(v, :dev₊Vbase, :default, 5.0)
    set_default!(v, :busbar₊Vbase, 2.0)
    resolve_default_from!(v)
    # the mismatched tag is treated as an explicit default and kept
    @test get_default(v, :dev₊Vbase) == 5.0
    @test !NetworkDynamics.has_metadata(v, :dev₊Vbase, :default_from_value)
end

@testset "missing source parameter always warns" begin
    v = ivertex(:bus, [:dev₊Vbase => 0.0];
                inherits=(:busbar₊Vbase => :dev₊Vbase,), nodefault=(:dev₊Vbase,))
    # warns even though verbose is false (the default)
    @test_logs (:warn,) resolve_default_from!(v)
    @test !has_default(v, :dev₊Vbase)
end

@testset "chained same-component inherit" begin
    # c inherits from b, b inherits from a; resolved iteratively in one call
    v = ivertex(:bus, [:a₊V => 3.0, :b₊V => 0.0, :c₊V => 0.0];
                inherits=(:a₊V => :b₊V, :b₊V => :c₊V), nodefault=(:b₊V, :c₊V))
    @test resolve_default_from!(v) == 2
    @test get_default(v, :b₊V) == 3.0
    @test get_default(v, :c₊V) == 3.0
end

@testset "inherit on non-parameter errors" begin
    v = ivertex(:bus, [:busbar₊Vbase => 1.0])
    # attach inherit to a state symbol -> structural misuse
    set_default_from!(v, :x, :busbar₊Vbase)
    @test_throws ArgumentError resolve_default_from!(v)
end

@testset "src/dst inherit on vertex parameter errors" begin
    v = ivertex(:bus, [:dev₊Vbase => 0.0];
                inherits=((:src, :busbar₊Vbase) => :dev₊Vbase,), nodefault=(:dev₊Vbase,))
    @test_throws ArgumentError resolve_default_from!(v)
end

@testset "cross-component inherit_src/inherit_dst" begin
    g = path_graph(2)
    v1 = ivertex(:b1, [:busbar₊Vbase => 10.0])
    v2 = ivertex(:b2, [:busbar₊Vbase => 20.0])

    e = EdgeModel(g=AntiSymmetric((e, vs, vd, p, t) -> (e[1] = 0.0)), outdim=1,
                  psym=[:srcend₊Vbase, :dstend₊Vbase], name=:line)
    set_default_from!(e, :srcend₊Vbase, (:src, :busbar₊Vbase))
    set_default_from!(e, :dstend₊Vbase, (:dst, :busbar₊Vbase))

    nw = Network(g, [v1, v2], e)
    ef = nw.im.edgem[1]
    # transformer-like: different per-end defaults
    @test get_default(ef, :srcend₊Vbase) == 10.0
    @test get_default(ef, :dstend₊Vbase) == 20.0
end

@testset "cross-component inherit on aliased edges errors / needs dealias" begin
    # star graph: one edge template shared across two edges with different vertices.
    # cross-component inheritance would assign different values to the shared instance.
    g = star_graph(3)
    mkvs() = [ivertex(Symbol(:b, i), [:busbar₊Vbase => Float64(i)]) for i in 1:3]
    mke() = (e = EdgeModel(g=AntiSymmetric((e, vs, vd, p, t) -> (e[1] = 0.0)), outdim=1,
                           psym=[:srcend₊Vbase], name=:line);
             set_default_from!(e, :srcend₊Vbase, (:dst, :busbar₊Vbase)); e)

    # aliased edge model -> not supported, must dealias
    @test_throws ArgumentError Network(g, mkvs(), mke())

    # with dealias=true each edge gets its own instance and value
    nw = Network(g, mkvs(), mke(); dealias=true)
    # the two edges go from vertex 1 to vertices 2 and 3 respectively
    @test get_default(nw.im.edgem[1], :srcend₊Vbase) == 2.0
    @test get_default(nw.im.edgem[2], :srcend₊Vbase) == 3.0
    @test nw.im.edgem[1] !== nw.im.edgem[2]
end

@testset "cyclic cross/local inheritance resolves" begin
    # edge: srcend inherits from src vertex busbar; relay inherits from srcend (local)
    g = path_graph(2)
    v1 = ivertex(:b1, [:busbar₊Vbase => 7.0])
    v2 = ivertex(:b2, [:busbar₊Vbase => 8.0])
    e = EdgeModel(g=AntiSymmetric((e, vs, vd, p, t) -> (e[1] = 0.0)), outdim=1,
                  psym=[:srcend₊Vbase, :relay₊Vbase], name=:line)
    set_default_from!(e, :srcend₊Vbase, (:src, :busbar₊Vbase))
    set_default_from!(e, :relay₊Vbase, :srcend₊Vbase)  # local, depends on cross-resolved value

    nw = Network(g, [v1, v2], e)
    ef = nw.im.edgem[1]
    @test get_default(ef, :srcend₊Vbase) == 7.0
    @test get_default(ef, :relay₊Vbase) == 7.0
end

@testset "MTK round-trip" begin
    @mtkmodel IBusbar begin
        @variables begin
            o(t), [output=true]
            i(t), [input=true]
        end
        @parameters begin
            Vbase = 5.0
        end
        @equations begin
            o ~ Vbase * i
        end
    end
    @mtkmodel IDevice begin
        @variables begin
            o(t), [output=true]
            i(t), [input=true]
        end
        @parameters begin
            Vbase, [default_from = :busbar₊Vbase]
        end
        @equations begin
            o ~ Vbase * i
        end
    end
    @mtkmodel IBus begin
        @components begin
            busbar = IBusbar()
            dev = IDevice()
        end
        @variables begin
            θ(t) = 0.0, [output=true]
            P(t), [input=true]
        end
        @equations begin
            busbar.i ~ P
            dev.i ~ P
            Dt(θ) ~ busbar.o + dev.o
        end
    end

    vm = VertexModel(IBus(name=:bus), [:P], [:θ])
    @test get_default_from(vm, :dev₊Vbase) == :busbar₊Vbase
    # resolved already in the VertexModel constructor
    @test get_default(vm, :dev₊Vbase) == 5.0
end

@testset "interaction with chk_global_parameters" begin
    g = path_graph(2)
    v1 = ivertex(:b1, [:busbar₊Vbase => 1.0, :dev₊Vbase => 0.0];
                 inherits=(:busbar₊Vbase => :dev₊Vbase,), nodefault=(:dev₊Vbase,))
    v2 = ivertex(:b2, [:busbar₊Vbase => 1.0, :dev₊Vbase => 0.0];
                 inherits=(:busbar₊Vbase => :dev₊Vbase,), nodefault=(:dev₊Vbase,))
    # mark the device Vbase as global so it is checked across the network
    set_scope!(v1, :dev₊Vbase, :global)
    set_scope!(v2, :dev₊Vbase, :global)
    e = EdgeModel(g=AntiSymmetric((e, vs, vd, p, t) -> (e[1] = 0.0)), outdim=1, pdim=0, name=:e)

    nw = Network(g, [v1, v2], e)
    # both inherited 1.0 -> consistent
    @test get_default(nw.im.vertexm[1], :dev₊Vbase) == 1.0
    @test get_default(nw.im.vertexm[2], :dev₊Vbase) == 1.0
    @test chk_global_parameters(nw) == true
end
