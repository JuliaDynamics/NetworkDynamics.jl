using NetworkDynamics
using ModelingToolkitBase
using ModelingToolkitBase: t_nounits as t, D_nounits as Dt
using SciCompDSL
using Graphs
using Test

@testset "scope metadata basics" begin
    # plain VertexModel + set_scope!
    v = VertexModel(f=(dv, v, esum, p, t) -> (dv[1] = 0.0), g=1,
                    sym=[:x => 0.0], psym=[:Vbase => 1.0], name=:plainv)
    @test !has_scope(v, :Vbase)
    set_scope!(v, :Vbase, :global)
    @test has_scope(v, :Vbase)
    @test get_scope(v, :Vbase) == :global
    delete_scope!(v, :Vbase)
    @test !has_scope(v, :Vbase)
    set_scope!(v, :Vbase, :global)
    strip_scopes!(v)
    @test !has_scope(v, :Vbase)
end

@testset "scope from @mtkmodel definition" begin
    @mtkmodel ScopedNode begin
        @variables begin
            θ(t) = 0.0, [output = true]
            P(t), [input = true]
        end
        @parameters begin
            Vbase = 1.0, [scope = :global]
            kp = 2.0, [scope = :device]
            loc = 3.0
        end
        @equations begin
            Dt(θ) ~ kp * (Vbase + P + loc)
        end
    end
    vm = VertexModel(ScopedNode(name=:scoped), [:P], [:θ])

    @test get_scope(vm, :Vbase) == :global
    @test get_scope(vm, :kp) == :device
    @test !has_scope(vm, :loc)
end

# helper to build a simple vertex with given scoped parameters
function scopedvertex(name, psympairs, scopes)
    v = VertexModel(f=(dv, v, esum, p, t) -> (dv[1] = 0.0), g=1,
                    sym=[:x => 0.0], psym=psympairs, name=name)
    for (s, sc) in scopes
        set_scope!(v, s, sc)
    end
    v
end

@testset "global scope consistency on Network" begin
    g = path_graph(2)
    v1 = scopedvertex(:v1, [:Vbase => 1.0, :loc => 5.0], [:Vbase => :global, :loc => :local])
    v2 = scopedvertex(:v2, [:Vbase => 1.0, :loc => 7.0], [:Vbase => :global, :loc => :local])
    e = EdgeModel(g=AntiSymmetric((e, vs, vd, p, t) -> (e[1] = 0.0)), outdim=1, pdim=0, name=:e)

    nw = Network(g, [v1, v2], e)
    @test chk_global_parameters(nw) == true

    # make the global parameter inconsistent
    set_default!(v2, :Vbase, 2.0)
    nw2 = Network(g, [v1, v2], e)
    @test (@test_logs (:warn,) chk_global_parameters(nw2)) == false
    @test chk_global_parameters(nw2; verbose=false) == false
end

@testset "device scope consistency within component" begin
    g = path_graph(2)
    # namespaced parameters sharing the basename `kp`
    vok = scopedvertex(:vok, [Symbol("a₊kp") => 1.0, Symbol("b₊kp") => 1.0],
                       [Symbol("a₊kp") => :device, Symbol("b₊kp") => :device])
    vbad = scopedvertex(:vbad, [Symbol("a₊kp") => 1.0, Symbol("b₊kp") => 9.0],
                        [Symbol("a₊kp") => :device, Symbol("b₊kp") => :device])
    e = EdgeModel(g=AntiSymmetric((e, vs, vd, p, t) -> (e[1] = 0.0)), outdim=1, pdim=0, name=:e)

    nwok = Network(g, [vok, vok], e)
    @test chk_global_parameters(nwok) == true

    nwbad = Network(g, [vok, vbad], e)
    @test chk_global_parameters(nwbad; verbose=false) == false
end

@testset "device scope from @mtkmodel subcomponents" begin
    # two subcomponents each carrying a device-scoped parameter `kp`, which become
    # the namespaced parameters `a₊kp` and `b₊kp` on the VertexModel level
    @mtkmodel DevSub begin
        @variables begin
            o(t), [output = true]
            i(t), [input = true]
        end
        @parameters begin
            kp = 1.0, [scope = :device]
        end
        @equations begin
            o ~ kp * i
        end
    end
    @mtkmodel DevNode begin
        @components begin
            a = DevSub()
            b = DevSub()
        end
        @variables begin
            θ(t) = 0.0, [output = true]
            P(t), [input = true]
        end
        @equations begin
            a.i ~ P
            b.i ~ P
            Dt(θ) ~ a.o + b.o
        end
    end

    vmok = VertexModel(DevNode(name=:dev), [:P], [:θ])
    @test get_scope(vmok, Symbol("a₊kp")) == :device
    @test get_scope(vmok, Symbol("b₊kp")) == :device

    vmbad = VertexModel(DevNode(name=:dev), [:P], [:θ])
    set_default!(vmbad, Symbol("b₊kp"), 9.0)

    g = path_graph(2)
    e = EdgeModel(g=AntiSymmetric((e, vs, vd, p, t) -> (e[1] = 0.0)), outdim=1, pdim=0, name=:e)

    # consistent device parameters within each component
    nwok = Network(g, [vmok, vmok], e)
    @test chk_global_parameters(nwok) == true

    # inconsistent device parameters within the second component
    nwbad = Network(g, [vmok, vmbad], e)
    @test chk_global_parameters(nwbad; verbose=false) == false

    # device inconsistency is also caught automatically on ODEProblem construction
    s0 = NWState(nwbad)
    s0.v[:, :θ] .= 0.0
    @test_logs (:warn,) ODEProblem(nwbad, s0, (0.0, 1.0))
end

@testset "consistency check on NWState / NWParameter (current values)" begin
    g = path_graph(2)
    v1 = scopedvertex(:v1, [:Vbase => 1.0], [:Vbase => :global])
    v2 = scopedvertex(:v2, [:Vbase => 1.0], [:Vbase => :global])
    e = EdgeModel(g=AntiSymmetric((e, vs, vd, p, t) -> (e[1] = 0.0)), outdim=1, pdim=0, name=:e)
    nw = Network(g, [v1, v2], e)

    # defaults are consistent
    @test chk_global_parameters(nw) == true

    p = NWParameter(nw)
    @test chk_global_parameters(p) == true

    # tamper the current value of one global parameter
    p[VIndex(2, :Vbase)] = 42.0
    @test chk_global_parameters(p; verbose=false) == false

    s = NWState(nw)
    @test chk_global_parameters(s) == true
    s.p[VIndex(2, :Vbase)] = 42.0
    @test chk_global_parameters(s; verbose=false) == false
end

@testset "automatic check on ODEProblem construction" begin
    g = path_graph(2)
    v1 = scopedvertex(:v1, [:Vbase => 1.0], [:Vbase => :global])
    v2 = scopedvertex(:v2, [:Vbase => 1.0], [:Vbase => :global])
    e = EdgeModel(g=AntiSymmetric((e, vs, vd, p, t) -> (e[1] = 0.0)), outdim=1, pdim=0, name=:e)
    nw = Network(g, [v1, v2], e)

    s0 = NWState(nw)
    s0.v[:, :x] .= 0.0

    # consistent parameters -> no warning
    @test_logs ODEProblem(nw, s0, (0.0, 1.0))

    # inconsistent parameters -> warning emitted during construction
    s0.p[VIndex(2, :Vbase)] = 42.0
    @test_logs (:warn,) ODEProblem(nw, s0, (0.0, 1.0))

    # check can be disabled globally
    NetworkDynamics.CHECK_GLOBAL_PARAMETERS[] = false
    @test_logs ODEProblem(nw, s0, (0.0, 1.0))
    NetworkDynamics.CHECK_GLOBAL_PARAMETERS[] = true
end
