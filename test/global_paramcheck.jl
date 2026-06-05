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
            kp = 2.0, [scope = :component]
            loc = 3.0
        end
        @equations begin
            Dt(θ) ~ kp * (Vbase + P + loc)
        end
    end
    vm = VertexModel(ScopedNode(name=:scoped), [:P], [:θ])

    @test get_scope(vm, :Vbase) == :global
    @test get_scope(vm, :kp) == :component
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
    @test (@test_logs (:warn, r"Inconsistent global parameter :Vbase") chk_global_parameters(nw2)) == false
    @test chk_global_parameters(nw2; verbose=false) == false
end

@testset "component scope consistency within component" begin
    g = path_graph(2)
    # namespaced parameters sharing the basename `kp`
    vok = scopedvertex(:vok, [Symbol("a₊kp") => 1.0, Symbol("b₊kp") => 1.0],
                       [Symbol("a₊kp") => :component, Symbol("b₊kp") => :component])
    vbad = scopedvertex(:vbad, [Symbol("a₊kp") => 1.0, Symbol("b₊kp") => 9.0],
                        [Symbol("a₊kp") => :component, Symbol("b₊kp") => :component])
    e = EdgeModel(g=AntiSymmetric((e, vs, vd, p, t) -> (e[1] = 0.0)), outdim=1, pdim=0, name=:e)

    nwok = Network(g, [vok, vok], e)
    @test chk_global_parameters(nwok) == true

    nwbad = Network(g, [vok, vbad], e)
    @test chk_global_parameters(nwbad; verbose=false) == false
end

@testset "component scope from @mtkmodel subcomponents" begin
    # two subcomponents each carrying a component-scoped parameter `kp`, which become
    # the namespaced parameters `a₊kp` and `b₊kp` on the VertexModel level
    @mtkmodel DevSub begin
        @variables begin
            o(t), [output = true]
            i(t), [input = true]
        end
        @parameters begin
            kp = 1.0, [scope = :component]
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
    @test get_scope(vmok, Symbol("a₊kp")) == :component
    @test get_scope(vmok, Symbol("b₊kp")) == :component

    vmbad = VertexModel(DevNode(name=:dev), [:P], [:θ])
    set_default!(vmbad, Symbol("b₊kp"), 9.0)

    g = path_graph(2)
    e = EdgeModel(g=AntiSymmetric((e, vs, vd, p, t) -> (e[1] = 0.0)), outdim=1, pdim=0, name=:e)

    # consistent component parameters within each component
    nwok = Network(g, [vmok, vmok], e)
    @test chk_global_parameters(nwok) == true

    # inconsistent component parameters within the second component
    nwbad = Network(g, [vmok, vmbad], e)
    @test chk_global_parameters(nwbad; verbose=false) == false

    # component inconsistency is also caught automatically on ODEProblem construction
    s0 = NWState(nwbad)
    s0.v[:, :θ] .= 0.0
    @test_logs (:warn, r"Inconsistent component parameter :kp") match_mode = :any begin
        ODEProblem(nwbad, s0, (0.0, 1.0))
    end
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

@testset "scope on non-parameter symbols warns" begin
    g = path_graph(2)
    e = EdgeModel(g=AntiSymmetric((e, vs, vd, p, t) -> (e[1] = 0.0)), outdim=1, pdim=0, name=:e)
    v1 = scopedvertex(:v1, [:Vbase => 1.0], [:Vbase => :global])
    v2 = scopedvertex(:v2, [:Vbase => 1.0], [:Vbase => :global])
    # attach scope to a *state* symbol, which is meaningless and should warn
    set_scope!(v1, :x, :global)

    nw = Network(g, [v1, v2], e)
    @test_logs (:warn, r"Symbol :x .* not a parameter") match_mode = :any begin
        chk_global_parameters(nw)
    end
    # still consistent for the actual parameters
    @test chk_global_parameters(nw; verbose=false) == true
end

@testset "mixed scopes for the same basename warn" begin
    e = EdgeModel(g=AntiSymmetric((e, vs, vd, p, t) -> (e[1] = 0.0)), outdim=1, pdim=0, name=:e)

    # network-wide collision: `Vbase` is :global on v1 but :local on v2
    g = path_graph(2)
    v1 = scopedvertex(:v1, [:Vbase => 1.0], [:Vbase => :global])
    v2 = scopedvertex(:v2, [:Vbase => 1.0], [:Vbase => :local])
    nw = Network(g, [v1, v2], e)
    @test_logs (:warn, r"Vbase.*mixed scopes across the network") match_mode = :any begin
        chk_global_parameters(nw)
    end
    # values are not compared across scopes -> still "consistent"
    @test chk_global_parameters(nw; verbose=false) == true

    # within-component collision: `a₊kp` is :component but `b₊kp` is :local
    vmix = scopedvertex(:vmix, [Symbol("a₊kp") => 1.0, Symbol("b₊kp") => 9.0],
                        [Symbol("a₊kp") => :component, Symbol("b₊kp") => :local])
    nwc = Network(g, [vmix, vmix], e)
    @test_logs (:warn, r"kp.*mixed scopes within") match_mode = :any begin
        chk_global_parameters(nwc)
    end

    # network-wide collision between :global and :component for the same name
    vg = scopedvertex(:vg, [:y => 1.0], [:y => :global])
    vk = scopedvertex(:vk, [:y => 1.0], [:y => :component])
    nwgk = Network(g, [vg, vk], e)
    @test_logs (:warn, r"y.*mixed scopes across the network") match_mode = :any begin
        chk_global_parameters(nwgk)
    end

    # `:component` here vs `:local` in another vertex is fine -> no warning.
    # `gx` is :global (uniform), `c` is :component in vca but :local in vcb.
    vca = scopedvertex(:vca, [:gx => 1.0, :c => 1.0], [:gx => :global, :c => :component])
    vcb = scopedvertex(:vcb, [:gx => 1.0, :c => 2.0], [:gx => :global, :c => :local])
    nwcomp = Network(g, [vca, vcb], e)
    @test_logs chk_global_parameters(nwcomp)

    # uniform scopes -> no mixed-scope warning
    vu1 = scopedvertex(:vu1, [:Vbase => 1.0], [:Vbase => :global])
    vu2 = scopedvertex(:vu2, [:Vbase => 1.0], [:Vbase => :global])
    nwu = Network(g, [vu1, vu2], e)
    @test_logs chk_global_parameters(nwu)
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
    @test_logs (:warn, r"Inconsistent global parameter :Vbase") match_mode = :any begin
        ODEProblem(nw, s0, (0.0, 1.0))
    end

    # check can be disabled globally
    NetworkDynamics.CHECK_GLOBAL_PARAMETERS[] = false
    @test_logs ODEProblem(nw, s0, (0.0, 1.0))
    NetworkDynamics.CHECK_GLOBAL_PARAMETERS[] = true
end
