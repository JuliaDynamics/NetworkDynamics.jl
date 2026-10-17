using NetworkDynamics
using NetworkDynamics: wrap_component_callbacks, get_callbacks, getcondition, getaffect, getaffect_neg,
                       condition_dim, condition_urange, condition_outrange, affect_dim,
                       collect_c_or_a_indices, AffectAccess, shortrepr,
                       _batch_condition, _scalar_member_affect, _batch_vector_affect, _gather,
                       IterativeBatches
using Graphs
using OrdinaryDiffEqTsit5
using OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqNonlinearSolve
using LinearAlgebra: Diagonal
using SciMLBase
using Chairmarks
using Test
using ModelingToolkitBase
using DiffEqCallbacks

@__MODULE__()==Main ? includet(joinpath(pkgdir(NetworkDynamics), "test", "ComponentLibrary.jl")) : (const Lib = Main.Lib)

function basenetwork()
    g = SimpleGraph([0 1 1 0 1;
                     1 0 1 1 0;
                     1 1 0 1 0;
                     0 1 1 0 1;
                     1 0 0 1 0])

    vs = [Lib.swing_mtk() for _ in 1:5];
    set_default!(vs[1], :Pmech, -1)
    set_default!(vs[2], :Pmech, 1.5)
    set_default!(vs[3], :Pmech, -1)
    set_default!(vs[4], :Pmech, -1)
    set_default!(vs[5], :Pmech, 1.5)

    ls = [Lib.line_mtk() for _ in 1:7];
    nw = Network(g, vs, ls)
    sinit = NWState(nw)
    s0 = find_fixpoint(nw)
    set_defaults!(nw, s0)
    nw
end

@testset "continuous callback batch tests" begin
    nw = basenetwork()

    tript = Float64[]
    tripi = Int[]
    cond = ComponentCondition([:P, :₋P, :srcθ], [:limit, :K]) do u, p, t
        abs(u[:P]) - p[:limit]
    end
    affect = ComponentAffect([],[:active]) do u, p, ctx
        @info "Trip line $(ctx.eidx) between $(ctx.src) and $(ctx.dst) at t=$(ctx.t)"
        push!(tript, ctx.t)
        push!(tripi, ctx.eidx)
        p[:active] = 0
    end
    cb = ContinuousComponentCallback(cond, affect)
    set_callback!.(nw.im.edgem, Ref(cb))

    batches = wrap_component_callbacks(nw);
    @test length(batches) == 1
    cbb = only(batches);

    # internally sym and psym form one list
    @test condition_dim(cbb) == 5
    @test all(affect_dim.(Ref(cbb), getaffect, 1:7) .== 1)

    @test condition_urange.(Ref(cbb), 1:length(cbb)) == [1:5,6:10,11:15,16:20,21:25,26:30,31:35]
    @test AffectAccess(cbb, getaffect).ranges == [1:1,2:2,3:3,4:4,5:5,6:6,7:7]
    @test condition_outrange.(Ref(cbb), 1:length(cbb)) == [1:1,2:2,3:3,4:4,5:5,6:6,7:7]

    @test collect_c_or_a_indices(cbb, getcondition) == collect(Iterators.flatten(collect(EIndex(i, [:P, :₋P, :srcθ, :limit, :K])) for i in 1:7))
    @test collect_c_or_a_indices(cbb, getaffect) == collect(EIndex(i, :active) for i in 1:7)

    batchcond = _batch_condition(cbb)
    out = zeros(7)
    fill!(out, NaN)
    s0 = NWState(nw)
    u = uflat(s0)
    integrator = (; p = pflat(s0))
    b = @b $batchcond($out, $u, NaN, $integrator)
    @test b.allocs == 0

    # test the preste time callback
    tripfirst = PresetTimeComponentCallback(1.0, affect) # reuse the same affect
    add_callback!(nw[EIndex(5)], tripfirst)

    # add a useless discrete callback
    useless_triggertime = Ref{Float64}(0.0)
    usless_cond = ComponentCondition([:P, :₋P, :srcθ], [:limit, :K]) do u, p, t
        t > 0.1 && iszero(useless_triggertime[])
    end
    usless_affect = ComponentAffect([], [:limit, :K]) do u, p, ctx
        @info "Usless effect triggered at $(ctx.t)"
        useless_triggertime[] = ctx.t
    end
    useless_cb = DiscreteComponentCallback(usless_cond, usless_affect)
    add_callback!(nw[EIndex(1)], useless_cb)

    s0 = NWState(nw)
    prob = ODEProblem(nw, uflat(s0), (0,6), copy(pflat(s0)))
    sol = solve(prob, Tsit5());

    @test 0.1 < useless_triggertime[] <= 1.0

    @test tripi == [5,7,4,1,3,2]
    tref = [1, 2.247676397005474, 2.502523192233235, 3.1947647115093654, 3.3380530127462587, 3.4042696241577888]
    @test maximum(abs.(tript - tref)) < 1e-5
end

@testset "show functions for callbacks" begin
    nop = (args...) -> nothing
    new = ContinuousComponentCallback(ComponentCondition(nop, [:P, :limit]), ComponentAffect(nop, [:P, :active]))
    legacy = ContinuousComponentCallback(ComponentCondition(nop, [:P], [:limit]), ComponentAffect(nop, [:P], [:active]))
    @test shortrepr(new) == shortrepr(legacy) == "(:P, :limit) affecting (:P, :active)"
    @test repr("text/plain", new) == "ContinuousComponentCallback((:P, :limit) affecting (:P, :active))"
    updown = ContinuousComponentCallback(ComponentCondition(nop, [:P]), ComponentAffect(nop, [:active]);
                                         affect_neg! = ComponentAffect(nop, [:P, :limit]))
    @test shortrepr(updown) == "(:P) affecting (:active, :P, :limit)"
    vec = VectorContinuousComponentCallback(ComponentCondition(nop, [:θ]), ComponentAffect(nop, [:ω, :D]), 2)
    @test repr("text/plain", vec) == "VectorContinuousComponentCallback((:θ) affecting (:ω, :D), len=2)"
    pt = PresetTimeComponentCallback(1.0, ComponentAffect(nop, [], [:Pmech]))
    @test shortrepr(pt) == "(:Pmech) affected at t=1.0"

    nw = basenetwork()
    e = nw.im.edgem[1]
    set_callback!(e, new)
    add_callback!(e, legacy)
    @test count("(:P, :limit) affecting (:P, :active)", repr("text/plain", e)) == 2
    @test occursin("callback", repr("text/plain", nw))

    @test delete_callbacks!(e)
    @test !has_callback(e)
    @test !delete_callbacks!(e)
end

@testset "vector callbacks" begin
    nw = basenetwork()
    u0 = zeros(dim(nw))
    p0 = NWParameter(nw)
    p0.v[:, :D] .= 1
    prob = ODEProblem(nw, u0, (0, 10.0), pflat(p0))
    sol = solve(prob, Tsit5())

    events = []
    cond = ComponentCondition([:θ, :ω], []) do out, u, p ,t
        out[1] = 0.18 - abs(u[:θ])
        out[2] = -0.2 - u[:ω]
    end
    affect = ComponentAffect([:θ, :ω],[]) do u, p, event_signs, ctx
        for event_idx in eachindex(event_signs)
            event_signs[event_idx] == 0 && continue
            push!(events, (;θ=u[:θ], ω=u[:ω], t=ctx.t, vidx=ctx.vidx, event_idx=event_idx))
            @info "Triggered event_idx $event_idx (dir $(event_signs[event_idx])) t=$(ctx.t) on $(ctx.vidx)"
        end
    end
    ccb = VectorContinuousComponentCallback(cond, affect, 2)
    set_callback!(nw.im.vertexm[1], ccb)
    set_callback!(nw.im.vertexm[2], ccb)

    cbbs = wrap_component_callbacks(nw);
    @test length(cbbs) == 1

    nwcb = get_callbacks(nw);
    prob = remake(prob, callback=nwcb);
    sol = solve(prob, Tsit5());

    # plot for interacive inspection
    # let
    #     fig = Figure();
    #     ax1 = Axis(fig[1,1])
    #     ax2 = Axis(fig[2,1])
    #     lines!(ax1, sol; idxs=vidxs(sol, 1:2, :θ))
    #     hlines!(ax1, [0.18], color=:black)
    #     hlines!(ax1, [-0.18], color=:black)
    #     lines!(ax2, sol; idxs=vidxs(sol, 1:2, :ω))
    #     hlines!(ax2, [-0.2], color=:black)
    #     for e in events
    #         color = CairoMakie.Makie.wong_colors()[e.vidx]
    #         @show color
    #         if e.event_idx ==1
    #             scatter!(ax1, [e.t], [e.θ], color=color)
    #         else
    #             scatter!(ax2, [e.t], [e.ω], color=color)
    #         end
    #     end
    #     fig
    # end

    ref_events = [
        (θ=-0.026601727368181737, ω=-0.19999999999999982, t=0.2449490689620347, vidx=1, event_idx=2)
        (θ=0.17999999999999963, ω=0.3881117850566674, t=0.6113543225776816, vidx=2, event_idx=1)
        (θ=-0.1600980341412601, ω=-0.20000000000000004, t=0.781905478960556, vidx=1, event_idx=2)
        (θ=-0.18000000000000002, ω=-0.1403736920664344, t=0.8982900942874107, vidx=1, event_idx=1)
        (θ=0.26148399109049375, ω=-0.19999999999999998, t=1.3593250291764198, vidx=2, event_idx=2)
        (θ=-0.18, ω=0.11335697475352356, t=1.4312072751759022, vidx=1, event_idx=1)
        (θ=0.18000000000000033, ω=-0.3092908878282453, t=1.6621217378351785, vidx=2, event_idx=1)
        (θ=0.06997679994525381, ω=-0.2000000000000001, t=2.0617667253585683, vidx=2, event_idx=2)
        (θ=0.17999999999999997, ω=0.13070953683042436, t=3.3585329163801663, vidx=2, event_idx=1)
        (θ=0.18000000000000013, ω=-0.09014601420085934, t=4.213717530716612, vidx=2, event_idx=1)
    ]
    function check_events(events)
        @test length(events) == length(ref_events)
        for (e, re) in zip(events, ref_events)
            @test e.vidx == re.vidx
            @test e.event_idx == re.event_idx
            @test abs(e.θ - re.θ) < 1e-4
            @test abs(e.ω - re.ω) < 1e-4
            @test abs(e.t - re.t) < 1e-4
        end
    end
    check_events(events)

    # the same scenario in the single list form, the affect also reads a parameter and an observed
    empty!(events)
    cond2 = ComponentCondition([:θ, :ω, :D]) do out, u, t
        out[1] = 0.18*u[:D] - abs(u[:θ])
        out[2] = -0.2*u[:D] - u[:ω]
    end
    affect2 = ComponentAffect([:θ, :ω, :D, :Pdamping]) do u, event_signs, ctx
        @test u[:Pdamping] ≈ -u[:D] * u[:ω]
        for i in eachindex(event_signs)
            event_signs[i] == 0 && continue
            push!(events, (;θ=u[:θ], ω=u[:ω], t=ctx.t, vidx=ctx.vidx, event_idx=i))
        end
    end
    vcb = VectorContinuousComponentCallback(cond2, affect2, 2)
    set_callback!(nw.im.vertexm[1], vcb)
    set_callback!(nw.im.vertexm[2], vcb)
    cbb = only(wrap_component_callbacks(nw))
    @test condition_outrange.(Ref(cbb), 1:2) == [1:2, 3:4]
    solve(ODEProblem(nw, u0, (0, 10.0), pflat(p0)), Tsit5())
    check_events(events)

    # the vector affect writes a state and a parameter
    nw = basenetwork()
    wcond = ComponentCondition([:ω]) do out, u, t
        out[1] = u[:ω] - 0.05
        out[2] = -u[:ω] - 0.05
    end
    waffect = ComponentAffect([:ω, :D]) do u, event_signs, ctx
        u[:ω] = 0
        u[:D] = 2
    end
    set_callback!(nw.im.vertexm[2], VectorContinuousComponentCallback(wcond, waffect, 2))
    s0 = NWState(nw)
    s0.p.v[2, :Pmech] = 2.5 # accelerate vertex 2
    sol = solve(ODEProblem(nw, s0, (0, 5)), Tsit5())
    D = sol[VPIndex(2, :D)]
    @test length(D) ≥ 2
    @test D[1] == 0.1 && all(==(2), D[2:end])
    tfire = sol.discretes[1].t[2]
    @test sol(tfire; idxs=VIndex(2, :ω), continuity=:left) ≈ 0.05
    @test sol(tfire; idxs=VIndex(2, :ω), continuity=:right) == 0

    # the batched vector affect is allocation free
    cbb = only(wrap_component_callbacks(nw))
    integ = init(ODEProblem(nw, s0, (0, 5)), Tsit5())
    vbatch = _batch_vector_affect(cbb)
    signs = Int8[1, 0]
    vbatch(integ, signs)
    @test (@b $vbatch($integ, $signs)).allocs == 0
end

@testset "single symbol list callbacks" begin
    # the line tripping scenario from above, with a lower limit on two lines
    function solve_tripping(cb)
        nw = basenetwork()
        set_callback!.(nw.im.edgem, Ref(cb))
        s0 = NWState(nw)
        s0.p.e[5, :limit] = 0.7
        s0.p.e[6, :limit] = 0.7
        s0.p.v[1, :Pmech] = 0.5
        nw, solve(ODEProblem(nw, s0, (0, 10)), Tsit5())
    end

    @testset "continuous callback reading observed, input and parameter" begin
        seen = []
        cond = ComponentCondition([:P, :limit]) do u, t
            abs(u[:P]) - u[:limit]
        end
        affect = ComponentAffect([:P, :active, :limit, :srcθ]) do u, ctx
            @test abs(u[:P]) ≈ u[:limit] atol=1e-6
            @test u[:srcθ] ≈ ctx.integrator[VIndex(ctx.src, :θ)]
            push!(seen, (ctx.t, ctx.eidx))
            u[:active] = 0
            @test u[:active] == 0 # the snapshot follows the write
            @test_throws ArgumentError u[:P] = 0.0    # output
            @test_throws ArgumentError u[:srcθ] = 0.0 # input
        end
        nw, sol = solve_tripping(ContinuousComponentCallback(cond, affect))

        # reference trips from the legacy form of the same callback
        @test last.(seen) == [5, 7, 6, 4, 3, 2, 1]
        @test first.(seen) ≈ [1.4839602700694718, 1.5682604084132876, 2.296184166113089, 2.3841068948960435,
                              2.438238941620517, 3.1356685492339857, 3.2051417467390775] atol=1e-4
        @test sol(10; idxs=EIndex(5, :P)) == 0

        cbb = only(wrap_component_callbacks(nw))
        batchcond = _batch_condition(cbb)
        out = fill(NaN, 7)
        s0 = NWState(nw)
        b = @b $batchcond($out, $(uflat(s0)), NaN, $((; p=pflat(s0))))
        @test b.allocs == 0
        @test out ≈ abs.(s0.e[1:7, :P]) .- s0.p.e[1:7, :limit]

        # gather and affect are allocation free, in the new and in the legacy form
        function affect_allocs(affect)
            nw = basenetwork()
            set_callback!.(nw.im.edgem, Ref(ContinuousComponentCallback(cond, affect)))
            cbb = only(wrap_component_callbacks(nw))
            integ = init(ODEProblem(nw, NWState(nw), (0, 1.0)), Tsit5())
            acc, batchaff = _scalar_member_affect(cbb, getaffect)
            scratch = _gather(acc, integ)
            batchaff(integ, scratch, 3) # afterwards the write changes nothing
            (@b _gather($acc, $integ)).allocs, (@b $batchaff($integ, $scratch, 3)).allocs
        end
        @test affect_allocs(ComponentAffect([:P, :active, :srcθ]) do u, ctx; u[:active] = 0 end) == (0, 0)
        @test affect_allocs(ComponentAffect([:P, :srcθ], [:active]) do u, p, ctx; p[:active] = 0 end) == (0, 0)
    end

    @testset "vertex affect: read input and observed, write state and parameter" begin
        nw = basenetwork()
        seen = []
        affect = ComponentAffect([:θ, :ω, :Pmech, :P, :Pdamping, :D]) do u, ctx
            push!(seen, (; P=u[:P], Pdamping=u[:Pdamping], ω=u[:ω], D=u[:D]))
            @test u[:P] ≈ ctx.integrator[VIndex(ctx.vidx, :P)]
            u[:ω] = 0.1
            u[:Pmech] = 0.0
        end
        set_callback!(nw.im.vertexm[2], PresetTimeComponentCallback(1.0, affect))
        s0 = NWState(nw)
        s0.v[:, :ω] .= 0.01 # get some damping
        sol = solve(ODEProblem(nw, s0, (0, 2)), Tsit5())

        e = only(seen)
        @test e.Pdamping ≈ -e.D * e.ω
        @test e.D == 0.1

        # both writes went into the integrator
        i1 = findlast(==(1.0), sol.t)
        @test sol[VIndex(2, :ω)][i1] == 0.1
        @test sol[VIndex(2, :ω)][i1-1] ≈ e.ω
        @test sol(1.5; idxs=VPIndex(2, :Pmech)) == 0.0
    end

    @testset "discrete callback with parameter in condition" begin
        nw = basenetwork()
        fired = []
        cond = ComponentCondition([:ω, :Pmech]) do u, t
            t > 0.5 && u[:Pmech] > 0
        end
        affect = ComponentAffect([:Pmech, :ω]) do u, ctx
            push!(fired, (ctx.t, ctx.vidx))
            u[:Pmech] = 0
            u[:ω] = 0.1
        end
        cb = DiscreteComponentCallback(cond, affect)
        set_callback!.(nw.im.vertexm, Ref(cb))
        @test length(wrap_component_callbacks(nw)) == 1
        sol = solve(ODEProblem(nw, NWState(nw), (0, 2)), Tsit5())

        # vertex 2 and 5 fire in the same step
        @test last.(fired) == [2, 5]
        @test fired[1][1] == fired[2][1] > 0.5
        @test sol[VPIndex(2, :Pmech)] == [1.5, 0]
        @test sol[VPIndex(5, :Pmech)] == [1.5, 0]
        i = findlast(==(fired[1][1]), sol.t)
        @test sol[VIndex(2, :ω)][i] == sol[VIndex(5, :ω)][i] == 0.1
    end

    @testset "continuous callback with affect_neg" begin
        f = (du, u, in, p, t) -> begin
            du[1] = -sin(t)
            du[2] = cos(t)
        end
        vm = VertexModel(; f, g=1:2, dim=2, indim=2, sym=[:cos=>1, :sin=>0], psym=[:ups=>0, :downs=>0])
        em = EdgeModel(; g=AntiSymmetric((out, in) -> out .= 0), outdim=2, indim=2)
        nw = Network(path_graph(2), vm, em; dealias=true)

        cond = ComponentCondition((u, t) -> u[:sin], [:sin])
        up = ComponentAffect([:ups]) do u, ctx; u[:ups] += 1 end
        down = ComponentAffect([:downs, :cos]) do u, ctx; u[:downs] += 1 end
        cb = ContinuousComponentCallback(cond, up; affect_neg! = down)
        s0 = NWState(nw)
        sol = solve(ODEProblem(nw, s0, (0, 4π+0.1); add_comp_cb=VIndex(1)=>cb), Tsit5())
        # sin starts at 0, so the crossings are down at π, 3π and up at 2π, 4π
        @test sol[VPIndex(1, :ups)] == [0, 0, 1, 1, 2]
        @test sol[VPIndex(1, :downs)] == [0, 1, 1, 2, 2]
        @test sol.discretes[1].t[2:end] ≈ [π, 2π, 3π, 4π] atol=1e-3
        # the saved parameters only record the two counters of vertex 1
        plog = sol.discretes[1].u
        @test plog isa NetworkDynamics.ParameterLog
        @test count(j -> isassigned(plog.tracks, j), eachindex(plog.tracks)) == 2
        s = NWState(sol, 3.5π)
        @test s.p.v[1, :ups] == 1 && s.p.v[1, :downs] == 2
    end

    @testset "batching" begin
        nw = basenetwork()
        # legacy conditions with the same function batch, even if built separately
        legacyf = (u, p, t) -> u[1] - p[1]
        affect = ComponentAffect((u, p, ctx) -> nothing, [], [:active])
        for e in nw.im.edgem
            set_callback!(e, ContinuousComponentCallback(ComponentCondition(legacyf, [:P], [:limit]), affect))
        end
        @test length(wrap_component_callbacks(nw)) == 1
        # ... also if the legacy symbols differ, each member resolves its own names
        set_callback!(nw.im.edgem[1], ContinuousComponentCallback(ComponentCondition(legacyf, [:₋P], [:limit]), affect))
        cbb = only(wrap_component_callbacks(nw))
        @test collect_c_or_a_indices(cbb, getcondition)[1:4] == [EIndex(1, :₋P), EIndex(1, :limit), EIndex(2, :P), EIndex(2, :limit)]
        s0 = NWState(nw)
        out = zeros(7)
        _batch_condition(cbb)(out, uflat(s0), 0.0, (; p=pflat(s0)))
        @test out[1] ≈ s0.e[1, :₋P] - s0.p.e[1, :limit]
        @test out[2] ≈ s0.e[2, :P] - s0.p.e[2, :limit]
        # different lengths still split
        set_callback!(nw.im.edgem[1], ContinuousComponentCallback(ComponentCondition(legacyf, [:₋P, :P], [:limit]), affect))
        @test length(wrap_component_callbacks(nw)) == 2

        # single list conditions with the same function and length batch, the symbols may differ
        f = (u, t) -> u[1] - u[2]
        aff = ComponentAffect((u, ctx) -> nothing, [:active])
        for (i, e) in pairs(nw.im.edgem)
            syms = isodd(i) ? [:P, :limit] : [:₋P, :K]
            set_callback!(e, ContinuousComponentCallback(ComponentCondition(f, syms), aff))
        end
        cbb = only(wrap_component_callbacks(nw))
        @test collect_c_or_a_indices(cbb, getcondition)[1:4] == [EIndex(1, :P), EIndex(1, :limit), EIndex(2, :₋P), EIndex(2, :K)]
        s0 = NWState(nw)
        out = zeros(7)
        _batch_condition(cbb)(out, uflat(s0), 0.0, (; p=pflat(s0)))
        @test out[1] ≈ s0.e[1, :P] - s0.p.e[1, :limit]
        @test out[2] ≈ s0.e[2, :₋P] - s0.p.e[2, :K]

        # closures from one place batch, each member keeps its own captures
        mkcond(c) = ComponentCondition((u, t) -> u[1] - c, [:P])
        for (i, e) in pairs(nw.im.edgem)
            set_callback!(e, ContinuousComponentCallback(mkcond(i), aff))
        end
        cbb = only(wrap_component_callbacks(nw))
        _batch_condition(cbb)(out, uflat(s0), 0.0, (; p=pflat(s0)))
        @test out ≈ [s0.e[i, :P] - i for i in 1:7]
    end

    @testset "closures capturing the namespace batch" begin
        # like a block in a component library: the callback is built per instance and captures the
        # namespaced symbols, so every instance has its own closure, all of the same type
        function limit_callback(ns, lim)
            x = Symbol(ns, :₊x)
            hit = Symbol(ns, :₊hit)
            cond = ComponentCondition([x, hit]) do u, t
                iszero(u[hit]) && u[x] > lim
            end
            affect = ComponentAffect([hit]) do u, ctx
                u[hit] = ctx.t
            end
            DiscreteComponentCallback(cond, affect)
        end
        function block(ns, lim)
            v = VertexModel(; f=(dx, x, ein, p, t) -> (dx[1] = 1.0; nothing), g=1,
                sym=[Symbol(ns, :₊x)=>0], psym=[Symbol(ns, :₊hit)=>0])
            set_callback!(v, limit_callback(ns, lim))
            v
        end
        nw = Network(SimpleGraph(2), [block(:a, 0.5), block(:b, 1.5)], EdgeModel[])
        cba = only(get_callbacks(nw.im.vertexm[1]))
        cbb = only(get_callbacks(nw.im.vertexm[2]))
        @test cba.condition.f !== cbb.condition.f
        @test NetworkDynamics._batchequal(cba, cbb)
        @test length(wrap_component_callbacks(nw)) == 1

        # each member still evaluates with its own captures
        sol = solve(ODEProblem(nw, NWState(nw), (0, 2)), Tsit5(); dtmax=0.1)
        ta = sol[VPIndex(1, :a₊hit)][end]
        tb = sol[VPIndex(2, :b₊hit)][end]
        @test 0.5 < ta < 0.61
        @test 1.5 < tb < 1.61
    end
end

@testset "wrong symboltype test" begin
    nw = basenetwork()
    nop = (args...) -> nothing

    # every named symbol is fine anywhere: param in condition sym, observed/input/output in
    # affect sym, state or observed in the legacy psym
    cond = ComponentCondition(nop, [:P, :₋P, :srcθ, :limit], [:limit, :K, :P])
    affect = ComponentAffect(nop, [:₋P, :srcθ],[:active, :Δθ])
    affect_neg! = ComponentAffect(nop, [:P], [:limit])
    cb = ContinuousComponentCallback(cond, affect; affect_neg!)
    set_callback!.(nw.im.edgem, Ref(cb))
    cbb = only(wrap_component_callbacks(nw))
    @test _batch_condition(cbb) isa Function
    @test _scalar_member_affect(cbb, getaffect) isa Tuple{AffectAccess, Function}
    @test _scalar_member_affect(cbb, getaffect_neg) isa Tuple{AffectAccess, Function}
    @test get_callbacks(nw) isa SciMLBase.DECallback

    cond = ComponentCondition(nop, [:θ, :ω, :P, :Pdamping, :M])
    affect = ComponentAffect(nop, [:θ, :P, :Pdamping, :Pmech])
    set_callback!(nw.im.vertexm[1], DiscreteComponentCallback(cond, affect))
    set_callback!(nw.im.vertexm[2], PresetTimeComponentCallback(1.0, affect))
    @test has_callback(nw.im.vertexm[1]) && has_callback(nw.im.vertexm[2])

    # undefined symbols throw on set_callback!
    good_cond = ComponentCondition(nop, [:P, :limit])
    good_aff = ComponentAffect(nop, [:active])
    bad_cond = ComponentCondition(nop, [:P, :doesnotexist])
    bad_aff = ComponentAffect(nop, [:active, :doesnotexist])
    e1 = nw.im.edgem[1]
    @test_throws ArgumentError set_callback!(e1, ContinuousComponentCallback(bad_cond, good_aff))
    @test_throws ArgumentError set_callback!(e1, ContinuousComponentCallback(good_cond, bad_aff))
    @test_throws ArgumentError set_callback!(e1, ContinuousComponentCallback(good_cond, good_aff; affect_neg! = bad_aff))
    # legacy form
    @test_throws ArgumentError set_callback!(e1, ContinuousComponentCallback(
        ComponentCondition(nop, [:P], [:doesnotexist]), ComponentAffect(nop, [], [:active])))
    # a symbol of a vertex is not a symbol of the edge
    @test_throws ArgumentError set_callback!(e1, PresetTimeComponentCallback(1.0, ComponentAffect(nop, [:ω])))

    # without the check, undefined symbols throw when the batch is built
    nw = basenetwork()
    cb = ContinuousComponentCallback(bad_cond, bad_aff; affect_neg! = bad_aff)
    set_callback!.(nw.im.edgem, Ref(cb); check=false)
    cbb = only(wrap_component_callbacks(nw))
    @test_throws ArgumentError _batch_condition(cbb)
    @test_throws ArgumentError _scalar_member_affect(cbb, getaffect)
    @test_throws ArgumentError _scalar_member_affect(cbb, getaffect_neg)
    @test_throws ArgumentError get_callbacks(nw)

    @test_throws ArgumentError AffectAccess(nw, [EIndex(1, :active), EIndex(1, :doesnotexist)])
    @test_throws ArgumentError AffectAccess(nw, [VIndex(1, :θ), VPIndex(1, :doesnotexist)])
    acc = AffectAccess(nw, [VIndex(1, :θ), VIndex(1, :Pmech), VPIndex(1, :D), VIndex(1, :P), VIndex(1, :Pdamping)])
    @test acc.uidx[1] > 0 && acc.pidx[1] == 0 # state
    @test acc.uidx[2] == 0 && acc.pidx[2] > 0 # param through VIndex
    @test acc.uidx[3] == 0 && acc.pidx[3] > 0 # param through VPIndex
    @test acc.uidx[4:5] == [0, 0] && acc.pidx[4:5] == [0, 0] # input and observed
end

@testset "check callbacks with different affneg" begin
    f = (du, u, in, p, t) -> begin
        du[1] = -sin(t)
        du[2] = cos(t)
    end
    vm = VertexModel(; f, g=1:2, dim=2, indim=2, sym=[:cos=>1, :sin=>0])

    ge = (out, in) -> begin
        out[1] = 0
        out[2] = 0
    end
    em = EdgeModel(; g=AntiSymmetric(ge), outdim=2, indim=2)

    g = path_graph(2)
    nw = Network(g, vm, em; dealias=true)

    # upcrossing no vector
    cond_vec = ComponentCondition([:sin, :cos], []) do out, u, p, t
        out[1] = u[:sin]
        out[2] = u[:cos]
    end
    vec_sin_up = []
    vec_cos_up = []
    vec_sin_down = []
    vec_cos_down = []
    # single sign-aware affect handles both crossing directions (up: +1, down: -1)
    affect_vec = ComponentAffect([], []) do u, p, event_signs, ctx
        pos = ctx.t/pi
        for i in eachindex(event_signs)
            s = event_signs[i]
            s == 0 && continue
            if i == 1
                println("sin $(s > 0 ? "up" : "down") at $pos")
                push!(s > 0 ? vec_sin_up : vec_sin_down, pos)
            else
                println("cos $(s > 0 ? "up" : "down") at $pos")
                push!(s > 0 ? vec_cos_up : vec_cos_down, pos)
            end
        end
    end
    cb_vec = VectorContinuousComponentCallback(cond_vec, affect_vec, 2)
    # set_callback!(nw[VIndex(1)], cb_vec)
    # set_callback!(nw[VIndex(1)], cb_vec)
    s0 = NWState(nw)
    prob = ODEProblem(nw, uflat(s0), (0, 4π+0.1), pflat(s0), add_comp_cb=VIndex(1)=>cb_vec)
    sol = solve(prob, Tsit5());
    @assert SciMLBase.successful_retcode(sol)

    @test vec_sin_up ≈ [2, 4] atol=1e-3
    @test vec_sin_down ≈ [1, 3] atol=1e-3
    @test vec_cos_up ≈ [1.5, 3.5] atol=1e-3
    @test vec_cos_down ≈ [0.5, 2.5] atol=1e-3

    ####
    #### only react to upcrossings (previously expressed via affect_neg! = nothing)
    ####
    empty!(vec_sin_up)
    empty!(vec_cos_up)
    empty!(vec_sin_down)
    empty!(vec_cos_down)
    affect_vec_uponly = ComponentAffect([], []) do u, p, event_signs, ctx
        pos = ctx.t/pi
        for i in eachindex(event_signs)
            event_signs[i] > 0 || continue  # skip no-event (0) and downcrossings (-1)
            i == 1 ? push!(vec_sin_up, pos) : push!(vec_cos_up, pos)
        end
    end
    cb_vec = VectorContinuousComponentCallback(cond_vec, affect_vec_uponly, 2)
    prob = ODEProblem(nw, uflat(s0), (0, 4π+0.1), pflat(s0), add_comp_cb=Dict(VIndex(1)=>cb_vec))
    sol = solve(prob, Tsit5());
    @assert SciMLBase.successful_retcode(sol)
    @test vec_sin_up ≈ [2, 4] atol=1e-3
    @test isempty(vec_sin_down)
    @test vec_cos_up ≈ [1.5, 3.5] atol=1e-3
    @test isempty(vec_cos_down)

    ####
    #### non-vector callback
    ####
    sin_up = []
    sin_down = []
    cond = ComponentCondition([:sin], []) do u, p ,t
        u[:sin]
    end
    affect = ComponentAffect([], []) do u, p, ctx
        pos = ctx.t/pi
        println("sin up at $pos")
        push!(sin_up, pos)
    end
    affect_neg = ComponentAffect([], []) do u, p, ctx
        pos = ctx.t/pi
        println("sin down at $pos")
        push!(sin_down, pos)
    end
    cb = ContinuousComponentCallback(cond, affect; affect_neg! = affect_neg)
    prob = ODEProblem(nw, uflat(s0), (0, 4π+0.1), pflat(s0), add_comp_cb=Dict(VIndex(1)=>cb))
    sol = solve(prob, Tsit5());
    @test sin_up ≈ [2, 4] atol=1e-3
    @test sin_down ≈ [1, 3] atol=1e-3

    empty!(sin_up)
    empty!(sin_down)
    cb = ContinuousComponentCallback(cond, affect; affect_neg! = nothing)
    prob = ODEProblem(nw, uflat(s0), (0, 4π+0.1), pflat(s0), add_comp_cb=Dict(VIndex(1)=>cb))
    sol = solve(prob, Tsit5());
    @test sin_up ≈ [2, 4] atol=1e-3
    @test isempty(sin_down)
end

@testset "callback constructors reject wrong arity" begin
    # a legacy (u, p, t) function passes the ComponentCondition arity check as f!(out, u, t),
    # so the callback constructors have to catch it
    f3 = (u, p, t) -> u[1] - p[1]
    a3 = (u, p, ctx) -> nothing
    c = ComponentCondition(f3, [:P, :limit])
    a = ComponentAffect(a3, [:active])
    @test_throws ArgumentError ContinuousComponentCallback(c, a)
    @test_throws ArgumentError DiscreteComponentCallback(c, a)
    @test_throws ArgumentError PresetTimeComponentCallback([1.0], a)
    # the other way round: scalar signatures in a vector callback
    cs = ComponentCondition((u, t) -> 0.0, [:P])
    as = ComponentAffect((u, ctx) -> nothing, [:P])
    @test_throws ArgumentError VectorContinuousComponentCallback(cs, as, 2)
    @test VectorContinuousComponentCallback(ComponentCondition((out, u, t) -> nothing, [:P]),
                                            ComponentAffect((u, signs, ctx) -> nothing, [:P]), 2) isa VectorContinuousComponentCallback
    # legacy forms define both arities and pass everywhere
    cl = ComponentCondition(f3, [:P], [:limit])
    al = ComponentAffect(a3, [], [:active])
    @test ContinuousComponentCallback(cl, al) isa ContinuousComponentCallback
    @test DiscreteComponentCallback(cl, al) isa DiscreteComponentCallback
    @test PresetTimeComponentCallback([1.0], al) isa PresetTimeComponentCallback
    @test VectorContinuousComponentCallback(cl, al, 2) isa VectorContinuousComponentCallback
end

@testset "ctx.dt_reset opt out" begin
    counter = ComponentAffect([:Pmech]) do u, ctx
        u[:Pmech] = u[:Pmech] + 1e-3 # bookkeeping change, keep the step
        ctx.dt_reset[] = false
    end
    kick = ComponentAffect([:ω]) do u, ctx
        u[:ω] = u[:ω] + 0.05
    end
    # an opted out affect alone keeps the step, a second affect asking for the reset wins
    nw1 = basenetwork()
    set_callback!(nw1.im.vertexm[1], PresetTimeComponentCallback([1.0], counter))
    sol1 = solve(ODEProblem(nw1, NWState(nw1), (0, 2.0)), Tsit5())
    nw2 = basenetwork()
    set_callback!(nw2.im.vertexm[1], PresetTimeComponentCallback([1.0], counter))
    set_callback!(nw2.im.vertexm[2], PresetTimeComponentCallback([1.0], kick))
    sol2 = solve(ODEProblem(nw2, NWState(nw2), (0, 2.0)), Tsit5())
    i1 = findlast(==(1.0), sol1.t); i2 = findlast(==(1.0), sol2.t) # the event time is saved twice
    step1 = sol1.t[i1+1] - 1.0
    step2 = sol2.t[i2+1] - 1.0
    @test step1 > 10*step2 # with the reset the first step after the event is tiny
    @test sol1[VPIndex(1, :Pmech)][end] ≈ -1 + 1e-3 # parameter change was still saved
    @test length(sol1[VPIndex(1, :Pmech)]) == 2

    # discrete batch: two members fire together, one opts out, the other one asks for the reset
    nw3 = basenetwork()
    fired = Int[]
    dc_out = DiscreteComponentCallback(ComponentCondition([:θ]) do u, t; t > 1.0 && u[:θ] < 1e3 end,
        ComponentAffect([:Pmech]) do u, ctx
            push!(fired, ctx.vidx)
            u[:Pmech] = 100.0 # stop firing again
            ctx.vidx == 1 && (ctx.dt_reset[] = false)
        end)
    set_callback!(nw3.im.vertexm[1], dc_out)
    set_callback!(nw3.im.vertexm[2], dc_out)
    @test length(wrap_component_callbacks(nw3)) == 1
    sol3 = solve(ODEProblem(nw3, NWState(nw3), (0, 2.0)), Tsit5())
    @test sort(fired) == [1, 2]
    it = findfirst(t -> t > 1.0, sol3.t)
    @test sol3.t[it+1] - sol3.t[it] < 1e-3 # vertex 2 forced the reset

    # continuous batch: both members cross at t=1, reset and parameter save happen once per event
    cc_cond = ComponentCondition([:θ]) do u, t; t - 1.0 end
    function solve_cc(optout)
        nw = basenetwork()
        aff = ComponentAffect([:Pmech]) do u, ctx
            u[:Pmech] = u[:Pmech] + 1e-3
            ctx.vidx in optout && (ctx.dt_reset[] = false)
        end
        cb = ContinuousComponentCallback(cc_cond, aff)
        set_callback!(nw.im.vertexm[1], cb)
        set_callback!(nw.im.vertexm[2], cb)
        @test length(wrap_component_callbacks(nw)) == 1
        solve(ODEProblem(nw, NWState(nw), (0, 2.0)), Tsit5())
    end
    sol_one = solve_cc([1])     # only vertex 2 asks for the reset
    sol_both = solve_cc([1, 2]) # nobody does
    i_one = findlast(==(1.0), sol_one.t); i_both = findlast(==(1.0), sol_both.t)
    # one member is enough for the reset, without it the step continues undisturbed
    @test sol_both.t[i_both+1] > sol_one.t[i_one+1]
    for sol in (sol_one, sol_both), v in 1:2
        ts = sol[VPIndex(v, :Pmech)]
        @test length(ts) == 2 # initial value and one save at the event
        @test ts[end] ≈ ts[1] + 1e-3
    end
end

@testset "ctx.derivative_discontinuity opt out" begin
    # Both vertices fire in one discrete batch at the tstop. The affect rewrites a parameter with
    # its own value, so the rhs stays the same.
    function solve_noop(optout)
        nw = basenetwork()
        cond = ComponentCondition((u, t) -> t == 1.0, [:Pmech])
        aff = ComponentAffect([:Pmech]) do u, ctx
            u[:Pmech] = u[:Pmech]
            ctx.vidx in optout && (ctx.derivative_discontinuity[] = false)
        end
        if !isnothing(optout)
            foreach(v -> set_callback!(nw.im.vertexm[v], DiscreteComponentCallback(cond, aff)), 1:2)
            @test length(wrap_component_callbacks(nw)) == 1
        end
        s0 = NWState(nw)
        s0.v[1, :ω] += 0.1
        solve(ODEProblem(nw, s0, (0, 2.0)), Tsit5(); tstops=[1.0])
    end
    sol_plain = solve_noop(nothing)
    sol_out = solve_noop([1, 2])
    sol_one = solve_noop([1])
    # without a discontinuity the solver continues as if the event was a plain tstop
    @test unique(sol_out.t) == sol_plain.t
    @test sol_out.stats.nf == sol_plain.stats.nf
    # one member reporting a discontinuity is enough
    @test sol_one.stats.nf > sol_plain.stats.nf
end

@testset "batched preset time callbacks" begin
    mkkick(dω) = ComponentAffect([:ω]) do u, ctx
        u[:ω] = u[:ω] + dω
    end
    function solve_with(cbs; kwargs...)
        nw = basenetwork()
        for (v, cb) in cbs
            add_callback!(nw.im.vertexm[v], cb)
        end
        nw, solve(ODEProblem(nw, NWState(nw), (0, 2.0)), Tsit5(); reltol=1e-10, abstol=1e-10, kwargs...)
    end

    # one closure per component on the same time grid ends up in one batch
    ts = [0.5, 1.0]
    nw, sol = solve_with([v => PresetTimeComponentCallback(ts, mkkick(0.01v)) for v in 1:5])
    @test length(wrap_component_callbacks(nw)) == 1
    @test get_callbacks(nw) isa DiscreteCallback # a single DiffEq callback

    # reference: the same kicks from a single network level callback
    nwref = basenetwork()
    ωidx = [NetworkDynamics.SII.variable_index(nwref, VIndex(v, :ω)) for v in 1:5]
    refkick = PresetTimeCallback(ts, integrator -> begin
        for v in 1:5
            integrator.u[ωidx[v]] += 0.01v
        end
        SciMLBase.auto_dt_reset!(integrator)
    end)
    solref = solve(ODEProblem(nwref, NWState(nwref), (0, 2.0); add_nw_cb=refkick), Tsit5();
        reltol=1e-10, abstol=1e-10)
    @test sol.u[end] ≈ solref.u[end] rtol=1e-8
    @test sol(0.75) ≈ solref(0.75) rtol=1e-8

    # different time grids form different batches
    nw, _ = solve_with([1 => PresetTimeComponentCallback(ts, mkkick(0.01)),
                        2 => PresetTimeComponentCallback(ts, mkkick(0.02)),
                        3 => PresetTimeComponentCallback([0.7], mkkick(0.03))])
    @test length(wrap_component_callbacks(nw)) == 2

    # different affect functions share a batch
    setp = ComponentAffect([:Pmech]) do u, ctx
        u[:Pmech] = 0.5
    end
    nw, sol = solve_with([1 => PresetTimeComponentCallback(ts, mkkick(0.1)),
                          2 => PresetTimeComponentCallback(ts, setp)])
    @test length(wrap_component_callbacks(nw)) == 1
    @test sol[VPIndex(2, :Pmech)][end] == 0.5
    iev = findall(==(0.5), sol.t) # saved before and after the event
    @test sol[VIndex(1, :ω)][last(iev)] ≈ sol[VIndex(1, :ω)][first(iev)] + 0.1

    # two callbacks on one component read the same snapshot, the last write wins
    _, sol = solve_with([1 => PresetTimeComponentCallback(ts, mkkick(0.1)),
                         1 => PresetTimeComponentCallback(ts, mkkick(0.2))])
    iev = findall(==(0.5), sol.t)
    @test sol[VIndex(1, :ω)][last(iev)] ≈ sol[VIndex(1, :ω)][first(iev)] + 0.2

    # the grid is compared like PresetTimeCallback stores it, the kwargs have to match
    function batchcount(cbs)
        nw = basenetwork()
        for (v, cb) in cbs
            add_callback!(nw.im.vertexm[v], cb)
        end
        length(wrap_component_callbacks(nw))
    end
    @test batchcount([1 => PresetTimeComponentCallback(1.0, setp),
                      2 => PresetTimeComponentCallback([1.0], setp)]) == 1
    @test batchcount([1 => PresetTimeComponentCallback(0.5:0.5:1.0, setp),
                      2 => PresetTimeComponentCallback([1.0, 0.5], setp)]) == 1
    @test batchcount([1 => PresetTimeComponentCallback(ts, setp),
                      2 => PresetTimeComponentCallback(ts, setp; save_positions=(false, false))]) == 2

    # vertex and edge members share a batch
    nw = basenetwork()
    deactivate = ComponentAffect([:active]) do u, ctx
        u[:active] = 0
    end
    add_callback!(nw.im.vertexm[1], PresetTimeComponentCallback(ts, setp))
    add_callback!(nw.im.edgem[1], PresetTimeComponentCallback(ts, deactivate))
    @test length(wrap_component_callbacks(nw)) == 1
    sol = solve(ODEProblem(nw, NWState(nw), (0, 2.0)), Tsit5())
    @test sol[VPIndex(1, :Pmech)][end] == 0.5
    @test sol[EPIndex(1, :active)][end] == 0

    # flags are combined over the members, one asking for the discontinuity is enough
    function noop(optout)
        ComponentAffect([:Pmech]) do u, ctx
            u[:Pmech] = u[:Pmech]
            optout && (ctx.derivative_discontinuity[] = false)
        end
    end
    _, sol_plain = solve_with([]; tstops=[1.0])
    nw, sol_out = solve_with([v => PresetTimeComponentCallback([1.0], noop(true)) for v in 1:2])
    @test length(wrap_component_callbacks(nw)) == 1
    @test sol_out.stats.nf == sol_plain.stats.nf
    _, sol_one = solve_with([1 => PresetTimeComponentCallback([1.0], noop(true)),
                             2 => PresetTimeComponentCallback([1.0], noop(false))])
    @test sol_one.stats.nf > sol_plain.stats.nf
end

@testset "symbolic view test" begin
    a = collect(1:10)
    v = SymbolicView(view(a,1:3), (:a,:b,:c))
    @test v[:a] == 1
    @test v[:b] == 2
    @test v[:c] == 3
    v[:c] = 7
    @test a == [1,2,7,4,5,6,7,8,9,10]
end

@testset "ODEProblem callback keywords" begin
    nw = basenetwork()

    # Create reusable callbacks
    # Component callback (similar to those used in batch tests)
    triggered_comp = Int[]
    cond = ComponentCondition([:P, :₋P, :srcθ], [:limit, :K]) do u, p, t
        t - 1 # trigger at t=1
    end
    affect = ComponentAffect([],[:active]) do u, p, ctx
        push!(triggered_comp, ctx.eidx)
        p[:active] = 0
    end
    comp_cb = ContinuousComponentCallback(cond, affect)

    # Network-level callback (PresetTimeCallback)
    triggered_nw = Float64[]
    nw_cb = DiffEqCallbacks.PresetTimeCallback([1.0, 2.0], integrator -> push!(triggered_nw, integrator.t))

    s0 = NWState(nw)
    tspan = (0.0, 3.0)

    @testset "add_comp_cb: add component callbacks" begin
        empty!(triggered_comp)
        # Test adding component callback via keyword
        prob = ODEProblem(nw, s0, tspan; add_comp_cb=Dict(EIndex(3)=>comp_cb))
        # Verify callback was added (by checking it's in the callback structure)
        @test !isnothing(prob.kwargs[:callback])

        # Solve to verify the component callback actually triggers
        sol = solve(prob, Tsit5())
        @test 3 in triggered_comp  # EIndex(3) should trigger
    end

    @testset "add_nw_cb: add network-level callback" begin
        empty!(triggered_nw)
        # Test adding network-level callback via keyword
        prob = ODEProblem(nw, s0, tspan; add_nw_cb=nw_cb)
        @test prob isa ODEProblem
        @test !isnothing(prob.kwargs[:callback])

        # Solve to verify the network callback actually triggers
        sol = solve(prob, Tsit5())
        @test length(triggered_nw) == 2
        @test triggered_nw ≈ [1.0, 2.0] atol=1e-10
    end

    @testset "combined: add_comp_cb + add_nw_cb" begin
        empty!(triggered_comp)
        empty!(triggered_nw)

        # Test both keywords together
        prob = ODEProblem(nw, s0, tspan;
                         add_comp_cb=Dict(EIndex(3)=>comp_cb),
                         add_nw_cb=nw_cb)
        @test prob isa ODEProblem
        @test !isnothing(prob.kwargs[:callback])

        # Solve to verify both callbacks work
        sol = solve(prob, Tsit5())
        @test length(triggered_nw) == 2
        @test triggered_nw ≈ [1.0, 2.0] atol=1e-10
    end

    @testset "override_cb: completely replace callbacks" begin
        # Create a simple override callback
        override_triggered = Float64[]
        override_cb = DiffEqCallbacks.PresetTimeCallback([0.5], integrator -> push!(override_triggered, integrator.t))

        # Test override completely replaces network callbacks
        prob = ODEProblem(nw, s0, tspan; override_cb=override_cb)
        @test prob isa ODEProblem
        @test !isnothing(prob.kwargs[:callback])
        @test prob.kwargs[:callback] isa DiscreteCallback

        # Solve and verify only override callback triggered
        sol = solve(prob, Tsit5())
        @test length(override_triggered) == 1
        @test override_triggered[1] ≈ 0.5 atol=1e-10
    end

    @testset "error cases" begin
        # Test that passing callback directly throws error
        @test_throws ArgumentError ODEProblem(nw, s0, tspan; callback=nw_cb)
        override_cb = DiffEqCallbacks.PresetTimeCallback([0.5], integrator -> nothing)
        # Test that combining override_cb with add_comp_cb throws error
        @test_throws ArgumentError ODEProblem(nw, s0, tspan;
                                             override_cb=override_cb,
                                             add_comp_cb=Dict(EIndex(3)=>comp_cb))
        # Test that combining override_cb with add_nw_cb throws error
        @test_throws ArgumentError ODEProblem(nw, s0, tspan;
                                             override_cb=override_cb,
                                             add_nw_cb=nw_cb)
    end

    # Create a network with embedded callbacks for testing callback merging
    nw_with_cb = basenetwork()
    embedded_triggered = Ref{Float64}(0.0)
    embedded_cond = ComponentCondition([:P, :₋P, :srcθ], [:limit, :K]) do u, p, t
        t > 0.5 && iszero(embedded_triggered[])
    end
    embedded_affect = ComponentAffect([], [:limit]) do u, p, ctx
        embedded_triggered[] = ctx.t
    end
    embedded_cb = DiscreteComponentCallback(embedded_cond, embedded_affect)
    add_callback!(nw_with_cb[EIndex(1)], embedded_cb)
    s0_with_cb = NWState(nw_with_cb)

    @testset "combined with existing network callbacks" begin
        embedded_triggered[] = 0.0
        empty!(triggered_comp)
        empty!(triggered_nw)

        # Network has embedded callback, add both component and network callbacks
        prob = ODEProblem(nw_with_cb, s0_with_cb, tspan;
                         add_comp_cb=Dict(EIndex(3)=>comp_cb),
                         add_nw_cb=nw_cb)
        @test prob isa ODEProblem
        @test !isnothing(prob.kwargs[:callback])

        # Solve and verify all three callback types trigger
        sol = solve(prob, Tsit5())
        @test embedded_triggered[] > 0.5  # Embedded callback triggered
        @test !isempty(triggered_comp)  # Additional component callback triggered
        @test 3 in triggered_comp
        @test length(triggered_nw) == 2  # Network callback triggered
        @test triggered_nw ≈ [1.0, 2.0] atol=1e-10
    end

    @testset "specialize keyword and network type" begin
        # the network specific types (component functions) live in the core
        corestr = string(typeof(nw_with_cb.core))
        prob_full = ODEProblem(nw_with_cb, s0_with_cb, tspan)
        prob_auto = ODEProblem(nw_with_cb, s0_with_cb, tspan; specialize=SciMLBase.AutoSpecialize)
        @test SciMLBase.specialization(prob_full.f) == SciMLBase.FullSpecialize
        @test SciMLBase.specialization(prob_auto.f) == SciMLBase.AutoSpecialize
        @test prob_full.f.mass_matrix == prob_auto.f.mass_matrix

        # symbolic indexing goes through the network itself
        @test prob_full.f.sys === prob_auto.f.sys === nw_with_cb
        @test extract_nw(prob_full) === extract_nw(prob_auto) === nw_with_cb

        # the core type stays out of the callbacks and the integrator
        @test !occursin(corestr, string(typeof(prob_full.kwargs[:callback])))
        integ_full = init(prob_full, Tsit5())
        integ_auto = init(prob_auto, Tsit5())
        @test !occursin(corestr, string(typeof(integ_full)))
        @test !occursin(corestr, string(typeof(integ_auto)))
        @test extract_nw(integ_auto) === nw_with_cb

        # a fully typed network puts it back, so the check above can fail
        nw_typed = Network(nw_with_cb; fullytyped=true)
        @test string(typeof(nw_typed.core)) == corestr
        @test occursin(corestr, string(typeof(init(ODEProblem(nw_typed, uflat(s0_with_cb), tspan, pflat(s0_with_cb)), Tsit5()))))

        embedded_triggered[] = 0.0
        sol_full = solve(prob_full, Tsit5())
        @test embedded_triggered[] > 0
        embedded_triggered[] = 0.0
        sol_auto = solve(prob_auto, Tsit5())
        @test embedded_triggered[] > 0
        @test sol_full.t == sol_auto.t
        @test sol_full.u == sol_auto.u
        @test extract_nw(sol_auto) === nw_with_cb
        @test sol_auto(1.0; idxs=VIndex(1, :θ)) == sol_full(1.0; idxs=VIndex(1, :θ))
    end
end

@testset "iterative discrete callbacks" begin
    # a ramp x' = 1 from x(0) = 0, the parameters are the memory of the discrete blocks
    function ramp_vertex(; psym=[:flag=>0, :held=>-1, :above=>0])
        VertexModel(; f=(dx, x, ein, p, t) -> (dx[1] = 1.0; nothing), g=1, sym=[:x=>0], psym)
    end
    function one_vertex_nw(v)
        Network(SimpleGraph(1), [v], EdgeModel[])
    end
    isolve(nw; alg=Tsit5(), kwargs...) = solve(ODEProblem(nw, NWState(nw), (0, 2); kwargs...), alg; dtmax=0.1)

    @testset "order independence" begin
        # a hysteresis switches `flag` when x crosses 1, a sample-and-hold stores `flag` on the same crossing
        blocks(iterative) = (DiscreteComponentCallback(
            ComponentCondition([:x, :flag]) do u, t
                iszero(u[:flag]) ? u[:x] > 1 : u[:x] < 1
            end,
            ComponentAffect([:flag]) do u, ctx
                u[:flag] = 1 - u[:flag]
            end; iterative),
        DiscreteComponentCallback(
            ComponentCondition([:x, :above]) do u, t
                iszero(u[:above]) ? u[:x] > 1 : u[:x] < 1
            end,
            ComponentAffect([:held, :above, :flag]) do u, ctx
                u[:held] = u[:flag]
                u[:above] = 1 - u[:above]
            end; iterative))
        function held(iterative, hystfirst)
            hyst, hold = blocks(iterative)
            v = ramp_vertex()
            set_callback!(v, hystfirst ? (hyst, hold) : (hold, hyst))
            nw = one_vertex_nw(v)
            cbb = wrap_component_callbacks(nw)
            iterative && @test only(cbb) isa IterativeBatches && length(only(cbb).batches) == 2
            sol = isolve(nw)
            @test sol[VPIndex(1, :flag)][end] == 1
            sol[VPIndex(1, :held)][end]
        end
        # synchronous: the hold stores the flag from before the switch, whatever the order
        @test held(true, true) == held(true, false) == 0
        # as plain discrete callbacks the batches run one after the other
        @test held(false, true) == 1
        @test held(false, false) == 0
    end

    @testset "cascade within one instant" begin
        # a block sets the limit to zero, a limited state above the limit is put onto it at the same t
        fired = []
        block = DiscreteComponentCallback(
            ComponentCondition([:x, :blocked]) do u, t
                iszero(u[:blocked]) && u[:x] > 1
            end,
            ComponentAffect([:blocked, :lim]) do u, ctx
                push!(fired, (:block, ctx.t))
                u[:blocked] = 1
                u[:lim] = 0.5
            end; iterative=true)
        limit = DiscreteComponentCallback(
            ComponentCondition([:x, :lim]) do u, t
                u[:x] > u[:lim]
            end,
            ComponentAffect([:x, :lim]) do u, ctx
                push!(fired, (:limit, ctx.t))
                u[:x] = u[:lim]
            end; iterative=true)
        for cbs in ((block, limit), (limit, block))
            empty!(fired)
            # x stops at its limit, p[2] = lim
            f = (dx, x, ein, p, t) -> (dx[1] = x[1] < p[2] ? 1.0 : 0.0; nothing)
            v = VertexModel(; f, g=1, sym=[:x=>0], psym=[:blocked=>0, :lim=>10])
            set_callback!(v, cbs)
            nw = one_vertex_nw(v)
            sol = isolve(nw)
            @test first.(fired) == [:block, :limit]
            @test fired[1][2] == fired[2][2]
            it = findlast(==(fired[1][2]), sol.t)
            @test sol[VIndex(1, :x)][it] == 0.5
            @test length(sol[VPIndex(1, :lim)]) == 2 # one parameter save per instant
        end
    end

    @testset "jump from a network level callback" begin
        # a plain discrete callback steps `inp` in the first step after t=0.5, it is sorted behind the
        # component callbacks. The iterative block still reacts at the same instant.
        function run(iterative; preset=false)
            fired = Float64[]
            cb = DiscreteComponentCallback(
                ComponentCondition([:inp, :flag]) do u, t
                    iszero(u[:flag]) && u[:inp] > 0.5
                end,
                ComponentAffect([:flag]) do u, ctx
                    push!(fired, ctx.t)
                    u[:flag] = 1
                    ctx.dt_reset[] = false # bookkeeping only, avoids a dt reset right at the tstop
                end; iterative)
            v = ramp_vertex(psym=[:inp=>0, :flag=>0])
            set_callback!(v, cb)
            nw = one_vertex_nw(v)
            pidx = NetworkDynamics.SII.parameter_index(nw, VPIndex(1, :inp))
            tstep = Ref(NaN)
            stepaffect = integrator -> begin
                tstep[] = integrator.t
                integrator.p[pidx] = 1
                save_parameters!(integrator)
            end
            step = if preset
                PresetTimeCallback(0.5, stepaffect)
            else
                DiscreteCallback((u, t, integrator) -> t ≥ 0.5 && isnan(tstep[]), stepaffect)
            end
            isolve(nw; add_nw_cb=step)
            only(fired), tstep[]
        end
        fired, tstep = run(true)
        @test fired == tstep
        fired, tstep = run(false)
        @test fired > tstep # runs before the network level callback, sees the jump one step late
        # a preset-time callback from `add_nw_cb` is sorted to the front, so everybody sees it
        @test run(true; preset=true) == (0.5, 0.5)
        @test run(false; preset=true) == (0.5, 0.5)
    end

    @testset "preset-time component callbacks run first" begin
        fired = Float64[]
        cb = DiscreteComponentCallback(
            ComponentCondition([:inp, :flag]) do u, t
                iszero(u[:flag]) && u[:inp] > 0.5
            end,
            ComponentAffect([:flag]) do u, ctx
                push!(fired, ctx.t)
                u[:flag] = 1
                ctx.dt_reset[] = false
            end)
        step = PresetTimeComponentCallback(0.5, ComponentAffect([:inp]) do u, ctx
            u[:inp] = 1
            ctx.dt_reset[] = false
        end)
        v = ramp_vertex(psym=[:inp=>0, :flag=>0])
        set_callback!(v, (cb, step)) # registered after the discrete callback
        nw = one_vertex_nw(v)
        isolve(nw)
        @test only(fired) == 0.5
    end

    @testset "DAE reinit between rounds" begin
        # 0 = k⋅x - y, the first block flips k, the second one reacts to the new algebraic y
        fired = []
        flip = DiscreteComponentCallback(
            ComponentCondition([:x, :flag1]) do u, t
                iszero(u[:flag1]) && u[:x] > 1
            end,
            ComponentAffect([:k, :flag1]) do u, ctx
                push!(fired, (:flip, ctx.t))
                u[:k] = -1
                u[:flag1] = 1
            end; iterative=true)
        react = DiscreteComponentCallback(
            ComponentCondition([:y, :flag2]) do u, t
                iszero(u[:flag2]) && u[:y] < 0
            end,
            ComponentAffect([:flag2]) do u, ctx
                push!(fired, (:react, ctx.t))
                u[:flag2] = 1
            end; iterative=true)
        f = (du, u, ein, p, t) -> begin
            du[1] = 1.0
            du[2] = p[1] * u[1] - u[2]
            nothing
        end
        v = VertexModel(; f, g=1, sym=[:x=>0, :y=>0], psym=[:k=>1, :flag1=>0, :flag2=>0],
            mass_matrix=Diagonal([1, 0]))
        set_callback!(v, (flip, react))
        nw = one_vertex_nw(v)
        sol = isolve(nw; alg=Rodas5P())
        @test SciMLBase.successful_retcode(sol)
        @test first.(fired) == [:flip, :react]
        @test fired[1][2] == fired[2][2]
        it = findlast(==(fired[1][2]), sol.t)
        @test sol[VIndex(1, :y)][it] ≈ -sol[VIndex(1, :x)][it]
    end

    @testset "bounded iteration" begin
        toggle = DiscreteComponentCallback(
            ComponentCondition([:x]) do u, t
                u[:x] > 1
            end,
            ComponentAffect([:flag]) do u, ctx
                u[:flag] = 1 - u[:flag]
            end; iterative=true)
        v = ramp_vertex(psym=[:flag=>0])
        set_callback!(v, toggle)
        nw = one_vertex_nw(v)
        @test_logs (:warn, r"after 3 of at most 3 rounds.*VIndex\(1\)") match_mode=:any isolve(nw; event_maxiter=3)
        @test_throws ErrorException isolve(nw; event_failure=:error)
        @test_throws ArgumentError isolve(nw; event_failure=:foo)

        # an affect that writes nothing would loop forever, it is reported after the first round
        noop = DiscreteComponentCallback(
            ComponentCondition([:x]) do u, t
                u[:x] > 1
            end,
            ComponentAffect([:flag]) do u, ctx
                u[:flag] = 0
            end; iterative=true)
        set_callback!(v, noop)
        nw = one_vertex_nw(v)
        @test_logs (:warn, r"after 1 of at most 10 rounds") match_mode=:any isolve(nw)
    end

    @testset "assembly" begin
        c = ComponentCondition([:x]) do u, t; u[:x] > 1 end
        a = ComponentAffect([:flag]) do u, ctx; u[:flag] = 1 end
        @test_throws ArgumentError DiscreteComponentCallback(c, a; iterative=true, save_positions=(false, false))
        @test_throws ArgumentError DiscreteComponentCallback(c, a; iterative=true, initializealg=SciMLBase.NoInit())
        @test contains(repr(MIME"text/plain"(), DiscreteComponentCallback(c, a; iterative=true)), "iterative")

        # closures from one place share a batch, also with different captures
        mkcond(lim) = ComponentCondition([:x]) do u, t; u[:x] > lim end
        v1 = ramp_vertex(psym=[:flag=>0]); v2 = ramp_vertex(psym=[:flag=>0])
        set_callback!(v1, DiscreteComponentCallback(mkcond(1), a; iterative=true))
        set_callback!(v2, DiscreteComponentCallback(mkcond(2), a; iterative=true))
        nw = Network(SimpleGraph(2), [v1, v2], EdgeModel[])
        @test length(only(wrap_component_callbacks(nw)).batches) == 1
        set_callback!(v1, DiscreteComponentCallback(mkcond(1), a))
        set_callback!(v2, DiscreteComponentCallback(mkcond(2), a))
        nw = Network(SimpleGraph(2), [v1, v2], EdgeModel[])
        @test length(wrap_component_callbacks(nw)) == 1

        # order at an instant: preset first, then other discrete, then add_nw_cb, the iterative set last
        v = ramp_vertex(psym=[:flag=>0, :inp=>0])
        disc = DiscreteComponentCallback(ComponentCondition([:x]) do u, t; false end, a)
        iter = DiscreteComponentCallback(c, a; iterative=true)
        cont = ContinuousComponentCallback(ComponentCondition([:x]) do u, t; u[:x] - 1 end, a)
        preset = PresetTimeComponentCallback(0.5, ComponentAffect([:inp]) do u, ctx; u[:inp] = 1 end)
        set_callback!(v, (iter, disc, cont, preset))
        nw = one_vertex_nw(v)
        cbs = get_callbacks(nw)
        @test length(cbs.continuous_callbacks) == 1
        dcbs = cbs.discrete_callbacks
        @test length(dcbs) == 3
        @test dcbs[1].condition isa DiffEqCallbacks.PresetTimeFunction
        @test dcbs[end].initializealg isa SciMLBase.NoInit
        user = DiscreteCallback((u, t, integrator) -> false, integrator -> nothing)
        userpreset = PresetTimeCallback(0.7, integrator -> nothing)
        prob = ODEProblem(nw, NWState(nw), (0, 1); add_nw_cb=CallbackSet(user, userpreset))
        dcbs = prob.kwargs[:callback].discrete_callbacks
        @test length(dcbs) == 5
        @test dcbs[1].condition isa DiffEqCallbacks.PresetTimeFunction # the component one
        @test dcbs[2] === userpreset
        @test dcbs[4] === user
        @test dcbs[5].initializealg isa SciMLBase.NoInit
    end

    @testset "quiet step does not allocate" begin
        a = ComponentAffect([:flag]) do u, ctx; u[:flag] = 1 end
        mkcond(lim) = ComponentCondition([:x]) do u, t; u[:x] > lim end
        v = ramp_vertex(psym=[:flag=>0])
        set_callback!(v, (DiscreteComponentCallback(mkcond(1), a; iterative=true),
                          DiscreteComponentCallback(ComponentCondition([:x]) do u, t; u[:x] < -1 end, a; iterative=true)))
        nw = one_vertex_nw(v)
        cb = get_callbacks(nw)
        integrator = SciMLBase.init(ODEProblem(nw, NWState(nw), (0, 1)), Tsit5())
        function count_allocs(cond, integrator)
            cond(integrator.u, integrator.t, integrator)
            @allocated cond(integrator.u, integrator.t, integrator)
        end
        @test !cb.condition(integrator.u, integrator.t, integrator)
        @test count_allocs(cb.condition, integrator) == 0
    end
end
