using NetworkDynamics
using Graphs
using OrdinaryDiffEqTsit5
using Chairmarks
using Test
using ModelingToolkit
using DiffEqCallbacks

(isinteractive() && @__MODULE__()==Main ? includet : include)("ComponentLibrary.jl")

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

    batches = NetworkDynamics.wrap_component_callbacks(nw);
    @test length(batches) == 1
    cbb = only(batches);

    @test NetworkDynamics.condition_dim(cbb) == 3
    @test NetworkDynamics.condition_pdim(cbb) == 2
    @test all(NetworkDynamics.affect_dim.(Ref(cbb), NetworkDynamics.getaffect, 1:7) .== 0)
    @test all(NetworkDynamics.affect_pdim.(Ref(cbb), NetworkDynamics.getaffect, 1:7) .== 1)

    @test NetworkDynamics.condition_urange.(Ref(cbb), 1:length(cbb)) == [1:3,4:6,7:9,10:12,13:15,16:18,19:21]
    @test NetworkDynamics.condition_prange.(Ref(cbb), 1:length(cbb)) == [1:2,3:4,5:6,7:8,9:10,11:12,13:14]
    @test NetworkDynamics.affect_urange.(Ref(cbb), NetworkDynamics.getaffect, 1:length(cbb)) == [1:0 for i in 1:7]
    @test NetworkDynamics.affect_prange.(Ref(cbb), NetworkDynamics.getaffect, 1:length(cbb)) == [1:1,2:2,3:3,4:4,5:5,6:6,7:7]
    @test NetworkDynamics.condition_outrange.(Ref(cbb), 1:length(cbb)) == [1:1,2:2,3:3,4:4,5:5,6:6,7:7]

    @test NetworkDynamics.collect_c_or_a_indices(cbb, NetworkDynamics.getcondition, :sym) == collect(Iterators.flatten(collect(EIndex(i, [:P, :₋P, :srcθ])) for i in 1:7))
    @test NetworkDynamics.collect_c_or_a_indices(cbb, NetworkDynamics.getcondition, :psym) == collect(Iterators.flatten(collect(EPIndex(i, [:limit, :K])) for i in 1:7))
    @test NetworkDynamics.collect_c_or_a_indices(cbb, NetworkDynamics.getaffect, :sym) == []
    @test NetworkDynamics.collect_c_or_a_indices(cbb, NetworkDynamics.getaffect, :psym) == collect(EPIndex(i, :active) for i in 1:7)

    batchcond = NetworkDynamics._batch_condition(cbb)
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

    nwcb = NetworkDynamics.get_callbacks(nw);
    s0 = NWState(nw)
    prob = ODEProblem(nw, uflat(s0), (0,6), copy(pflat(s0)), callback=nwcb)
    sol = solve(prob, Tsit5());

    @test 0.1 < useless_triggertime[] <= 1.0

    @test tripi == [5,7,4,1,3,2]
    tref = [1, 2.247676397005474, 2.502523192233235, 3.1947647115093654, 3.3380530127462587, 3.4042696241577888]
    @test maximum(abs.(tript - tref)) < 1e-5
end

@testset "show functions for callbacks" begin
    nw = basenetwork()
    v = nw.im.vertexm[1]

    empty_function = (args...) -> nothing
    cond = ComponentCondition(empty_function, [:θ, :ω], [])
    affect = ComponentAffect(empty_function, [:θ, :ω],[])
    cb = VectorContinuousComponentCallback(cond, affect, 2)
    add_callback!(v, cb)
    show(stdout, MIME"text/plain"(), v)
    cond = ComponentCondition(empty_function, [:θ], [])
    affect = ComponentAffect(empty_function, [],[:M])
    cb2 = ContinuousComponentCallback(cond, affect)
    add_callback!(v, cb2)
    show(stdout, MIME"text/plain"(), v)
    show(stdout, MIME"text/plain"(), nw)

    cond = ComponentCondition(empty_function, [:θ], [])
    affect = ComponentAffect(empty_function, [:ω],[])
    affect_neg = ComponentAffect(empty_function, [:θ], [:M])
    cb3 = ContinuousComponentCallback(cond, affect; affect_neg! = affect_neg)
    show(stdout, MIME"text/plain"(), cb3)

    cond = ComponentCondition(empty_function, [:θ], [])
    affect = ComponentAffect(empty_function, [:ω],[])
    cb3 = ContinuousComponentCallback(cond, affect; affect_neg! = nothing)
    show(stdout, MIME"text/plain"(), cb3)
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
    affect = ComponentAffect([:θ, :ω],[]) do u, p, event_idx, ctx
        push!(events, (;θ=u[:θ], ω=u[:ω], t=ctx.t, vidx=ctx.vidx, event_idx=event_idx))
        @info "Triggered event_idx $event_idx t=$(ctx.t) on $(ctx.vidx)"
    end
    ccb = VectorContinuousComponentCallback(cond, affect, 2)
    set_callback!(nw.im.vertexm[1], ccb)
    set_callback!(nw.im.vertexm[2], ccb)

    cbbs = NetworkDynamics.wrap_component_callbacks(nw);
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
    for (e, re) in zip(events, ref_events)
        @test e.vidx == re.vidx
        @test e.event_idx == re.event_idx
        @test abs(e.θ - re.θ) < 1e-5
        @test abs(e.ω - re.ω) < 1e-5
        @test abs(e.t - re.t) < 1e-5
    end
end

@testset "wrong symboltype test" begin
    nw = basenetwork()

    # invalid pram in condition u
    # invalid obs in affect u
    cond = ComponentCondition((args...)->nothing, [:P, :₋P, :srcθ, :limit], [:limit, :K])
    affect = ComponentAffect((args...)->nothing, [:₋P],[:active])
    cb = ContinuousComponentCallback(cond, affect)
    set_callback!.(nw.im.edgem, Ref(cb); check=false);
    cbb = only(NetworkDynamics.wrap_component_callbacks(nw));
    # oserved can handle parameters now!
    # @test_throws ArgumentError NetworkDynamics._batch_condition(cbb)
    @test_throws ArgumentError NetworkDynamics._batch_affect(cbb, NetworkDynamics.getaffect)

    # invalid state in condition p
    cond = ComponentCondition((args...)->nothing, [:P, :₋P, :srcθ], [:limit, :K, :P])
    affect = ComponentAffect((args...)->nothing, [],[:active, :₋P])
    cb = ContinuousComponentCallback(cond, affect)
    set_callback!.(nw.im.edgem, Ref(cb); check=false);
    cbb = only(NetworkDynamics.wrap_component_callbacks(nw));
    @test_throws ArgumentError NetworkDynamics._batch_condition(cbb)
    @test_throws ArgumentError NetworkDynamics._batch_affect(cbb, NetworkDynamics.getaffect)

    # test on set_callback
    cond = ComponentCondition((args...)->nothing, [:P, :₋P, :srcθ, :limit], [:limit, :K])
    affect = ComponentAffect((args...)->nothing, [],[:active])
    cb = ContinuousComponentCallback(cond, affect)
    @test_throws ArgumentError set_callback!(nw.im.edgem[1], cb)
    cond = ComponentCondition((args...)->nothing, [:P, :₋P, :srcθ], [:limit, :K, :P])
    affect = ComponentAffect((args...)->nothing, [],[:active])
    cb = ContinuousComponentCallback(cond, affect)
    @test_throws ArgumentError set_callback!(nw.im.edgem[1], cb)
    cond = ComponentCondition((args...)->nothing, [:P, :₋P, :srcθ], [:limit, :K])
    affect = ComponentAffect((args...)->nothing, [:₋P],[:active])
    cb = ContinuousComponentCallback(cond, affect)
    @test_throws ArgumentError set_callback!(nw.im.edgem[1], cb)
    cond = ComponentCondition((args...)->nothing, [:P, :₋P, :srcθ], [:limit, :K])
    affect = ComponentAffect((args...)->nothing, [],[:active, :P])
    cb = ContinuousComponentCallback(cond, affect)
    @test_throws ArgumentError set_callback!(nw.im.edgem[1], cb)

    # test conditions for negative affect
    affect  = ComponentAffect((args...)->nothing, [],[])
    cond = ComponentCondition((args...)->nothing, [:P, :₋P, :srcθ], [:limit, :K])
    affect_neg! = ComponentAffect((args...)->nothing, [:₋P],[:active])
    cb = ContinuousComponentCallback(cond, affect; affect_neg!)
    @test_throws ArgumentError set_callback!(nw.im.edgem[1], cb)
    cond = ComponentCondition((args...)->nothing, [:P, :₋P, :srcθ], [:limit, :K])
    affect_neg! = ComponentAffect((args...)->nothing, [],[:active, :P])
    cb = ContinuousComponentCallback(cond, affect; affect_neg!)
    @test_throws ArgumentError set_callback!(nw.im.edgem[1], cb)
end

# @testset "check callbacks with different affneg" begin
begin
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
    affect_vec = ComponentAffect([], []) do u, p, event_idx, ctx
        pos = ctx.t/pi
        if event_idx == 1
            println("sin up at $pos")
            push!(vec_sin_up, pos)
        else
            println("cos up at $pos")
            push!(vec_cos_up, pos)
        end
    end
    affect_vec_neg = ComponentAffect([], []) do u, p, event_idx, ctx
        pos = ctx.t/pi
        if event_idx == 1
            println("sin down at $pos")
            push!(vec_sin_down, pos)
        else
            println("cos down at $pos")
            push!(vec_cos_down, pos)
        end
    end
    cb_vec = VectorContinuousComponentCallback(cond_vec, affect_vec, 2; affect_neg! = affect_vec_neg)
    # set_callback!(nw[VIndex(1)], cb_vec)
    # set_callback!(nw[VIndex(1)], cb_vec)
    s0 = NWState(nw)
    prob = ODEProblem(nw, uflat(s0), (0, 4π+0.1), pflat(s0), callback=get_callbacks(nw, VIndex(1)=>cb_vec))
    sol = solve(prob, Tsit5());
    @assert SciMLBase.successful_retcode(sol)

    @test vec_sin_up ≈ [2, 4] atol=1e-3
    @test vec_sin_down ≈ [1, 3] atol=1e-3
    @test vec_cos_up ≈ [1.5, 3.5] atol=1e-3
    @test vec_cos_down ≈ [0.5, 2.5] atol=1e-3

    ####
    #### no negative affect
    ####
    empty!(vec_sin_up)
    empty!(vec_cos_up)
    empty!(vec_sin_down)
    empty!(vec_cos_down)
    cb_vec = VectorContinuousComponentCallback(cond_vec, affect_vec, 2; affect_neg! = nothing)
    prob = ODEProblem(nw, uflat(s0), (0, 4π+0.1), pflat(s0), callback=get_callbacks(nw, Dict(VIndex(1)=>cb_vec)))
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
    prob = ODEProblem(nw, uflat(s0), (0, 4π+0.1), pflat(s0), callback=get_callbacks(nw, Dict(VIndex(1)=>cb)))
    sol = solve(prob, Tsit5());
    @test sin_up ≈ [2, 4] atol=1e-3
    @test sin_down ≈ [1, 3] atol=1e-3

    empty!(sin_up)
    empty!(sin_down)
    cb = ContinuousComponentCallback(cond, affect; affect_neg! = nothing)
    prob = ODEProblem(nw, uflat(s0), (0, 4π+0.1), pflat(s0), callback=get_callbacks(nw, Dict(VIndex(1)=>cb)))
    sol = solve(prob, Tsit5());
    @test sin_up ≈ [2, 4] atol=1e-3
    @test isempty(sin_down)
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
