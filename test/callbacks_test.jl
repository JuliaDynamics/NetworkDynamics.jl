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

@testset "basic callback tests" begin
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
    cb = ContinousComponentCallback(cond, affect)
    set_callback!.(nw.im.edgem, Ref(cb))

    batches = NetworkDynamics.collect_callbackbatches(nw);
    @test length(batches) == 1
    cbb = only(batches);

    @test NetworkDynamics.condition_dim(cbb) == 3
    @test NetworkDynamics.condition_pdim(cbb) == 2
    @test NetworkDynamics.affect_dim(cbb) == 0
    @test NetworkDynamics.affect_pdim(cbb) == 1

    @test NetworkDynamics.condition_urange.(Ref(cbb), 1:length(cbb)) == [1:3,4:6,7:9,10:12,13:15,16:18,19:21]
    @test NetworkDynamics.condition_prange.(Ref(cbb), 1:length(cbb)) == [1:2,3:4,5:6,7:8,9:10,11:12,13:14]
    @test NetworkDynamics.affect_urange.(Ref(cbb), 1:length(cbb)) == [1:0 for i in 1:7]
    @test NetworkDynamics.affect_prange.(Ref(cbb), 1:length(cbb)) == [1:1,2:2,3:3,4:4,5:5,6:6,7:7]
    @test NetworkDynamics.condition_outrange.(Ref(cbb), 1:length(cbb)) == [1:1,2:2,3:3,4:4,5:5,6:6,7:7]

    @test NetworkDynamics.collect_c_or_a_indices(cbb, :condition, :sym) == collect(Iterators.flatten(collect(EIndex(i, [:P, :₋P, :srcθ])) for i in 1:7))
    @test NetworkDynamics.collect_c_or_a_indices(cbb, :condition, :psym) == collect(Iterators.flatten(collect(EPIndex(i, [:limit, :K])) for i in 1:7))
    @test NetworkDynamics.collect_c_or_a_indices(cbb, :affect, :sym) == []
    @test NetworkDynamics.collect_c_or_a_indices(cbb, :affect, :psym) == collect(EPIndex(i, :active) for i in 1:7)

    batchcond = NetworkDynamics.batch_condition(cbb)
    out = zeros(7)
    fill!(out, NaN)
    s0 = NWState(nw)
    u = uflat(s0)
    integrator = (; p = pflat(s0))
    b = @b $batchcond($out, $u, NaN, $integrator)
    @test b.allocs == 0


    trip_first_cb = PresetTimeCallback(1.0, integrator -> begin
        i = 5
        @info "Trip initial line $i at t=$(integrator.t)"
        p = NWParameter(integrator)
        p.e[i,:active] = 0
        auto_dt_reset!(integrator)
        save_parameters!(integrator)
    end)
    nwcb = NetworkDynamics.get_callbacks(nw)
    s0 = NWState(nw)
    prob = ODEProblem(nw, uflat(s0), (0,6), copy(pflat(s0)), callback=CallbackSet(trip_first_cb, nwcb))
    sol = solve(prob, Tsit5());

    @test tripi == [7,4,1,3,2]
    tref = [2.247676397005474, 2.502523192233235, 3.1947647115093654, 3.3380530127462587, 3.4042696241577888]
    @test maximum(abs.(tript - tref)) < 1e-5
end

@testset "wrong symboltype test" begin
    nw = basenetwork()

    # invalid pram in condition u
    # invalid obs in affect u
    cond = ComponentCondition((args...)->nothing, [:P, :₋P, :srcθ, :limit], [:limit, :K])
    affect = ComponentAffect((args...)->nothing, [:₋P],[:active])
    cb = ContinousComponentCallback(cond, affect)
    set_callback!.(nw.im.edgem, Ref(cb); check=false);
    cbb = only(NetworkDynamics.collect_callbackbatches(nw));
    @test_throws ArgumentError NetworkDynamics.batch_condition(cbb)
    @test_throws ArgumentError NetworkDynamics.batch_affect(cbb)

    # invalid state in condition p
    cond = ComponentCondition((args...)->nothing, [:P, :₋P, :srcθ], [:limit, :K, :P])
    affect = ComponentAffect((args...)->nothing, [],[:active, :₋P])
    cb = ContinousComponentCallback(cond, affect)
    set_callback!.(nw.im.edgem, Ref(cb); check=false);
    cbb = only(NetworkDynamics.collect_callbackbatches(nw));
    @test_throws ArgumentError NetworkDynamics.batch_condition(cbb)
    @test_throws ArgumentError NetworkDynamics.batch_affect(cbb)

    # test on set_callback
    cond = ComponentCondition((args...)->nothing, [:P, :₋P, :srcθ, :limit], [:limit, :K])
    affect = ComponentAffect((args...)->nothing, [],[:active])
    cb = ContinousComponentCallback(cond, affect)
    @test_throws ArgumentError set_callback!(nw.im.edgem[1], cb)
    cond = ComponentCondition((args...)->nothing, [:P, :₋P, :srcθ], [:limit, :K, :P])
    affect = ComponentAffect((args...)->nothing, [],[:active])
    cb = ContinousComponentCallback(cond, affect)
    @test_throws ArgumentError set_callback!(nw.im.edgem[1], cb)
    cond = ComponentCondition((args...)->nothing, [:P, :₋P, :srcθ], [:limit, :K])
    affect = ComponentAffect((args...)->nothing, [:₋P],[:active])
    cb = ContinousComponentCallback(cond, affect)
    @test_throws ArgumentError set_callback!(nw.im.edgem[1], cb)
    cond = ComponentCondition((args...)->nothing, [:P, :₋P, :srcθ], [:limit, :K])
    affect = ComponentAffect((args...)->nothing, [],[:active, :P])
    cb = ContinousComponentCallback(cond, affect)
    @test_throws ArgumentError set_callback!(nw.im.edgem[1], cb)
end
