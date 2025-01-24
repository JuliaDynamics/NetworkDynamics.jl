using NetworkDynamics
using Graphs
using OrdinaryDiffEqTsit5
using Chairmarks
using Test
using ModelingToolkit

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

    cond = ComponentCondition([:P, :₋P, :srcθ], [:limit, :K]) do u, p, t
        u[:P] - p[:limit]
    end
    affect = ComponentAffect([],[:active]) do u, p, ctx
        @info "Trip line between $(ctx.src) and $(ctx.dst)"
        p[:active] = 0
    end
    cb = ContinousComponentCallback(cond, affect)
    for i in 1:5
        set_callback!(ls[i], cb)
    end

    batches = NetworkDynamics.collect_callbackbatches(nw);
    @test length(batches) == 1
    cbb = only(batches);

    @test NetworkDynamics.condition_dim(cbb) == 3
    @test NetworkDynamics.condition_pdim(cbb) == 2
    @test NetworkDynamics.affect_dim(cbb) == 0
    @test NetworkDynamics.affect_pdim(cbb) == 1

    @test NetworkDynamics.condition_urange.(Ref(cbb), 1:length(cbb)) == [1:3,4:6,7:9,10:12,13:15]
    @test NetworkDynamics.condition_prange.(Ref(cbb), 1:length(cbb)) == [1:2,3:4,5:6,7:8,9:10]
    @test NetworkDynamics.affect_urange.(Ref(cbb), 1:length(cbb)) == [1:0 for i in 1:5]
    @test NetworkDynamics.affect_prange.(Ref(cbb), 1:length(cbb)) == [1:1,2:2,3:3,4:4,5:5]
    @test NetworkDynamics.condition_outrange.(Ref(cbb), 1:length(cbb)) == [1:1,2:2,3:3,4:4,5:5]

    @test NetworkDynamics.collect_c_or_a_indices(cbb, :condition, :sym) == collect(Iterators.flatten(collect(EIndex(i, [:P, :₋P, :srcθ])) for i in 1:5))
    @test NetworkDynamics.collect_c_or_a_indices(cbb, :condition, :psym) == collect(Iterators.flatten(collect(EPIndex(i, [:limit, :K])) for i in 1:5))
    @test NetworkDynamics.collect_c_or_a_indices(cbb, :affect, :sym) == []
    @test NetworkDynamics.collect_c_or_a_indices(cbb, :affect, :psym) == collect(EPIndex(i, :active) for i in 1:5)

    batchcond = NetworkDynamics.batch_condition(cbb)
    out = zeros(5)
    fill!(out, NaN)
    u = uflat(s0)
    integrator = (; p = pflat(s0))
    b = @b $batchcond($out, $u, NaN, $integrator)
    @test b.allocs == 0

    batchaffect = NetworkDynamics.batch_affect(cbb)
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
