using NetworkDynamics
using SparseArrays
using SparseConnectivityTracer
using Graphs
using OrdinaryDiffEqRosenbrock
using ModelingToolkit
using ModelingToolkit: D_nounits as Dt, t_nounits as t
using InteractiveUtils: subtypes
SE = Base.get_extension(NetworkDynamics, :NetworkDynamicsSparsityExt)

(isinteractive() && @__MODULE__()==Main ? includet : include)("ComponentLibrary.jl")

@testset "basic tests" begin
    g = complete_graph(4)
    vf = [Lib.kuramoto_second(), Lib.diffusion_vertex(), Lib.kuramoto_second(), Lib.diffusion_vertex()]
    ef = [Lib.diffusion_odeedge(),
          Lib.kuramoto_edge(),
          Lib.kuramoto_edge(),
          Lib.diffusion_edge_fid(),
          Lib.diffusion_odeedge(),
          Lib.diffusion_edge_fid()]
    nw = Network(g, vf, ef)
    x0 = rand(dim(nw))
    _p0 = NWParameter(nw)
    _p0.e[2:3,:K] .= 1.0
    p0 = pflat(_p0)

    j1 = get_jac_prototype(nw)

    # prob = ODEProblem(nw, x0, (0.0, 1.0), p0)
    # _nw = ODEFunction(nw; jac_prototype=get_jac_prototype(nw))
    # prob_jac = ODEProblem(_nw, x0, (0.0, 1.0), p0)
    # @b solve($prob, $Rodas5P())
    # @b solve($prob_jac, $Rodas5P())

    @testset "test retrieval of jac prototype under different execution schemes" begin
        # KAAggregator is known to be broken
        _nw = Network(nw, aggregator=KAAggregator(+))
        @test_broken get_jac_prototype(_nw; make_compatible=false)

        styles = [KAExecution{true}(),
                  KAExecution{false}(),
                  SequentialExecution{true}(),
                  SequentialExecution{false}(),
                  PolyesterExecution{true}(),
                  PolyesterExecution{false}(),
                  ThreadedExecution{true}(),
                  ThreadedExecution{false}()]
        unmatchedstyles = filter(subtypes(NetworkDynamics.ExecutionStyle)) do abstractstyle
            !any(s -> s isa abstractstyle, styles)
        end
        @assert isempty(unmatchedstyles) "Some ExecutionStyle won't be tested: $unmatchedstyles"
        aggregators = [NetworkDynamics.NaiveAggregator,
                       KAAggregator,
                       SequentialAggregator,
                       PolyesterAggregator,
                       ThreadedAggregator,
                       SparseAggregator]
        unmatchedaggregators = filter(subtypes(NetworkDynamics.Aggregator)) do abstractaggregator
            !any(s -> s <: abstractaggregator, aggregators)
        end
        @assert isempty(unmatchedaggregators) "Some AggrgationStyle won't be tested: $unmatchedaggregators"

        exsaggs = [(ex, agg) for ex in styles for agg in aggregators]
        for (execution, aggregator) in exsaggs
            _nw = Network(nw; execution, aggregator=aggregator(+))
            j = get_jac_prototype(_nw)
            @test j == j1
        end
    end
end

@testset "test of remove conditional" begin
    @mtkmodel ConditionalNode begin
        @variables begin
            p(t)=1, [description="pressure at node"]
            q_nw(t), [description="flow from node to network"]
        end
        @parameters begin
            C=1, [description="capacitance of node"]
            q_external, [description="external flow into node"]
        end
        @equations begin
            C*Dt(p) ~ q_external + q_nw
        end
    end
    @named vmtk = ConditionalNode()

    @mtkmodel ValveToggle begin
        @variables begin
            p_src(t), [description="pressure at src"]
            p_dst(t), [description="pressure at dst"]
            q(t), [description="flow through valve"]
        end
        @parameters begin
            K=1, [description="conductance of valve"]
            active=1, [description="active state of valve"]
        end
        @equations begin
            q ~ ifelse(active > 0, K * (p_src - p_dst), 0)
        end
    end
    @named valvet_mtk = ValveToggle()

    g = wheel_graph(10)
    v = VertexModel(vmtk, [:q_nw], [:p])
    valvet = EdgeModel(valvet_mtk, [:p_src], [:p_dst], AntiSymmetric([:q]))

    nw = Network(g, v, valvet)
    j1 = get_jac_prototype(nw) # fails because of the conditional!

    nw_sparse = Network(nw)
    set_jac_prototype!(nw_sparse)

    s0 = NWState(nw_sparse)
    prob = ODEProblem(nw_sparse, uflat(s0), (0,10), pflat(s0))
    @test prob.f.jac_prototype === nw_sparse.jac_prototype # should be the same prototype

    function true_if_else_block(cond, t, f)
        if cond
            return t
        else
            return f
        end
    end
    ModelingToolkit.@register_symbolic true_if_else_block(cond, t, f)
    @mtkmodel ValveToggle2 begin
        @variables begin
            p_src(t), [description="pressure at src"]
            p_dst(t), [description="pressure at dst"]
            q(t), [description="flow through valve"]
        end
        @parameters begin
            K=1, [description="conductance of valve"]
            active=1, [description="active state of valve"]
        end
        @equations begin
            q ~ true_if_else_block(active > 0, K * (p_src - p_dst), 0)
        end
    end
    @named valvet2_mtk = ValveToggle2()
    valvet2 = EdgeModel(valvet2_mtk, [:p_src], [:p_dst], AntiSymmetric([:q]))
    nw2 = Network(g, v, valvet2)
    j2 = get_jac_prototype(nw2) # should fall bakc to dense
    @test j1 == j2 # no diff in that case
end


@testset "test filter conditionals" begin
    compare_expr(a, b) = Base.remove_linenums!(a) == Base.remove_linenums!(b)

    assigment = :(dest = if cond; truepath; else; falsepath; end)
    target = :(dest = ifelse(cond, begin
        truepath
    end, begin
        falsepath
    end))
    @test compare_expr(SE.filter_conditionals_expr(assigment), target)

    with_elseif = :(if cond; truepath; elseif cond2; true2; else; falsepath; end)
    @test_throws SE.RemainingConditionalsException SE.filter_conditionals_expr(with_elseif)
end
