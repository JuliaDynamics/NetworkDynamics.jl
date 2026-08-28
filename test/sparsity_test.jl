using NetworkDynamics
using SparseArrays
using SparseConnectivityTracer
using Graphs
using OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqNonlinearSolve
using NonlinearSolve: NewtonRaphson
using LinearAlgebra: Diagonal
using SteadyStateDiffEq: SteadyStateDiffEq, SteadyStateProblem
using SciMLBase: solve
using ModelingToolkitBase
using ModelingToolkitBase: D_nounits as Dt, t_nounits as t
using SciCompDSL
using InteractiveUtils: subtypes
SE = Base.get_extension(NetworkDynamics, :NetworkDynamicsSparsityExt)

@__MODULE__()==Main ? includet(joinpath(pkgdir(NetworkDynamics), "test", "ComponentLibrary.jl")) : (const Lib = Main.Lib)

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
    Symbolics.@register_symbolic true_if_else_block(cond, t, f)
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
    # target = :(dest = ifelse(cond, begin
    #     truepath
    # end, begin
    #     falsepath
    # end))
    target = :(dest = begin
        truepath
    end + begin
        falsepath
    end)
    @test compare_expr(SE.filter_conditionals_expr(assigment), target)

    with_elseif = :(if cond; truepath; elseif cond2; true2; else; falsepath; end)
    @test_throws SE.RemainingConditionalsException SE.filter_conditionals_expr(with_elseif)
end

@testset "fixpoint solver selection" begin
    g = watts_strogatz(100, 4, 0.1; seed=1)
    v = VertexModel(; f=(dv, x, esum, p, t) -> (dv[1] = p[1] - x[1]^3 + esum[1]; nothing),
                    g=1:1, dim=1, pdim=1, sym=[:x], psym=[:p => 0.1])
    e = EdgeModel(; g=AntiSymmetric((e, vs, vd, p, t) -> (e[1] = p[1]*(vs[1] - vd[1]); nothing)),
                  outdim=1, pdim=1, psym=[:k => 1.0])
    nw = Network(g, v, e)

    @test NetworkDynamics.default_fixpoint_alg(nw) isa SteadyStateDiffEq.SSRootfind
    s_ref = find_fixpoint(nw)

    set_jac_prototype!(nw; verbose=false)
    @test !(NetworkDynamics.default_fixpoint_alg(nw) isa SteadyStateDiffEq.SSRootfind)
    @test uflat(find_fixpoint(nw)) ≈ uflat(s_ref) rtol=1e-6

    # nothing means "you decide", an explicit alg still wins
    @test uflat(find_fixpoint(nw; alg=nothing)) ≈ uflat(s_ref) rtol=1e-6
    @test uflat(find_fixpoint(nw; alg=SteadyStateDiffEq.SSRootfind())) ≈ uflat(s_ref) rtol=1e-6

    # the chosen alg must actually build jacobians rather than falling into Broyden
    s0 = NWState(nw; ufill=0.1)
    prob = SteadyStateProblem(NetworkDynamics.NetworkFixedT(nw, NaN), uflat(s0), pflat(s0))
    @test solve(prob, NetworkDynamics.default_fixpoint_alg(nw)).stats.njacs > 0
    @test solve(prob, SteadyStateDiffEq.SSRootfind()).stats.njacs == 0

    # ...and it must do so with sparse coloring, not dense AD
    @test_logs solve(prob, NetworkDynamics.default_fixpoint_alg(nw))
end


@testset "DAE initialization solver selection" begin
    # one differential state x and one algebraic state y per vertex
    g = watts_strogatz(100, 4, 0.1; seed=1)
    vf = function(dv, v, esum, p, t)
        dv[1] = -v[1] + v[2] + esum[1]
        dv[2] = v[2]^3 + v[2] - v[1] - p[1]
        nothing
    end
    v = VertexModel(; f=vf, g=1:1, dim=2, pdim=1, sym=[:x, :y], psym=[:p => 0.5],
                    mass_matrix=Diagonal([1.0, 0.0]))
    e = EdgeModel(; g=AntiSymmetric((e, vs, vd, p, t) -> (e[1] = p[1]*(vs[1] - vd[1]); nothing)),
                  outdim=1, pdim=1, psym=[:k => 1.0])
    nw = Network(g, v, e)

    @test NetworkDynamics.default_dae_init_alg(nw) isa BrownFullBasicInit
    @test isnothing(NetworkDynamics.default_dae_init_alg(nw).nlsolve)

    set_jac_prototype!(nw; verbose=false)
    @test !isnothing(NetworkDynamics.default_dae_init_alg(nw).nlsolve)

    # y=0 is inconsistent, so the initialization actually has to solve
    s0 = NWState(nw)
    s0.v[:, :x] .= 0.1
    s0.v[:, :y] .= 0.0

    integ = init(ODEProblem(nw, uflat(s0), (0.0, 1.0), pflat(s0)), Rodas5P())
    y = integ.u[2:2:end]
    x = integ.u[1:2:end]
    @test maximum(abs, y.^3 .+ y .- x .- 0.5) < 1e-12   # constraints satisfied

    # Upstream (OrdinaryDiffEqNonlinearSolve): `_initialize_dae!` decides whether to
    # allocate dual work buffers with `alg_autodiff(alg) isa AutoForwardDiff`, which is
    # false once `prepare_user_sparsity` has wrapped the AD type in `AutoSparse` -- i.e.
    # exactly when a `jac_prototype` is present. A ForwardDiff-based `nlsolve` then writes
    # Duals into Float64 buffers. `default_dae_init_alg` pins `AutoFiniteDiff` because of
    # this. When this test stops being broken, drop that pin.
    @test_broken begin
        prob = ODEProblem(nw, uflat(s0), (0.0, 1.0), pflat(s0);
                          initializealg=BrownFullBasicInit(nlsolve=NewtonRaphson()))
        init(prob, Rodas5P())
        true
    end
end
