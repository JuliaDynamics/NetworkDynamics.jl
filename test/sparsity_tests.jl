using NetworkDynamics
using SparseConnectivityTracer
using Graphs
using OrdinaryDiffEqRosenbrock
using ModelingToolkit
using ModelingToolkit: D_nounits as Dt, t_nounits as t
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
    j2 = get_jac_prototype(nw; dense=[EIndex(1), VIndex(1)])
    j3 = get_jac_prototype(nw; dense=true)
    j4 = get_jac_prototype(nw; dense=vcat([EIndex(i) for i in 1:ne(nw)], [VIndex(i) for i in 1:nv(nw)]))

    # test that they become strictly more dense
    @test j3 == j4
    for idx in eachindex(j1)
        if j1[idx] != 0
            @test j1[idx] == j2[idx] == j3[idx]
        end
    end
    for idx in eachindex(j2)
        if j2[idx] != 0
            @test j2[idx] == j3[idx]
        end
    end

    # prob = ODEProblem(nw, x0, (0.0, 1.0), p0)
    # _nw = ODEFunction(nw; jac_prototype=get_jac_prototype(nw))
    # prob_jac = ODEProblem(_nw, x0, (0.0, 1.0), p0)
    # @b solve($prob, $Rodas5P())
    # @b solve($prob_jac, $Rodas5P())
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
    @test_throws ErrorException get_jac_prototype(nw) # fails because of the conditional!
    get_jac_prototype(nw; remove_conditions=true) # should work now
    get_jac_prototype(nw; dense=true) # should work now
end


@testset "test filter conditionals" begin
    compare_expr(a, b) = Base.remove_linenums!(a) == Base.remove_linenums!(b)

    assigment = :(dest = if cond; truepath; else; falsepath; end)
    target = :(dest = begin
        truepath
    end + begin
        falsepath
    end)
    @test compare_expr(SE.filter_conditionals_expr(assigment), target)

    with_elseif = :(if cond; truepath; elseif cond2; true2; else; falsepath; end)
    @test_throws SE.RemainingConditionalsException SE.filter_conditionals_expr(with_elseif)
end
