using NetworkDynamics
using NetworkDynamics: ComponentGraph, legacy_graph_type, edge_multiplicity
using Graphs
using Graphs: SimpleEdge
using DataFrames
using InteractiveUtils: subtypes
using Test

@__MODULE__()==Main ? includet(joinpath(pkgdir(NetworkDynamics), "test", "ComponentLibrary.jl")) : (const Lib = Main.Lib)

dvertex(i) = VertexModel(f=Lib.diffusionvertex!, dim=1, g=1:1, vidx=i)
function dedge(src, dst; K=1.0, name=:diff_edge)
    EdgeModel(; g=AntiSymmetric(Lib.diffusionedge!), outdim=1, psym=[:K=>K], name, src, dst)
end
function odeedge(src, dst; name=:ode_edge)
    EdgeModel(; f=Lib.diffusion_dedge!, dim=2, sym=[:e_dst, :e_src], psym=[:τ=>100],
              g=Fiducial(dst=1:1, src=2:2), name, src, dst)
end
function rhs(nw)
    du = zeros(dim(nw))
    nw(du, collect(1.0:dim(nw)) .^ 2, pflat(NWParameter(nw)), 0.0)
    du
end

@testset "ComponentGraph" begin
    g = ComponentGraph(3, SimpleEdge.([1=>2, 1=>2, 2=>1, 3=>3]))
    @test nv(g) == 3
    @test ne(g) == 4
    @test is_directed(g)
    @test has_edge(g, 1, 2) && has_edge(g, 2, 1) && !has_edge(g, 2, 3)
    @test outneighbors(g, 1) == [2]
    @test inneighbors(g, 2) == [1]
    @test all_neighbors(g, 1) == [2]
    @test outdegree(g, 1) == 2
    @test indegree(g, 1) == 1
    @test degree(g, 1) == 3
    @test degree(g, 3) == 2 # self-loop counts twice
    @test degree(g) == [3, 3, 2]
    @test outdegree(g) == [2, 1, 1]
    @test indegree(g) == [1, 2, 1]
    @test !has_edge(g, 1, 5) && !has_edge(g, 0, 1)
    @test has_edge(g, 3, 3)
    @test g == ComponentGraph(3, collect(edges(g)))
    @test g != ComponentGraph(3, reverse(collect(edges(g))))
    @test hash(g) == hash(ComponentGraph(3, collect(edges(g))))
    @test_throws ArgumentError ComponentGraph(2, [SimpleEdge(1, 3)])

    @test legacy_graph_type(g) == :none
    @test legacy_graph_type(ComponentGraph(3, SimpleEdge.([1=>2, 2=>3]))) == :simple
    @test legacy_graph_type(ComponentGraph(3, SimpleEdge.([1=>2, 2=>1]))) == :digraph
    @test edge_multiplicity(edges(g)) == [(1,2), (2,2), (1,1), (1,1)]
end

@testset "graphless constructor" begin
    vs = dvertex.(1:4)
    # legacy path is unchanged: sorted into SimpleDiGraph order with warning
    es = [dedge(3, 1), dedge(1, 2), dedge(4, 1)]
    nw = @test_logs (:warn, r"Order of edge models") Network(vs, es)
    @test nw.im.g isa SimpleDiGraph
    @test nw.im.edgem == es[[2, 1, 3]]
    nw = @test_logs Network(vs, es[[2, 1, 3]])
    @test nw.im.g isa SimpleDiGraph
    @test Network(vs, [dedge(1, 2), dedge(1, 3)]).im.g isa SimpleGraph

    # opting out keeps input order
    nw = @test_logs Network(vs, es; legacy_graph=false)
    @test nw.im.g isa ComponentGraph
    @test nw.im.edgem == es

    # parallel edges keep input order and orientation, no warning
    es = [dedge(1, 2), dedge(1, 2), dedge(3, 1), dedge(4, 1), dedge(2, 1)]
    nw = @test_logs Network(vs, es)
    @test nw.im.g isa ComponentGraph
    @test nw.im.edgem == es
    @test nw.im.edgevec == SimpleEdge.([1=>2, 1=>2, 3=>1, 4=>1, 2=>1])

    # copy and explicit graph constructor keep the ComponentGraph
    @test copy(nw).im.g == nw.im.g
    nw2 = Network(ComponentGraph(4, nw.im.edgevec), vs, es)
    @test nw2.im.g == nw.im.g
end

@testset "parallel edges are independent" begin
    vs = dvertex.(1:3)
    # 1=>2 twice plus reversed 2=>1 behaves like a single edge with summed conductance
    es = [dedge(1, 2; K=1.0), dedge(1, 2; K=2.0), dedge(2, 1; K=0.5), dedge(3, 1; K=1.0)]
    ref = [dedge(1, 2; K=3.5), dedge(3, 1; K=1.0)]
    @test rhs(Network(vs, es)) ≈ rhs(Network(vs, ref; warn_order=false))

    # same result for all aggregators
    for accT in subtypes(NetworkDynamics.Aggregator)
        @test rhs(Network(vs, es; aggregator=accT(+))) ≈ rhs(Network(vs, ref; warn_order=false))
    end

    # states and parameters per branch
    nw = Network(dvertex.(1:2), [odeedge(1, 2; name=:line_a), odeedge(1, 2; name=:line_b)])
    s = NWState(nw)
    s.e[1, :e_dst] = 1.0
    s.e[2, :e_dst] = 2.0
    s.p.e[:line_b, :τ] = 7.0
    @test s[EIndex(:line_a, :e_dst)] == 1.0
    @test s[EIndex(:line_b, :e_dst)] == 2.0
    @test s.p.e[1, :τ] == 100
    @test s.p.e[2, :τ] == 7.0

    # jacobian sparsity detection
    nw = Network(vs, es; sparse=true)
    @test !isnothing(get_jac_prototype(nw))
end

@testset "endpoint lookup" begin
    es = [odeedge(1, 2), odeedge(1, 2), odeedge(3, 1)]
    nw = Network(dvertex.(1:3), es)
    s = NWState(nw)
    @test_throws "edge indices [1, 2]" nw[EIndex(1=>2)]
    @test_throws "edge indices [1, 2]" s[EIndex(1=>2, :e_dst)]
    @test nw[EIndex(3=>1)] === es[3]
    @test nw[EIndex(2)] === es[2]
    @test_throws "reverse edge" nw[EIndex(1=>3)]
end

@testset "describe_edges marks parallel edges" begin
    nw = Network(dvertex.(1:3), [dedge(1, 2), dedge(1, 2), dedge(3, 1)])
    @test describe_edges(nw).parallel == ["1/2", "2/2", ""]
    nw = Network(dvertex.(1:3), [dedge(1, 2), dedge(3, 1)])
    @test :parallel ∉ propertynames(describe_edges(nw))
end

@testset "interface defaults need same edge order" begin
    vs = dvertex.(1:3)
    nw1 = Network(vs, [dedge(1, 2), dedge(3, 1)]; legacy_graph=false)
    nw2 = Network(vs, [dedge(3, 1), dedge(1, 2)]; legacy_graph=false)
    s = NWState(nw1, [1.0, 2.0, 3.0], pflat(NWParameter(nw1)))
    @test_throws ArgumentError set_interface_defaults!(nw2, s)
    @test NetworkDynamics._same_topology(nw1, copy(nw1))
    @test !NetworkDynamics._same_topology(nw1, nw2)
end

@testset "loopback vertex has exactly one edge" begin
    hub = dvertex(1)
    other = dvertex(2)
    sat = VertexModel(g=(out, ins, p, t) -> (out[1] = ins[1]), outsym=[:u], insym=[:i], vidx=3)
    lb() = LoopbackConnection(potential=[:u], flow=[:i], src=3, dst=1)

    Network([hub, other, sat], [dedge(1, 2), lb()]; warn_order=false)
    @test_throws ArgumentError Network([hub, other, sat], [dedge(1, 2), lb(), lb()])
    @test_throws ArgumentError Network([hub, other, sat], [dedge(1, 2), lb(), dedge(2, 3)]; warn_order=false)
end
