using Graphs

@testset "test syms containing" begin
    vert = ODEVertex(f = nothing, sym = [:a1, :b2], dim=2)
    edge = StaticEdge(f = nothing, dim=2)
    g = complete_graph(3)
    nd = network_dynamics(vert, edge, g)

    @test idx_containing(nd, "1") == [1, 2, 3, 5]
    @test syms_containing(nd, "1") == nd.syms[idx_containing(nd, "1")]

    @test idx_containing(nd, :a) == [1,3,5]
    @test syms_containing(nd, :a) == nd.syms[idx_containing(nd, :a)]

    @test idx_containing(nd, r"_1$") == [1,2]
    @test syms_containing(nd, r"_1$") == nd.syms[idx_containing(nd, r"_1$")]
end
