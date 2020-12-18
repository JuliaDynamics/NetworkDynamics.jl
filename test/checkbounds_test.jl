using Test
using LightGraphs
using NetworkDynamics

N = 10
g = erdos_renyi(10,10)

ov = ODEVertex(f! = (dv, v, edges, p, t) -> nothing, dim = 1)
se = StaticEdge(f! = (e, v_s, v_d, p, t) -> nothing, dim = 1, coupling = :undirected)
sde = StaticDelayEdge(se)
oe = ODEEdge(se)


nd! = network_dynamics(ov,se,g)

dx0 = zeros(N)
x0 = zeros(N)


@testset "Checking parameter bounds" begin
    @test nd!(dx0, x0, 1., 0.) === nothing
    @test nd!(dx0, x0, ones(N), 0.) === nothing


    @test nd!(dx0, x0, (1., nothing), 0.) === nothing
    @test nd!(dx0, x0, (ones(N), nothing), 0.) === nothing
    @test nd!(dx0, x0, (ones(N, N), nothing), 0.) === nothing
    @test nd!(dx0, x0, (ones(N + 1, N), nothing), 0.) === nothing
    @test nd!(dx0, x0, (ones(N - 1, N), nothing), 0.) === nothing

    @test_throws ErrorException nd!(dx0, x0, (ones(N + 1), nothing), 0.) === nothing
    @test_throws ErrorException nd!(dx0, x0, (ones(N - 1), nothing), 0.) === nothing
    @test_throws ErrorException nd!(dx0, x0, (ones(N, N + 1), nothing), 0.) === nothing
    @test_throws ErrorException nd!(dx0, x0, (ones(N, N - 1), nothing), 0.) === nothing

    @test nd!(dx0, x0, (nothing, 1.), 0.) === nothing
    @test nd!(dx0, x0, (nothing, ones(N)), 0.) === nothing
    @test nd!(dx0, x0, (nothing, ones(N, N)), 0.) === nothing
    @test nd!(dx0, x0, (nothing, ones(N + 1, N)), 0.) === nothing
    @test nd!(dx0, x0, (nothing, ones(N - 1, N)), 0.) === nothing

    @test_throws ErrorException nd!(dx0, x0, (nothing, ones(N + 1)), 0.) === nothing
    @test_throws ErrorException nd!(dx0, x0, (nothing, ones(N - 1)), 0.) === nothing
    @test_throws ErrorException nd!(dx0, x0, (nothing, ones(N, N + 1)), 0.) === nothing
    @test_throws ErrorException nd!(dx0, x0, (nothing, ones(N, N - 1)), 0.) === nothing
end

@testset "Checking variable arrays' bounds" begin
    @test_throws ErrorException nd!(ones(N), ones(N-1), 1., 0.)
    @test_throws ErrorException nd!(ones(N-1), ones(N), 1., 0.)

end
