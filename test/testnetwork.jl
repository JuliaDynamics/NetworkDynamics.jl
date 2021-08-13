using NDPrototype
using LightGraphs
using Random

#############
include("../benchmark/benchmarks.jl")
include("../benchmark/benchmark_utils.jl")

# diffusion test
N = 10
vertex = diffusion_vertex()
edge = diffusion_edge()
@benchmark begin
    nd(dx, x0, nothing, 0.0)
end setup = begin
    g = watts_strogatz($N, 3, 0.8, seed=1)
    nd = Network(g, $vertex, $edge)
    x0 = randn(dim(nd))
    dx = similar(x0)
    nd(dx, x0, nothing, 0.0)
end

g = watts_strogatz(N, 3, 0.8, seed=1)
nd = Network(g, vertex, edge);
x0 = randn(dim(nd1)); dx = zero(x0);
@btime $nd($dx,$x0,nothing,0.0)

# parallel
@benchmark begin
    nd(dx, x0, nothing, 0.0)
end setup = begin
    g = watts_strogatz($N, 3, 0.8, seed=1)
    nd = Network(g, $vertex, $edge, parallel=true)
    x0 = randn(dim(nd))
    dx = similar(x0)
    nd(dx, x0, nothing, 0.0)
end

# odeedge
@benchmark begin
    nd(dx, x0, nothing, 0.0)
end setup = begin
    g = watts_strogatz($N, 3, 0.8, seed=1)
    edge = diffusion_dedge()
    nd = Network(g, $vertex, edge, accdim=1, parallel=true)
    x0 = randn(dim(nd))
    dx = similar(x0)
    nd(dx, x0, nothing, 0.0)
end

N = 1000
g = watts_strogatz(N, 3, 0.8, seed=1)
edge = diffusion_dedge()
vertex = diffusion_vertex()
nd1 = Network(g, vertex, edge, execution=:seq);
nd2 = Network(g, vertex, edge, execution=:threaded);
x0 = randn(dim(nd1));
dx1 = zero(x0);
dx2 = zero(x0);

@btime $nd1($dx1, $x0, nothing, 0.0)
@btime $nd2($dx2, $x0, nothing, 0.0)

nd2(dx2, x0, nothing, 0.0)
@test dx1 ≈ dx2
dx1

# inhomgeneous network
function heterogeneous(N)
    rng = MersenneTwister(1)
    g = watts_strogatz(N, 3, 0.8, seed=1)
    edge = static_kuramoto_edge()
    vertex = [kuramoto_vertex_1d(), kuramoto_vertex_2d()]
    vertices = vertex[shuffle(rng, vcat([1 for _ in 1:N÷2], [2 for _ in 1:N÷2]))]
    p = vcat(randn(rng, nv(g)), randn(rng, ne(g)))
    (p, vertices, edge, g)
end

N = 10
@benchmark begin
    nd(dx, x0, p, 0.0)
end setup = begin
    (p, v, e, g) = heterogeneous($N)
    nd = Network(g, v, e)
    @assert length(p) == pdim(nd)
    x0 = randn(dim(nd))
    dx = similar(x0)
    nd(dx, x0, p, 0.0) # call to init caches, we don't want to benchmark this
end
@benchmark begin
    nd(dx, x0, p, 0.0)
end setup = begin
    (p, v, e, g) = heterogeneous($N)
    nd = Network(g, v, e, parallel=true)
    @assert length(p) == pdim(nd)
    x0 = randn(dim(nd))
    dx = similar(x0)
    nd(dx, x0, p, 0.0) # call to init caches, we don't want to benchmark this
end


begin
    N = 10000
    (p, v, e, g) = heterogeneous(N)
    nd1 = Network(g, v, e, execution=:seq)
    (p, v, e, g) = heterogeneous(N)
    nd2 = Network(g, v, e, execution=:threaded)
    x0 = randn(dim(nd1))
    dx1 = zeros(dim(nd1))
    dx2 = zeros(dim(nd1))
    @time nd1(dx1, x0, p, 0.0) # call to init caches, we don't want to benchmark this
    @time nd1(dx2, x0, p, 0.0) # call to init caches, we don't want to benchmark this
    @test dx1 ≈ dx2

    @time nd2(dx1, x0, p, 0.0) # call to init caches, we don't want to benchmark this
    @time nd2(dx2, x0, p, 0.0) # call to init caches, we don't want to benchmark this
    @test dx1 ≈ dx2

    @btime $nd1($dx1, $x0, $p, 0.0)
  # 509.375 μs (6 allocations: 320 bytes)
    @btime $nd2($dx2, $x0, $p, 0.0)
  # 211.208 μs (222 allocations: 17.84 KiB)
end

