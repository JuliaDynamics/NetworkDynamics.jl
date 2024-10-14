using CUDA
using CUDA.CUSPARSE
using Adapt
using NetworkDynamics
using StableRNGs
using KernelAbstractions
using Graphs
using Random
using Test
(isinteractive() && @__MODULE__()==Main ? includet : include)("ComponentLibrary.jl")

rng = StableRNG(1)
g = complete_graph(4)
vf = [Lib.kuramoto_second(), Lib.diffusion_vertex(), Lib.kuramoto_second(), Lib.diffusion_vertex()]
ef = [Lib.diffusion_odeedge(),
      Lib.kuramoto_edge(),
      Lib.kuramoto_edge(),
      Lib.diffusion_edge_fid(),
      Lib.diffusion_odeedge(),
      Lib.diffusion_edge_fid()]
nw = Network(g, vf, ef)
@test_throws ArgumentError adapt(CuArray{Float32}, nw) # wrong execution

nw = Network(g, vf, ef; execution=KAExecution{true}(), aggregator=KAAggregator(+))

x0 = rand(rng, dim(nw))
p  = rand(rng, pdim(nw))
dx = zeros(length(x0))
nw(dx, x0, p, NaN)

to = CUDABackend()
@test_throws ArgumentError adapt(to, nw)
to = :foo
@test_throws ArgumentError adapt(to, nw)
to = CuArray([1,2,3])
@test_throws ArgumentError adapt(to, nw)
to = cu(rand(3))
@test_throws ArgumentError adapt(to, nw)
to = CuArray
@test_throws ArgumentError adapt(to, nw)

to1 = CuArray{Float32}
to2 = CuArray{Float64}
nw1 = adapt(to1, nw)
nw2 = adapt(to2, nw)

for nw in (nw1, nw2)
    @test nw.vertexbatches[1].indices isa CuArray{Int}
    @test nw.layer.edgebatches[1].indices isa CuArray{Int}
    @test nw.gbufprovider.map isa CuArray{Int}
    @test nw.layer.aggregator.m.map isa CuArray{Int}
    @test nw.layer.aggregator.m.symmap isa CuArray{Int}
end

@test nw1.caches.state.du isa CuArray{Float32}
@test nw1.caches.aggregation.du isa CuArray{Float32}
@test nw2.caches.state.du isa CuArray{Float64}
@test nw2.caches.aggregation.du isa CuArray{Float64}

x0_d1 = adapt(to1, x0)
p_d1 = adapt(to1, p)
dx_d1 = adapt(to1, zeros(length(x0)))
nw1(dx_d1, x0_d1, p_d1, NaN)
@test Vector(dx_d1) ≈ dx
@test_throws ArgumentError nw2(dx_d1, x0_d1, p_d1, NaN) # wrong type for cache

x0_d2 = adapt(to2, x0)
p_d2 = adapt(to2, p)
dx_d2 = adapt(to2, zeros(length(x0)))
nw2(dx_d2, x0_d2, p_d2, NaN)
@test Vector(dx_d2) ≈ dx
@test_throws ArgumentError nw1(dx_d2, x0_d2, p_d2, NaN) # wrong type for cache


# try SparseAggregator
nw2 = Network(g, vf, ef; execution=KAExecution{true}(), aggregator=SparseAggregator(+))
nw2_d = adapt(CuArray{Float32}, nw2)

@test nw2_d.layer.aggregator.m isa CuSparseMatrixCSC
fill!(dx_d1, 0)
nw2_d(dx_d1, x0_d1, p_d1, NaN)
@test Vector(dx_d1) ≈ dx

# mini benchmark

#=
Ns = Int[1e3, 1e4, 1e5, 1e6, 1e7]
cput = Float64[]
gput = Float64[]
gput32 = Float64[]
for N in Ns
    rng = StableRNG(1)
    g = watts_strogatz(N, 3, 0.8; rng=StableRNG(1))
    edge = Lib.kuramoto_edge()
    vertex = [Lib.kuramoto_second(), Lib.kuramoto_first()]
    vertices = vertex[shuffle(rng, vcat([1 for _ in 1:N÷2], [2 for _ in 1:N÷2]))]

    nw = Network(g, vertices, edge; execution=KAExecution{true}(), aggregator=KAAggregator(+))

    @info "N = $N"
    x0 = rand(rng, dim(nw));
    p  = rand(rng, pdim(nw));
    dx = zeros(length(x0));
    b = @b nw(dx, x0, p, NaN) seconds=1
    push!(cput, b.time)

    to = CUDABackend();
    nw_d = adapt(to, nw);
    x0_d = adapt(to, x0);
    p_d = adapt(to, p);
    dx_d = adapt(to, zeros(length(x0)));
    b = @b nw_d(dx_d, x0_d, p_d, NaN) seconds=1
    push!(gput, b.time)

    to = CUDABackend()
    nw_d = adapt(to, nw)
    x0_d = adapt(to, Float32.(x0))
    p_d = adapt(to, Float32.(p))
    dx_d = adapt(to, zeros(Float32, length(x0)))
    b = @b nw_d(dx_d, x0_d, p_d, NaN) seconds=1
    push!(gput32, b.time)
    GC.gc()
end

using GLMakie
fig = Figure()
ax = Axis(fig[1,1]; xscale=log10, yscale=log10)
Makie.scatterlines!(ax, Ns, cput)
Makie.scatterlines!(ax, Ns, gput)
Makie.scatterlines!(ax, Ns, gput32)
=#
