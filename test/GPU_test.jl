using CUDA
using Adapt
using NetworkDynamics
using StableRNGs
using KernelAbstractions
using Graphs
using Random
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
nw = Network(g, vf, ef; execution=KAExecution{true}(), aggregator=KAAggregator(+))

x0 = rand(rng, dim(nw))
p  = rand(rng, pdim(nw))
dx = zeros(length(x0))
nw(dx, x0, p, NaN)

to = CUDABackend()
nw_d = adapt(to, nw)
@test nw_d.vertexbatches[1].indices isa CuArray
@test nw_d.layer.edgebatches[1].indices isa CuArray
@test nw_d.layer.gather_map isa CuArray
@test nw_d.layer.aggregator.m.map isa CuArray
@test nw_d.layer.aggregator.m.symmap isa CuArray
x0_d = adapt(to, x0)
p_d = adapt(to, p)
dx_d = adapt(to, zeros(length(x0)))
nw_d(dx_d, x0_d, p_d, NaN)
@test Vector(dx_d) ≈ dx

# mini benchmark

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
