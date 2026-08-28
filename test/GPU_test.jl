using CUDA
using CUDA.CUSPARSE
using CUDSS # LinearSolve only offers a direct sparse device solve once this is loaded
using Adapt
using NetworkDynamics
using StableRNGs
using KernelAbstractions
using Graphs
using Random
using Test
using SparseConnectivityTracer
using OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqNonlinearSolve
using NonlinearSolve # for AutoFiniteDiff
using NonlinearSolve: NewtonRaphson, KrylovJL_GMRES, LUFactorization, QRFactorization
using OrdinaryDiffEqSDIRK: KenCarp4, TRBDF2
using OrdinaryDiffEqNonlinearSolve: BrownFullBasicInit
using SciMLBase
using SparseArrays
@__MODULE__()==Main ? includet(joinpath(pkgdir(NetworkDynamics), "test", "ComponentLibrary.jl")) : (const Lib = Main.Lib)

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
end

@test nw1.caches.output.du isa CuArray{Float32}
@test nw1.caches.aggregation.du isa CuArray{Float32}
@test nw2.caches.output.du isa CuArray{Float64}
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

@testset "MTK based model" begin
    g = cycle_graph(5) # 5-node cycle graph
    v1 = Lib.dqbus_swing_and_load()
    v2 = Lib.dqbus_swing()
    v3 = Lib.dqbus_pq()
    v4 = Lib.dqbus_pq()
    v5 = Lib.dqbus_pq()
    e = Lib.dqline(X=0.1, R=0.01)
    nw = Network(g, [v1, v2, v3, v4, v5], e; dealias=true, execution=KAExecution{true}(), aggregator=KAAggregator(+))
    nw_d = adapt(CuArray{Float32}, nw)
    @test nw_d.vertexbatches[1].compf.body == nothing
    du = zeros(dim(nw))
    u = ones(dim(nw))
    p = ones(pdim(nw))
    du_d = adapt(CuArray{Float32}, du)
    u_d = adapt(CuArray{Float32}, u)
    p_d = adapt(CuArray{Float32}, p)
    t = NaN
    nw_d(du_d, u_d, p_d, t)
    nw(du, u, p, t)
    @test maximum(abs.(du_d .- cu(du))) < 1e-6

    # test state adaption
    s0 = NWState(nw)
    adapted_s0 = adapt(CuArray{Float32}, s0)
    @test uflat(adapted_s0) isa CuArray{Float32}
    @test pflat(adapted_s0) isa CuArray{Float32}
end

# One device solve, so that every configuration below is a single test line rather than a
# commented-out block. Returns whether the solve succeeded; a configuration that errors is
# reported by `@test_broken` just like one that returns false.
#
# `layout` picks the Jacobian prototype: `:dense` sets none at all, `:csr` is what `adapt`
# gives and the only layout with a direct sparse solver on the device, `:csc` is the other
# cuSPARSE layout, which always falls back to Krylov.
#
# `maxiters` is far below the solver default on purpose. A configuration that crawls should
# come back as a failure rather than eat the whole CI budget.
function run_on_gpu(nw, s0; T, layout, solver, linsolve=nothing, initializealg=nothing,
                    tspan=(0.0, 1.0), maxiters=2000)
    nw_d = adapt(CuArray{T}, layout === :dense ? copy(nw) : set_jac_prototype!(copy(nw)))
    layout === :csc && set_jac_prototype!(nw_d, CuSparseMatrixCSC(nw_d.jac_prototype))
    prob = ODEProblem(nw_d, adapt(CuArray{T}, s0), tspan)
    alg = isnothing(linsolve) ? solver() : solver(; linsolve)
    kwargs = isnothing(initializealg) ? (;) : (; initializealg)
    sol = solve(prob, alg; maxiters, kwargs...)
    SciMLBase.successful_retcode(sol)
end

# Every entry is one test line; `broken` marks what is known not to work today. The list is
# deliberately shorter than the full product of the options — each new combination compiles a
# fresh solver specialisation, which is what the runtime here is spent on.
gpu_configs = [
    (; name="Float64 dense Rodas5P",     T=Float64, layout=:dense, solver=Rodas5P,  broken=false),
    (; name="Float64 dense Rodas5P GMRES", T=Float64, layout=:dense, solver=Rodas5P, broken=false,
       linsolve=KrylovJL_GMRES()),
    (; name="Float64 dense KenCarp4",    T=Float64, layout=:dense, solver=KenCarp4, broken=false),
    (; name="Float64 dense TRBDF2 GMRES", T=Float64, layout=:dense, solver=TRBDF2,  broken=false,
       linsolve=KrylovJL_GMRES()),
    # CSC has no direct solver on the device at all: LinearSolve warns and hands over to Krylov.
    (; name="Float64 CSC Rodas5P",       T=Float64, layout=:csc,   solver=Rodas5P,  broken=false),
    (; name="Float64 CSC KenCarp4",      T=Float64, layout=:csc,   solver=KenCarp4, broken=false),
    # CSR is what `adapt` produces. With CUDSS loaded the default linsolve is a real cuDSS LU,
    # and `QRFactorization` reaches cuSOLVER's sparse QR.
    (; name="Float64 CSR Rodas5P cuDSS", T=Float64, layout=:csr,   solver=Rodas5P,  broken=false),
    (; name="Float64 CSR Rodas5P cuDSS LU", T=Float64, layout=:csr, solver=Rodas5P, broken=false,
       linsolve=LUFactorization()),
    (; name="Float64 CSR Rodas5P cuSOLVER QR", T=Float64, layout=:csr, solver=Rodas5P, broken=false,
       linsolve=QRFactorization()),
    (; name="Float64 CSR KenCarp4",      T=Float64, layout=:csr,   solver=KenCarp4, broken=false),
    (; name="Float64 CSR TRBDF2 GMRES",  T=Float64, layout=:csr,   solver=TRBDF2,   broken=false,
       linsolve=KrylovJL_GMRES()),
    # A sparse prototype plus an explicit init nlsolve dies compiling a device broadcast that
    # stores a `Dual` into a float view. Independent of the eltype.
    (; name="Float64 CSR Rodas5P NewtonRaphson init", T=Float64, layout=:csr, solver=Rodas5P,
       broken=true, initializealg=BrownFullBasicInit(nlsolve=NewtonRaphson())),
    # Float32 never gets past DAE initialization, whatever the layout or the init algorithm.
    (; name="Float32 dense Rodas5P",     T=Float32, layout=:dense, solver=Rodas5P,  broken=true),
    (; name="Float32 CSR Rodas5P",       T=Float32, layout=:csr,   solver=Rodas5P,  broken=true),
]

@testset "actual GPU solve" begin
    nw, s0 = Lib.powergridlike_network(; execution=KAExecution{true}(), aggregator=SparseAggregator(+))
    s0.v[3, :Pset] = -1.2 # force some dynamics
    s0.v[5, :u_i] = -1 # force reinit

    # host reference: initialization actually has to do something
    prob = ODEProblem(nw, s0, (0.0, 10.0))
    @test prob.f.jac_prototype == nothing
    sol = solve(prob, Rodas5P())
    @test sol(0.0, idxs=VIndex(5, :u_i)) != -1 # reinit happend

    nw_sparse = set_jac_prototype!(copy(nw))
    @test ODEProblem(nw_sparse, s0, (0.0, 1.0)).f.jac_prototype isa SparseMatrixCSC
    @test adapt(CuArray{Float64}, nw_sparse).jac_prototype isa CuSparseMatrixCSR

    for cfg in gpu_configs
        kwargs = Base.structdiff(cfg, (; name=nothing, broken=nothing))
        @testset "$(cfg.name)" begin
            if cfg.broken
                @test_broken run_on_gpu(nw, s0; kwargs...)
            else
                @test run_on_gpu(nw, s0; kwargs...)
            end
        end
    end
end


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
