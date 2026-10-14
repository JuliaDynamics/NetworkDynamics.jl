#=
At which network size does a GPU solve pay off?

Run from the package root against the test environment, which has CUDA, CUDSS and the
component library:

    julia --project=test --threads=8 benchmark/gpu_scaling.jl

The network is a square lattice of swing machines and PQ loads coupled by static dq lines, so a
mass-matrix DAE in which half of the states are algebraic. The graph matters a lot here: the
sparse LU fills in moderately on a planar graph, while randomly rewired graphs (Watts-Strogatz)
fill in almost completely and make every direct solver slow. Each configuration solves it
from a flat start, which means the DAE initialization is part of every solve.

Besides the plain CPU and GPU runs there is a mixed one: the network stays on the CPU and only
the linear solves go to cuDSS. That works for models the GPU path cannot handle yet, such as
networks with callbacks.

Writes `gpu_scaling.png` and `gpu_scaling.txt` next to this script.
=#
using NetworkDynamics
using Graphs, SparseArrays, Printf
using CUDA, CUDA.CUSPARSE, Adapt
using CUDSS # gives LinearSolve a direct sparse solver on the device
using OrdinaryDiffEqRosenbrock, OrdinaryDiffEqBDF, OrdinaryDiffEqNonlinearSolve
using NonlinearSolve # DAE init polyalg
using NonlinearSolve: LinearSolve
using LinearAlgebra
using Chairmarks
using CairoMakie
isdefined(Main, :Lib) || include(joinpath(pkgdir(NetworkDynamics), "test", "ComponentLibrary.jl"))

# OrdinaryDiffEq sets up its solvers with `findall(!iszero, mass_matrix)`, which falls back to
# the generic O(n²) method for a Diagonal and takes minutes at 2e5 states. Patched here so the
# CPU runs measure the solve, see https://github.com/SciML/OrdinaryDiffEq.jl/pull/4711
Base.findall(::typeof(!iszero), D::Diagonal) = [CartesianIndex(i, i) for i in axes(D, 1) if !iszero(D.diag[i])]

# LU on the device for a problem that lives on the CPU. Each new W is copied over and
# refactorized, which reuses the symbolic analysis; each solve only moves b and x.
struct CudssOffloadLU <: LinearSolve.AbstractFactorization end
mutable struct CudssOffloadCache
    F::Any
    b::CuVector{Float64}
    x::CuVector{Float64}
end
function LinearSolve.init_cacheval(::CudssOffloadLU, A, b, u, Pl, Pr, maxiters::Int, abstol, reltol,
        verbose::Union{Bool,LinearSolve.LinearVerbosity}, assumptions::LinearSolve.OperatorAssumptions)
    CudssOffloadCache(nothing, CuVector{Float64}(undef, length(b)), CuVector{Float64}(undef, length(b)))
end
function LinearSolve.SciMLBase.solve!(cache::LinearSolve.LinearCache, alg::CudssOffloadLU; kwargs...)
    c = cache.cacheval
    if cache.isfresh
        A = CuSparseMatrixCSR(convert(AbstractMatrix, cache.A))
        c.F = isnothing(c.F) ? lu(A) : lu!(c.F, A)
        cache.isfresh = false
    end
    copyto!(c.b, cache.b)
    ldiv!(c.x, c.F, c.b)
    copyto!(cache.u, c.x)
    LinearSolve.SciMLBase.build_linear_solution(alg, cache.u, nothing, cache)
end

# `N` is rounded to the nearest square
function scaling_network(N; kwargs...)
    n = round(Int, sqrt(N))
    g = Graphs.grid((n, n))
    swing = Lib.dqbus_swing(; Pmech=1.0, V=1.0, M=1.0, D=1.0)
    load = Lib.dqbus_pq(; Pset=-1.0, Qset=0.0)
    vs = [isodd(i) ? swing : load for i in 1:nv(g)]
    Network(g, vs, Lib.dqline(R=0.01, X=0.1); kwargs...)
end

# every swing at θ = ω = 0, load voltages at their default 1 + 0im
function flat_start(nw)
    s0 = NWState(nw)
    uflat(s0)[isnan.(uflat(s0))] .= 0
    s0
end

# `linsolve=nothing` keeps the solver's default
configs = [
    (; name="CPU sequential", T=Float64, gpu=false, linsolve=nothing,
       execution=SequentialExecution{true}(), aggregator=SequentialAggregator(+)),
    (; name="CPU threaded",   T=Float64, gpu=false, linsolve=nothing,
       execution=ThreadedExecution{true}(), aggregator=ThreadedAggregator(+)),
    (; name="CPU threaded, cuDSS LU", T=Float64, gpu=false, linsolve=CudssOffloadLU(),
       execution=ThreadedExecution{true}(), aggregator=ThreadedAggregator(+)),
    (; name="GPU Float64",    T=Float64, gpu=true, linsolve=nothing,
       execution=KAExecution{true}(), aggregator=SparseAggregator(+)),
    (; name="GPU Float32",    T=Float32, gpu=true, linsolve=nothing,
       execution=KAExecution{true}(), aggregator=SparseAggregator(+)),
]
solvers = [(; name="Rodas5P", alg=Rodas5P), (; name="FBDF", alg=FBDF)]
Ns = [10, 18, 32, 56, 100, 178, 316] .^ 2
tend = 5.0

# Network, initial state and parameters for one configuration, moved to the device if needed.
# The sparsity pattern is the same for every configuration, so it is detected once per size.
function setup(cfg, N, jac)
    nw = set_jac_prototype!(scaling_network(N; cfg.execution, cfg.aggregator), jac)
    s0 = flat_start(nw)
    if cfg.gpu
        to = CuArray{cfg.T}
        nw, s0 = adapt(to, nw), adapt(to, s0)
    end
    nw, s0
end

# Wall time of `f()`, best of a few runs. Large sizes get a single run.
function besttime(f; budget=10.0)
    t = @elapsed f()
    n = 1
    while n < 3 && t * (n + 1) < budget
        t = min(t, @elapsed f())
        n += 1
    end
    t
end

# Each result goes to the console and to `io` right away, so a run that dies (out of memory,
# say) still leaves everything measured up to that point.
function record!(results, io, r)
    push!(results, r)
    line = @sprintf("%7d %7d  %-15s %-8s %12.6g s  %-8s  maxrss %5.2f GiB",
                    r.N, r.dim, r.config, r.solver, r.time, r.retcode, Sys.maxrss() / 2^30)
    for out in (stdout, io)
        println(out, line)
        flush(out)
    end
end

function run_scaling(configs, solvers, Ns; io=devnull)
    results = []
    for N in Ns
        jac = get_jac_prototype(scaling_network(N); verbose=false)
        for cfg in configs
            nw, s0 = setup(cfg, N, jac)
            u, p = uflat(s0), pflat(s0)
            du = similar(u)
            sync = cfg.gpu ? CUDA.synchronize : () -> nothing
            # the offload config shares its rhs with "CPU threaded"
            if isnothing(cfg.linsolve)
                rhs = @b (nw(du, u, p, 0.0); sync())
                record!(results, io, (; N, dim=dim(nw), config=cfg.name, solver="rhs", time=rhs.time, retcode=:none))
            end

            T = cfg.T
            prob = ODEProblem(nw, s0, (zero(T), T(tend)))
            for s in solvers
                alg = isnothing(cfg.linsolve) ? s.alg() : s.alg(; cfg.linsolve)
                local sol
                time = besttime() do
                    sol = solve(prob, alg)
                    sync()
                end
                record!(results, io, (; N, dim=dim(nw), config=cfg.name, solver=s.name, time, retcode=sol.retcode))
            end
            GC.gc(); cfg.gpu && CUDA.reclaim()
        end
    end
    results
end

function plot_scaling(results, configs, solvers)
    fig = Figure(size=(1200, 450))
    title = "$(Sys.cpu_info()[1].model) ($(Threads.nthreads()) threads) vs $(CUDA.name(CUDA.device()))"
    Label(fig[0, 1:3], title)
    panels = [("full solve, t ∈ [0, $tend]", [s.name for s in solvers]), ("one rhs call", ["rhs"])]
    colors = Makie.wong_colors()
    styles = [:solid, :dash]
    for (col, (paneltitle, rows)) in enumerate(panels)
        ax = Axis(fig[1, col]; title=paneltitle, xlabel="number of vertices", ylabel="time [s]",
                  xscale=log10, yscale=log10)
        for (ci, cfg) in enumerate(configs), (si, row) in enumerate(rows)
            rs = filter(r -> r.config == cfg.name && r.solver == row, results)
            isempty(rs) && continue
            scatterlines!(ax, [r.N for r in rs], [r.time for r in rs]; color=colors[ci], linestyle=styles[si])
        end
    end
    # color is the configuration, line style the solver
    config_entries = [LineElement(; color=colors[i]) for i in eachindex(configs)]
    solver_entries = [LineElement(; color=:gray, linestyle=styles[i]) for i in eachindex(solvers)]
    Legend(fig[1, 3], [config_entries, solver_entries],
           [[c.name for c in configs], [s.name for s in solvers]], ["configuration", "solver"];
           framevisible=false)
    fig
end

function main()
    # warm up every code path once, so the first size does not carry the compile time
    run_scaling(configs, solvers, [100])
    results = open(joinpath(@__DIR__, "gpu_scaling.txt"), "w") do io
        run_scaling(configs, solvers, Ns; io)
    end
    save(joinpath(@__DIR__, "gpu_scaling.png"), plot_scaling(results, configs, solvers))
    results
end

# `include` in a REPL only defines everything, call `main()` or the pieces by hand
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
