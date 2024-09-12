is_precompiling() = ccall(:jl_generating_output, Cint, ()) == 1

using TimerOutputs

to = TimerOutput();
@timeit to "precompile workload" begin
    @timeit to "imports" begin
        using NetworkDynamics
        using Graphs
    end

    @timeit to "create graph" begin
        g = complete_graph(3)
    end
    @timeit to "load library" begin
        include(joinpath(@__DIR__,"..","test","ComponentLibrary.jl"))
    end
    @timeit to "create network" begin
        edge = Lib.diffusion_edge()
        vert = Lib.diffusion_vertex()
        nd = Network(g, vert, edge)
    end
    u = rand(dim(nd))
    du = similar(u)
    p = rand(pdim(nd))
    @timeit to "coreloop" begin
        nd(du,u,p,0.0)
    end
end
is_precompiling() || print_timer(to)

# without precompilation
#  ────────────────────────────────────────────────────────────────────────────────
#                                         Time                    Allocations
#                                ───────────────────────   ────────────────────────
#        Tot / % measured:            8.61s /  96.8%            848MiB /  97.3%

#  Section               ncalls     time    %tot     avg     alloc    %tot      avg
#  ────────────────────────────────────────────────────────────────────────────────
#  precompile workload        1    8.33s  100.0%   8.33s    825MiB  100.0%   825MiB
#    create network           1    4.80s   57.5%   4.80s    535MiB   64.9%   535MiB
#    imports                  1    3.12s   37.5%   3.12s    226MiB   27.4%   226MiB
#    coreloop                 1    328ms    3.9%   328ms   58.3MiB    7.1%  58.3MiB
#    create graph             1   52.8ms    0.6%  52.8ms   4.82MiB    0.6%  4.82MiB
#    load library             1   14.9ms    0.2%  14.9ms    367KiB    0.0%   367KiB
#  ────────────────────────────────────────────────────────────────────────────────
