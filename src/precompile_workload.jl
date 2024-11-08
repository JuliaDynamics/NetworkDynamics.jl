is_precompiling() = ccall(:jl_generating_output, Cint, ()) == 1

using TimerOutputs: TimerOutput, @timeit, print_timer
using NetworkDynamics
using Graphs: complete_graph

Base.@propagate_inbounds function diffusionedge!(e, v_s, v_d, (p,), _)
    e .= p * (v_s[1] .- v_d[1])
end

Base.@propagate_inbounds function diffusionvertex!(dv, _, acc, _, _)
    dv[1] = acc[1]
    nothing
end

to = TimerOutput();
@timeit to "precompile workload" begin
    @timeit to "create graph" begin
        g = complete_graph(3)
    end
    @timeit to "create network" begin
        edge = EdgeFunction(g=AntiSymmetric(diffusionedge!), outdim=1, pdim=1, pdef=[1])
        vert = VertexFunction(f=diffusionvertex!, g=1, dim=1, pdim=0)
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
