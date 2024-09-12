is_precompiling() = ccall(:jl_generating_output, Cint, ()) == 1

using NetworkDynamics
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as Dt
using Graphs

@mtkmodel Bus begin
    @variables begin
        θ(t), [description = "voltage angle", output=true]
        P(t), [description = "Electical Powerflow into Network", input=true]
    end
end;

@mtkmodel SwingNode begin
    @extend Bus()
    @variables begin
        ω(t) = 0.0, [description = "Rotor frequency"]
    end
    @parameters begin
        M = 1, [description = "Inertia"]
        D = 0.1, [description = "Damping"]
        Pmech, [description = "Mechanical Power"]
    end
    @equations begin
        Dt(θ) ~ ω
        Dt(ω) ~ 1/M * (Pmech - D*ω + P)
    end
end;

@mtkmodel Line begin
    @variables begin
        srcθ(t), [description = "voltage angle at src end", input=true]
        srcP(t), [description = "flow towards not at src end", output=true]
        dstθ(t), [description = "voltage angle at dst end", input=true]
        dstP(t), [description = "flow towards not at dst end", output=true]
    end
end
@mtkmodel StaticPowerLine begin
    @extend Line()
    @parameters begin
        K = 100.0, [description = "Line conductance"]
    end
    @variables begin
        Δθ(t)
    end
    @equations begin
        Δθ ~ dstθ - srcθ
        srcP ~ -K*sin(Δθ)
        dstP ~ -srcP
    end
end

@named swing = SwingNode()
@named line = StaticPowerLine()

if is_precompiling()
    ODEVertex(swing, [:P], [:θ])
    StaticEdge(line, [:srcθ], [:dstθ], [:dstP, :srcP], Fiducial())
else
    @info "ODEVertex"
    @time @eval ODEVertex(swing, [:P], [:θ])
    @info "StaticEdge"
    @time @eval StaticEdge(line, [:srcθ], [:dstθ], [:dstP, :srcP], Fiducial())
end

# without precompile statements
# julia> include("/home/hw/.julia/dev/NetworkDynamics/ext/precompile_workload.jl");
# [ Info: ODEVertex
#  12.337674 seconds (18.87 M allocations: 960.629 MiB, 2.49% gc time, 99.61% compilation time: 5% of which was recompilation)
# [ Info: StaticEdge
#   7.627367 seconds (6.62 M allocations: 337.839 MiB, 0.64% gc time, 99.65% compilation time)

# with precompile statements
# julia> include("/home/hw/.julia/dev/NetworkDynamics/ext/precompile_workload.jl");
# [ Info: ODEVertex
#   2.570151 seconds (3.11 M allocations: 157.414 MiB, 1.24% gc time, 98.34% compilation time: 86% of which was recompilation)
# [ Info: StaticEdge
#   0.201824 seconds (23.27 k allocations: 1.269 MiB, 10.63% gc time, 92.13% compilation time)
