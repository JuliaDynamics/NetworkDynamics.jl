using NetworkDynamics
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as Dt
using StableRNGs

@mtkmodel Swing begin
    @variables begin
        θ(t) = 0.0, [description = "voltage angle", output=true]
        ω(t) = 0.0, [description = "Rotor frequency"]
        Pel(t), [description = "Electical power flowing from network into node", input=true]
        Pdamping(t), [description = "Damping Power (observed)"]
    end
    @parameters begin
        M = 1, [description = "Inertia"]
        D = 0.1, [description = "Damping"]
        Pmech=1.0, [description = "Mechanical Power"]
    end
    @equations begin
        Dt(θ) ~ ω
        Pdamping ~ - D * ω
        Dt(ω) ~ 1/M * (Pmech + Pdamping + Pel)
    end
end

@mtkmodel Load begin
    @parameters begin
        Pd = -1.0, [description = "Power Load"]
    end
    @variables begin
        θ(t) = 0.0, [description = "voltage angle", output=true]
        Pel(t), [description = "Electrical power flowing from network into node", input=true]
    end
    @equations begin
        0 ~ Pel + Pd
    end
end

@mtkmodel StaticPowerLine begin
    @parameters begin
        K = 1.63, [description = "line conductance parameter"]
    end
    @variables begin
        srcθ(t), [description = "voltage angle at src end", input=true]
        dstθ(t), [description = "voltage angle at dst end", input=true]
        P(t), [description = "flow towards node at dst end", output=true]
        Δθ(t), [description = "Voltage angle difference"]
    end
    @equations begin
        Δθ ~ srcθ - dstθ
        P ~ K*sin(Δθ)
    end
end

function CombinedModel(subsystems...; kwargs...)
    @variables begin
        θ(t) = 0.0, [description = "voltage angle", output=true]
        Pel(t), [description = "Electical power flowing from network into node", input=true]
    end
    eqs = [
        (θ ~ sub.θ for sub in subsystems)...,
        Pel ~ +((sub.Pel for sub in subsystems)...) # kirchoff constaint
    ]
    ODESystem(eqs, t, [θ, Pel], []; systems=collect(subsystems), kwargs...)
end

v1 = let
    @named swing = Swing(; Pmech=1.0, M=1.0, D=1.0)
    VertexModel(swing, [:Pel], [:θ], vidx=1)
end

v2 = let
    @named load = Load(; Pd=-1.0)
    VertexModel(load, [:Pel], [:θ], vidx=2)
end

v3 = let
    @named swing = Swing(; Pmech=1.0, M=1.0, D=1.0)
    @named load = Load(; Pd=-2.0)
    @named swingload = CombinedModel(swing, load)
    VertexModel(swingload, [:Pel], [:θ], vidx=3)
end

v4 = let
    @named swing = Swing(; Pmech=2.0, M=1.0, D=1.0)
    @named load = Load(; Pd=-1.0)
    @named swingload2 = CombinedModel(swing, load)
    VertexModel(swingload2, [:Pel], [:θ], vidx=4)
end

e1 = let
    @named line = StaticPowerLine(; K=0.5)
    EdgeModel(line, [:srcθ], [:dstθ], [:P]; annotation=:AntiSymmetric, src=1, dst=2)
end
e2 = let
    @named line = StaticPowerLine(; K=1.0)
    EdgeModel(line, [:srcθ], [:dstθ], [:P]; annotation=:AntiSymmetric, src=2, dst=3)
end
e3 = let
    @named line = StaticPowerLine()
    EdgeModel(line, [:srcθ], [:dstθ], [:P]; annotation=:AntiSymmetric, src=3, dst=4)
end

nw = Network([v1,v2,v3,v4],[e1,e2,e3])

rng = StableRNG(1)
u = randn(rng, dim(nw))
du1 = [NaN for _ in 1:dim(nw)]
nw(du1, u, pflat(NWState(nw)), NaN)

data = NetworkDynamics.parse_network("testgrid.yaml")
data["VertexModels"][1]

constructors = Dict("Swing"=>Swing, "VertexModel"=>VertexModel, "Load"=>Load)
NetworkDynamics.build_network(data, constructors)
