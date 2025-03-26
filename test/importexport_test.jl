using NetworkDynamics
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as Dt

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
    @structural_parameters begin
        K = 1.63 # structural parameter to keep them out of nw parameters
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
    θ_eqs   = [θ ~ sub.θ for sub in subsystes]
    Pel_eqs = [Pel ~ sub.Pel for sub in subsystems]
    ODESystem(; kwargs...)
end

mtk_models = [
    Swing(Pmech=1, M=1, D=1.0; name=:G1),
    Load(Pload; name=:L2),
    Swing(Pmech=1, M=1, D=1.0; name=:G2),
    Load(Pload; name=:L3),
    Swing(Pmech=1, M=1, D=1.0; name=:G3),
    Load(Pload; name=:L4),
];

vertexmodels = map(mtk_models) do mtk_model
    VertexModel(mtk_model, [:Pel], [:θ])
end;
vertexmodels2 = map(mtk_models2) do mtk_model
    VertexModel(mtk_model, [:Pel], [:θ])
end;


data = NetworkDynamics.parse_network("testgrid.yaml")

constructors = Dict("Swing"=>Swing, "VertexModel"=>VertexModel, "Load"=>Load)
NetworkDynamics.build_network(data, constructors)
