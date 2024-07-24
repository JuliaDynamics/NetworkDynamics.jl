import ModelingToolkit as MTK
using ModelingToolkit
using ModelingToolkitStandardLibrary
using ModelingToolkit, OrdinaryDiffEq, Plots
using ModelingToolkitStandardLibrary.Electrical
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkit: t_nounits as t, D_nounits as Dt
using NetworkDynamics

####
#### Try again with specific definition of input/output
####
@mtkmodel Bus begin
    @variables begin
        θ(t), [description = "voltage angle", output=true]
        P(t), [description = "Electical Powerflow into Network", input=true]
    end
end


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
end


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
    @equations begin
        srcP ~ -K*sin(dstθ - srcθ)
        dstP ~ -srcP
    end
end
@named line = StaticPowerLine()

@named swing = SwingNode(Pmech=1)

v = ODEVertex(swing, [:P], [:θ])
e = StaticEdge(line, [:dstθ, :srcθ], [:dstP, :srcP])

unknowns(line)






hasmethod

full_equations(swing)
equations(swing)
_sys = swing
inputs = MTK.unbound_inputs(swing)
outputs = MTK.unbound_outputs(swing)

swing_simp,inidx = structural_simplify(swing, (MTK.unbound_inputs(swing), MTK.unbound_outputs(swing)))
sys = swing
sys.P
MTK.getvar

inidx
inidx[1]

swing_simp[inidx]

full_equations(swing)
full_equations(swing_simp)
MTK.unbound_inputs(swing_simp)



MTK.calculate_massmatrix(swing)

@which show(swing_simp

state = MTK.get_tearing_state(swing_simp)
state.fullvars
state.structure.graph


full_equations(swing)
equation_dependencies(swing)
variable_dependencies(swing)
unknowns(swing)

@which SystemStructure
ODESystem

collect(edges(variable_dependencies(swing)))
unknowns(swing)


full_equations(line)
equation_dependencies(line)




MTK.generate_control_function(swing)

swing_simp,inidx = structural_simplify(swing, (MTK.unbound_inputs(swing), []))
equations(swing_simp)
full_equations(swing_simp)
observed(swing_simp)

MTK.inputs(swing)
MTK.generate_control_function(swing, MTK.inputs(swing))

MTK.is_bound

#=
- XXX: input cannot be declared on subsystem level as it leads to isbound_true on the system level
- works with extend though
=#
