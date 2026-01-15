using Test
using NetworkDynamics
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using OrdinaryDiffEqTsit5
using Graphs

# Shared components used by both implementations

@mtkmodel NWTerminal begin
    @variables begin
        v(t), [description="Voltage at node"]
        i(t), [description="Current flowing into node"]
    end
end

@mtkmodel VoltageSource begin
    @components begin
       p = NWTerminal()
    end
    @parameters begin
        V = 1.0
    end
    @variables begin
        i(t), [description="produced current by ideal voltage source (observable)"]
    end
    @equations begin
        i ~ -p.i
        p.v ~ V
    end
end
@named vs = VoltageSource()
vs_vertex = VertexModel(vs, [:p₊i], [:p₊v])

@mtkmodel Resistor begin
    @components begin
        src = NWTerminal()
        dst = NWTerminal()
    end
    @parameters begin
        R = 1
    end
    @equations begin
        dst.i ~ (src.v - dst.v)/R
        src.i ~ -dst.i
    end
end
@named resistor = Resistor()

# Part A: Modeling with Injector Nodes and LoopbackConnection

@mtkmodel Capacitor begin
    @components begin
        p = NWTerminal(;v=0)
    end
    @parameters begin
        C = 1.0
    end
    @equations begin
        D(p.v) ~ p.i / C
    end
end
@named cap = Capacitor()
hub_vertex = VertexModel(cap, [:p₊i], [:p₊v], name=:hub)

@mtkmodel ResistorInjector begin
    @components begin
        p = NWTerminal()
    end
    @parameters begin
        R = 100.0
    end
    @equations begin
        p.i ~ p.v/R
    end
end
@named resistor_inj = ResistorInjector()
Rinj_vertex = VertexModel(resistor_inj, [:p₊v], [:p₊i], name=:R_injector, ff_to_constraint=false)

@mtkmodel InductorInjector begin
    @components begin
        p = NWTerminal()
    end
    @parameters begin
        L = 0.1
    end
    @equations begin
        D(p.i) ~ p.v / L
    end
end
@named inductor_inj = InductorInjector()
Linj_vertex = VertexModel(inductor_inj, [:p₊v], [:p₊i], name=:L_injector, ff_to_constraint=false)

edges1 = [
    LoopbackConnection(potential=[:u], flow=[:i], src=:R_injector, dst=:hub, name=:R_loopback),
    LoopbackConnection(potential=[:u], flow=[:i], src=:L_injector, dst=:hub, name=:L_loopback),
    EdgeModel(resistor, [:src₊v], [:dst₊v], [:src₊i], [:dst₊i]; src=:vs, dst=:hub)
]
vertices1 = [vs_vertex, hub_vertex, Rinj_vertex, Linj_vertex]
nw1 = Network(vertices1, edges1; warn_order=false)

s0_1 = NWState(nw1)
s0_1.v[:hub, :p₊v] = 0.0
s0_1.v[:L_injector, :p₊i] = 0.0

prob1 = ODEProblem(nw1, s0_1, (0.0, 10.0))
Main.test_execution_styles(prob1) # testing all ex styles #src
sol1 = solve(prob1, Tsit5(), saveat=0.1)
@test SciMLBase.successful_retcode(sol1)

# Part B: Modeling with a Single VertexModel

@mtkmodel CRLModel begin
    @components begin
        cap = Capacitor()
        resistor = ResistorInjector()
        inductor = InductorInjector()
    end
    @variables begin
        v(t), [description="Voltage at node"]
        i(t)=0, [description="Current flowing into node"]
    end
    @equations begin
        0 ~ resistor.p.i + inductor.p.i + cap.p.i - i
        v ~ cap.p.v
        v ~ resistor.p.v
        v ~ inductor.p.v
    end
end
@named crl_model = CRLModel()
crl_vertex = VertexModel(crl_model, [:i], [:v], name=:CRL_vertex)

edges2 = [
    EdgeModel(resistor, [:src₊v], [:dst₊v], [:src₊i], [:dst₊i]; src=:vs, dst=:CRL_vertex)
]
vertices2 = [vs_vertex, crl_vertex]
nw2 = Network(vertices2, edges2; warn_order=false)

s0_2 = NWState(nw2)
s0_2[VIndex(:CRL_vertex, :v)] = 0.0
s0_2[VIndex(:CRL_vertex, :inductor₊p₊i)] = 0.0

prob2 = ODEProblem(nw2, s0_2, (0.0, 10.0))
Main.test_execution_styles(prob2) # testing all ex styles #src
sol2 = solve(prob2, Tsit5(), saveat=0.1)
@test SciMLBase.successful_retcode(sol2)

@test sol1[VIndex(2,:p₊v)] ≈ sol2[VIndex(2, :v)]
@test sol1[VIndex(4,:p₊i)] ≈ sol2[VIndex(2, :inductor₊p₊i)]
