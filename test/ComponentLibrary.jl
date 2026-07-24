module Lib
using Graphs
using NetworkDynamics
using ModelingToolkitBase
using ModelingToolkitBase: t_nounits as t, D_nounits as Dt
using SciCompDSL

Base.@propagate_inbounds function diffusionedge!(e, v_s, v_d, (p,), _)
    e .= p * (v_s[1] .- v_d[1])
end
function diffusion_edge()
    EdgeModel(g=AntiSymmetric(diffusionedge!), outdim=1, pdim=1, pdef=[1], name=:diff_edge)
end
function diffusion_edge_closure()
    force_closure = rand()
    gs = (e, v_s, v_d, (p,), _) -> begin
        e .= p * (v_s[1] .- v_d[1]) .+ 0 * force_closure
    end
    EdgeModel(;g=AntiSymmetric(gs), outdim=1, pdim=1, pdef=[1])
end

Base.@propagate_inbounds function diffusionedge_fid!(e_s, e_d, v_s, v_d, (p,), _)
    e_d[1] = p * (v_s[1] .- v_d[1])
    e_s[1] = -e_d[1]
end
function diffusion_edge_fid()
    EdgeModel(g=diffusionedge_fid!, outdim=1, pdim=1, pdef=[1], name=:diff_edge_fid)
end

Base.@propagate_inbounds function diffusion_dedge!(de, e, v_s, v_d, (τ,), _)
    de[1] = 1 / τ * (sin(v_s[1] - v_d[1]) - e[1])
    de[2] = 1 / τ * (sin(v_d[1] - v_s[1]) - e[2])
    nothing
end
function diffusion_odeedge()
    EdgeModel(f=diffusion_dedge!,
        dim=2, sym=[:e_dst, :e_src],
        pdim=1, pdef=[100], psym=[:τ],
        g=Fiducial(dst=1:1, src=2:2), name=:diff_edge_ode)
end

Base.@propagate_inbounds function diffusionvertex!(dv, _, acc, _, _)
    dv[1] = acc[1]
    nothing
end
diffusion_vertex() = VertexModel(f=diffusionvertex!, dim=1, g=1:1)

####
#### inhomogenious kuramoto system
####
Base.@propagate_inbounds function kuramoto_edge!(e, θ_s, θ_d, (K,), t)
    e .= K .* sin(θ_s[1] - θ_d[1])
end
function kuramoto_edge(; name=:kuramoto_edge, kwargs...)
    EdgeModel(;g=AntiSymmetric(kuramoto_edge!),
        outsym=[:P], psym=[:K], name, kwargs...)
end

Base.@propagate_inbounds function kuramoto_inertia!(dv, v, acc, p, t)
    M, D, Pm = p
    dv[1] = v[2]
    dv[2] = 1 / M * (Pm - D * v[2] + acc[1])
end
function kuramoto_second(; name=:kuramoto_second, kwargs...)
    VertexModel(; f=kuramoto_inertia!, sym=[:δ=>0, :ω=>0],
        psym=[:M=>1, :D=>0.1, :Pm=>1], g=StateMask(1), name, kwargs...)
end

Base.@propagate_inbounds function kuramoto_vertex!(dθ, θ, esum, (ω,), t)
    dθ[1] = ω + esum[1]
end
kuramoto_first() = VertexModel(; f=kuramoto_vertex!, sym=[:θ], psym=[:ω], g=1)

function swing_mtk(; kwargs...)
    @mtkmodel SwingNode begin
        @variables begin
            θ(t) = 0.0, [description = "voltage angle", output=true]
            P(t), [description = "Electical Powerflow into Network", input=true]
            ω(t) = 0.0, [description = "Rotor frequency"]
            Pdamping(t), [description = "Damping Power (observed)"]
        end
        @parameters begin
            M = 1, [description = "Inertia"]
            D = 0.1, [description = "Damping"]
            Pmech, [description = "Mechanical Power"]
        end
        @equations begin
            Dt(θ) ~ ω
            Pdamping ~ - D * ω
            Dt(ω) ~ 1/M * (Pmech + Pdamping + P)
        end
    end
    VertexModel(SwingNode(name=:swing_mtk), [:P], [:θ]; kwargs...)
end

function line_mtk(; kwargs...)
    @mtkmodel StaticPowerLine begin
        @variables begin
            srcθ(t), [description = "voltage angle at src end", input=true]
            dstθ(t), [description = "voltage angle at dst end", input=true]
            P(t), [description = "flow towards node at dst end", output=true]
            Δθ(t)
        end
        @parameters begin
            active = 1, [description = "line active"]
            K = 1.63, [description = "Line conductance"]
            limit = 1, [description = "Line limit"]
        end
        @equations begin
            Δθ ~ srcθ - dstθ
            P ~ active*K*sin(Δθ)
        end
    end
    EdgeModel(StaticPowerLine(name=:line_mtk), [:srcθ], [:dstθ], AntiSymmetric([:P]); kwargs...)
end

####
#### Powergrid Models
####
@mtkmodel BusBase begin
    @variables begin
        u_r(t)=1.0, [description = "voltage real part", output=true]
        u_i(t)=0.0, [description = "voltage imaginary part", output=true]
        i_r(t), [description = "current real part (coming from network)", input=true, guess=1.0]
        i_i(t), [description = "current imaginary part (coming from network)", input=true, guess=0.0]
        θ_meas(t), [description = "voltage angle"]
        u_mag(t), [description = "voltage magnitude"]
        Pinj(t), [description = "Active Power injection into network"]
        Qinj(t), [description = "Reactive Power injection into network"]
    end
    @equations begin
        θ_meas ~ atan(u_i, u_r)
        u_mag ~ sqrt(u_r^2 + u_i^2)
        Pinj ~ - u_r*i_r - u_i*i_i
        Qinj ~ - u_i*i_r + u_r*i_i
    end
end

@mtkmodel SwingDQ begin
    @extend BusBase()
    @variables begin
        θ(t), [description = "voltage angle", guess=0.0]
        ω(t), [description = "Rotor frequency", guess=0.0]
        Pel(t), [description = "Electical power flowing from network into node"]
        Pdamping(t), [description = "Damping Power (observed)"]
    end
    @parameters begin
        M = 1, [description = "Inertia"]
        D = 0.1, [description = "Damping"]
        Pmech, [description = "Mechanical Power", guess=1.0]
        V, [description = "Voltage magnitude", guess=1.0]
    end
    @equations begin
        Pel ~ u_r*i_r + u_i*i_i
        Pdamping ~ - D * ω
        Dt(θ) ~ ω
        Dt(ω) ~ 1/M * (Pmech + Pdamping + Pel)
        u_r ~ V*cos(θ)
        u_i ~ V*sin(θ)
    end
end

_no_simplify(x) = x
@register_symbolic _no_simplify(x)

@mtkmodel PQLoad begin
    @extend BusBase()
    @parameters begin
        Pset, [description = "Active Power setpoint", guess=0]
        Qset, [description = "Reactive Power setpoint", guess=0]
    end
    @equations begin
        0 ~ u_r*i_r + u_i*i_i + Pset
        0 ~ u_r*i_i - u_i*i_r + Qset
    end
end

@mtkmodel TimeDependentPQLoad begin
    @extend BusBase()
    @structural_parameters begin
        Pfun
    end
    @parameters begin
        Qset, [description = "Reactive Power setpoint", guess=0]
    end
    @equations begin
        0 ~ u_r*i_r + u_i*i_i + Pfun(t)
        0 ~ u_r*i_i - u_i*i_r + Qset
    end
end

@mtkmodel PVBus begin
    @extend BusBase()
    @parameters begin
        Pset, [description = "Active Power setpoint"]
        Vset = 1.0, [description = "Voltage setpoint"]
    end
    @equations begin
        0 ~ u_r*i_r + u_i*i_i + Pset
        Vset^2 ~ u_r^2 + u_i^2
    end
end

@mtkmodel SlackBus begin
    @extend BusBase()
    @equations begin
        u_r ~ 1
        u_i ~ 0
    end
end

@mtkmodel StaticPowerLineDQ begin
    @variables begin
        src_u_r(t), [description = "voltage real part", input=true]
        src_u_i(t), [description = "voltage imaginary part", input=true]
        src_i_r(t), [description = "current real part", output=true]
        src_i_i(t), [description = "current imaginary part", output=true]
        dst_u_r(t), [description = "voltage real part", input=true]
        dst_u_i(t), [description = "voltage imaginary part", input=true]
        dst_i_r(t), [description = "current real part", output=true]
        dst_i_i(t), [description = "current imaginary part", output=true]
        src_P(t), [description = "active power at src end"]
        dst_P(t), [description = "active power at dst end"]
    end
    @parameters begin
        R = 0.0, [description = "line resistance"]
        X = 1.0, [description = "line reactance"]
        active = 1, [description = "line active", guess=1]
    end
    begin
        Z = R + im*X
        Vsrc = src_u_r + im*src_u_i
        Vdst = dst_u_r + im*dst_u_i
        idst = active * 1/Z * (Vsrc - Vdst)
        isrc = -idst
    end
    @equations begin
        src_i_r ~ simplify(real(isrc))
        src_i_i ~ simplify(imag(isrc))
        dst_i_r ~ simplify(real(idst))
        dst_i_i ~ simplify(imag(idst))
        src_P ~ simplify(real(Vsrc * conj(isrc)))
        dst_P ~ simplify(real(Vdst * conj(idst)))
    end
end

@mtkmodel SwingAndLoadDQ begin
    @extend BusBase()
    @components begin
        swing = SwingDQ()
        load = PQLoad()
    end
    @equations begin
        swing.u_r ~ u_r
        swing.u_i ~ u_i
        load.u_r ~ u_r
        load.u_i ~ u_i
        0 ~ i_r - swing.i_r - load.i_r
        0 ~ i_i - swing.i_i - load.i_i
    end
end

function dqbus_swing(; kwargs...)
    @named swing = SwingDQ(; kwargs...)
    VertexModel(swing, [:i_r, :i_i], [:u_r, :u_i])
end
function dqbus_pq(; name=:pq, kwargs...)
    pq = PQLoad(; name, kwargs...)
    VertexModel(pq, [:i_r, :i_i], [:u_r, :u_i]; verbose=false)
end
function dqbus_timedeppq(; kwargs...)
    @named pq_timedep = TimeDependentPQLoad(; kwargs...)
    VertexModel(pq_timedep, [:i_r, :i_i], [:u_r, :u_i])
end
function dqbus_pv(; injector=false, name=:pv, kwargs...)
    pv = PVBus(; name, kwargs...)
    if injector
        vm = VertexModel(pv, [:u_r, :u_i], [:i_r, :i_i]; ff_to_constraint=false, assume_io_coupling=true)
        set_default!(vm, :i_r, 0.0)
        set_default!(vm, :i_i, 0.0)
        vm
    else
        VertexModel(pv, [:i_r, :i_i], [:u_r, :u_i]; verbose=false)
    end
end
function dqbus_slack(; injector=false, kwargs...)
    @named pv = SlackBus(; kwargs...)
    if injector
        vm = VertexModel(pv, [:u_r, :u_i], [:i_r, :i_i]; ff_to_constraint=false, assume_io_coupling=true)
        set_default!(vm, :i_r, 0.0)
        set_default!(vm, :i_i, 0.0)
        vm
    else
        VertexModel(pv, [:i_r, :i_i], [:u_r, :u_i])
    end
end
function dqbus_swing_and_load(; kwargs...)
    @named swing_and_load = SwingAndLoadDQ(; kwargs...)
    VertexModel(swing_and_load, [:i_r, :i_i], [:u_r, :u_i]; verbose=false)
end
function dqline(; name=:line, R, X, kwargs...)
    line = StaticPowerLineDQ(; name, R, X)
    EdgeModel(
        line,
        [:src_u_r, :src_u_i], [:dst_u_r, :dst_u_i],
        [:src_i_r, :src_i_i], [:dst_i_r, :dst_i_i];
        kwargs...
    )
end

function powergridlike_network()
    # Create a simple network with Kuramoto oscillators
    g = cycle_graph(5) # 5-node cycle graph
    v1s = dqbus_slack()
    v2s = dqbus_pv(Pset=1.5, Vset=1.0)
    v3s = dqbus_pq(Pset=-1.0, Qset=-0.1)
    v4s = dqbus_pq(Pset=-1.0, Qset=-0.1)
    v5s = dqbus_pq(Pset=-1.0, Qset=-0.1)
    e = dqline(X=0.1, R=0.01)

    nws = Network(g, [v1s, v2s, v3s, v4s, v5s], e)
    pf = find_fixpoint(nws)

    v1 = dqbus_swing_and_load()
    set_initconstraint!(v1, @initconstraint begin
        :load₊Pinj + 1.0
        - :u_r^2 - :u_i^2 + :swing₊V^2
    end)
    v2 = dqbus_swing()
    v3 = dqbus_pq()
    v4 = dqbus_pq()
    v5 = dqbus_pq()
    nw = Network(g, [v1, v2, v3, v4, v5], e; dealias=true)

    default_overrides = merge(
        interface_values(pf),
        Dict(EIndex(1:5, :active) .=> nothing))

    s0 = initialize_componentwise!(nw; subverbose=false, verbose=false, default_overrides)
    nw, s0
end

####
#### Injector node versions (reuse same MTK models, flip interface direction)
####
@mtkmodel ShuntHub begin
    @extend BusBase()
    @parameters begin
        G, [guess=1, description = "Shunt conductance"]
        B, [guess=1, description = "Shunt susceptance"]
    end
    begin
        Y = G + im*B
        ishunt = -Y * (u_r + im*u_i)
    end
    @equations begin
        i_r ~ simplify(real(ishunt))
        i_i ~ simplify(imag(ishunt))
    end
end

function dqbus_swing_injector(; kwargs...)
    @named swing = SwingDQ(; kwargs...)
    VertexModel(swing, [:u_r, :u_i], [:i_r, :i_i]; ff_to_constraint=false, assume_io_coupling=true)
end
function dqbus_pq_injector(; kwargs...)
    @named pq = PQLoad(; kwargs...)
    VertexModel(pq, [:u_r, :u_i], [:i_r, :i_i]; ff_to_constraint=false, assume_io_coupling=false)
end
function dqbus_shunt_hub(; kwargs...)
    @named shunt_hub = ShuntHub(; kwargs...)
    VertexModel(shunt_hub, [:i_r, :i_i], [:u_r, :u_i])
end

function powergridlike_injector_network()
    # Same 5-bus cycle topology as powergridlike_network, but buses 1 and 2
    # are split into hub + swing injector connected via LoopbackConnection.
    edges = [
        dqline(X=0.1, R=0.04, src=:hub1, dst=:hub2),
        dqline(X=0.15, R=0.1, src=:hub2, dst=:pq3),
        dqline(X=0.2, R=0.03, src=:pq3, dst=:pq4),
        dqline(X=0.02, R=0.015, src=:pq4, dst=:pq5),
        dqline(X=0.09, R=0.02, src=:pq5, dst=:hub1),
        LoopbackConnection(potential=[:u_r, :u_i], flow=[:i_r, :i_i], src=:inj1, dst=:hub1, name=:loopback1),
        LoopbackConnection(potential=[:u_r, :u_i], flow=[:i_r, :i_i], src=:inj2, dst=:hub2, name=:loopback2),
    ]

    # Powerflow: same topology with shunt hubs + injector PV/slack nodes
    pf_vertices = [
        dqbus_shunt_hub(G=0.1, B=0.01, name=:hub1),
        dqbus_shunt_hub(G=0.1, B=0.01, name=:hub2),
        dqbus_pq(Pset=-1.0, Qset=-0.1, name=:pq3),
        dqbus_pq(Pset=-1.0, Qset=-0.1, name=:pq4),
        dqbus_pq(Pset=-1.0, Qset=-0.1, name=:pq5),
        dqbus_slack(injector=true, name=:inj1),
        dqbus_pv(Pset=1.5, Vset=1.0, injector=true, name=:inj2),
    ]
    nws = Network(pf_vertices, edges)
    pf = find_fixpoint(nws)

    # Dynamic model: replace injectors with swing equations
    dyn_vertices = [
        dqbus_shunt_hub(name=:hub1),
        dqbus_shunt_hub(name=:hub2),
        dqbus_pq(name=:pq3),
        dqbus_pq(name=:pq4),
        dqbus_pq(name=:pq5),
        dqbus_swing_injector(name=:inj1),
        dqbus_swing_injector(name=:inj2),
    ]
    nw = Network(dyn_vertices, edges)

    default_overrides = interface_values(pf)
    s0 = initialize_componentwise!(nw; subverbose=false, verbose=false, default_overrides)
    nw, s0
end

####
#### Per-unit handling models (bound_to / weak / default_from)
####
# Anchored on the PowerDynamics use case: a bus carries a `busbar` sub-component with the base
# quantities `S_b`/`V_b` (the single source of truth), and devices on the bus must agree with
# that base. The physics is deliberate nonsense — only the structure is realistic — so the same
# models exercise the whole per-unit feature set:
#   (B) `bound_to`     — the injector's `S_b` is a hard alias of `busbar₊S_b` (implemented).
#   (A) `weak`         — a device rating that defaults to the base but stays settable (future).
#   (C) `default_from` — a line end picking up the base of its endpoint bus (future).
# The extension points for (A)/(C) are marked below so the models grow without being rewritten.

"The busbar sub-component: holds the base quantities and one trivial state that uses them so MTK
keeps the parameters. `busbar₊S_b` is the true parameter every device on the bus binds to."
function busbar_sys(; name)
    @parameters S_b=100.0 V_b=1.0
    @variables ub(t)=1.0
    System([Dt(ub) ~ S_b * V_b - ub], t; name)
end

"A bus vertex: a `busbar` plus an injector whose per-unit base `S_b` is `bound_to` the busbar's,
so it is eliminated from `psym` and can never desync from the bus base. The injector also carries
its own rating `Sn` which *weakly* defaults to the busbar base (`initf_weak`): it follows the base
unless the user pins it, the (A) `weak` use case — unlike `S_b`, `Sn` stays a live parameter."
function pu_bus(; name)
    @named busbar = busbar_sys()
    @variables u(t)=1.0 i(t) o(t)
    @parameters P=1.0
    @parameters S_b [bound_to = :busbar₊S_b]
    @parameters Sn [initf_weak = S_b]
    eqs = [Dt(u) ~ -u + i + P / S_b - u * Sn / busbar.S_b, o ~ u]
    System(eqs, t; name, systems=[busbar])
end

pu_bus_vertex(; kwargs...) = VertexModel(pu_bus(name=:bus), [:i], [:o]; kwargs...)

"A line edge connecting two buses. Physics is a trivial admittance flow."
function pu_line(; name)
    @variables src_u(t) dst_u(t) flow(t)
    @parameters Y=0.5
    # (C) `default_from` — the line-end base picked up from an endpoint bus, with an internal
    # quantity `bound_to` that line-end base — is exercised in `test/default_from_test.jl`.
    System([flow ~ Y * (src_u - dst_u)], t; name)
end

pu_line_edge(; kwargs...) = EdgeModel(pu_line(name=:line), [:src_u], [:dst_u], AntiSymmetric([:flow]); kwargs...)

"A minimal 2-bus network: two identical PU buses joined by one line."
function pu_network()
    v = pu_bus_vertex()
    Network(path_graph(2), [v, v], pu_line_edge())
end

####
#### Per-unit DQ injector network — all three features end-to-end
####
# A realistic-shaped grid slice (same hub/injector/line topology as
# `powergridlike_injector_network`) that exercises `bound_to`, `default_from` and the weak
# `initf_weak` together. The base rides on *observables only* — it never scales the dynamics — so
# it is inert to the physics yet still flows through every feature:
#   • hub — carries `busbar₊S_b` (the source of truth); a hub-local `S_b` is `bound_to` it.
#   • inj — its `S_b` is `default_from = (:hub, …)`; its rating `Sn` weakly follows `S_b`
#           (`initf_weak`) but stays independently settable.
#   • line — its end base is `default_from = (:src, …)`; an internal base is `bound_to` it.

"Stateless carrier of the system power base; kept alive by a per-unit observable on its host."
@mtkmodel PUBusBar begin
    @parameters begin
        S_b = 100.0, [description = "system power base [MVA]"]
    end
end

"A shunt hub carrying the busbar base. The hub-local `S_b` is `bound_to` `busbar₊S_b` (eliminated)."
@mtkmodel PUShuntHub begin
    @extend BusBase()
    @components begin
        busbar = PUBusBar()
    end
    @parameters begin
        G = 0.1, [description = "shunt conductance"]
        B = 0.01, [description = "shunt susceptance"]
        S_b, [bound_to = :busbar₊S_b, description = "hub base (alias of busbar₊S_b)"]
    end
    @variables begin
        S_shunt_pu(t), [description = "shunt apparent power in per-unit of the bus base"]
    end
    begin
        Y = G + im*B
        ishunt = -Y * (u_r + im*u_i)
    end
    @equations begin
        i_r ~ simplify(real(ishunt))
        i_i ~ simplify(imag(ishunt))
        S_shunt_pu ~ sqrt(Pinj^2 + Qinj^2) / S_b
    end
end

"A swing injector: base `default_from` the hub, rating `Sn` weakly follows the base (`initf_weak`)."
@mtkmodel PUSwingInjector begin
    @extend BusBase()
    @variables begin
        θ(t), [description = "voltage angle", guess = 0.0]
        ω(t), [description = "rotor frequency", guess = 0.0]
        Pel(t), [description = "electrical power from network into node"]
        loading(t), [description = "electrical loading in per-unit of the machine rating"]
    end
    @parameters begin
        M = 1, [description = "inertia"]
        D = 0.1, [description = "damping"]
        Pmech, [description = "mechanical power", guess = 1.0]
        V, [description = "voltage magnitude", guess = 1.0]
        S_b, [default_from = (:hub, :busbar₊S_b), description = "base copied weakly from the hub"]
        Sn, [initf_weak = S_b, description = "machine rating; weakly follows the base, stays settable"]
    end
    @equations begin
        Pel ~ u_r*i_r + u_i*i_i
        loading ~ Pel / Sn
        Dt(θ) ~ ω
        Dt(ω) ~ 1/M * (Pmech - D*ω + Pel)
        u_r ~ V*cos(θ)
        u_i ~ V*sin(θ)
    end
end

"A line whose end base is `default_from` its src bus, with an internal base `bound_to` that end base."
@mtkmodel PUStaticLine begin
    @variables begin
        src_u_r(t), [input=true]; src_u_i(t), [input=true]
        src_i_r(t), [output=true]; src_i_i(t), [output=true]
        dst_u_r(t), [input=true]; dst_u_i(t), [input=true]
        dst_i_r(t), [output=true]; dst_i_i(t), [output=true]
        src_P(t), [description = "active power at src end"]
        src_P_pu(t), [description = "src active power in per-unit of the line base"]
    end
    @parameters begin
        R = 0.04, [description = "line resistance"]
        X = 0.1, [description = "line reactance"]
        S_b, [default_from = (:src, :busbar₊S_b), description = "endpoint base copied from src bus"]
        S_b_int, [bound_to = :S_b, description = "internal base (alias of the endpoint base)"]
    end
    begin
        Z = R + im*X
        Vsrc = src_u_r + im*src_u_i
        Vdst = dst_u_r + im*dst_u_i
        idst = 1/Z * (Vsrc - Vdst)
        isrc = -idst
    end
    @equations begin
        src_i_r ~ simplify(real(isrc)); src_i_i ~ simplify(imag(isrc))
        dst_i_r ~ simplify(real(idst)); dst_i_i ~ simplify(imag(idst))
        src_P ~ simplify(real(Vsrc * conj(isrc)))
        src_P_pu ~ src_P / S_b_int
    end
end

pu_dqhub(; name, kwargs...) = VertexModel(PUShuntHub(; name, kwargs...), [:i_r,:i_i], [:u_r,:u_i]; verbose=false)
function pu_dqswing_injector(; name, kwargs...)
    VertexModel(PUSwingInjector(; name, kwargs...), [:u_r,:u_i], [:i_r,:i_i];
                ff_to_constraint=false, assume_io_coupling=true, verbose=false)
end
function pu_dqline(; name=:line, R=0.04, X=0.1, kwargs...)
    EdgeModel(PUStaticLine(; name, R, X),
              [:src_u_r,:src_u_i], [:dst_u_r,:dst_u_i], [:src_i_r,:src_i_i], [:dst_i_r,:dst_i_i];
              verbose=false, kwargs...)
end

"Two hubs joined by a line, each carrying a swing injector via a LoopbackConnection. The hubs hold
distinct system bases (100 / 200 MVA); the base flows into each injector (`default_from (:hub,…)`),
each injector's rating (`initf_weak`), and the line (`default_from (:src,…)`). Powerflow-then-
dynamic like [`powergridlike_injector_network`](@ref). Returns `(nw, s0)`."
function pu_injector_network()
    edges = [
        pu_dqline(R=0.04, X=0.1, src=:hub1, dst=:hub2),
        LoopbackConnection(potential=[:u_r,:u_i], flow=[:i_r,:i_i], src=:inj1, dst=:hub1, name=:loopback1),
        LoopbackConnection(potential=[:u_r,:u_i], flow=[:i_r,:i_i], src=:inj2, dst=:hub2, name=:loopback2),
    ]
    # powerflow: PU hubs + plain slack/PV injectors (the base is inert here)
    pf_vertices = [
        pu_dqhub(name=:hub1), pu_dqhub(name=:hub2),
        dqbus_slack(injector=true, name=:inj1),
        dqbus_pv(injector=true, Pset=1.0, Vset=1.0, name=:inj2),
    ]
    nws = Network(pf_vertices, edges; warn_order=false)
    pf = find_fixpoint(nws)

    # dynamic model: swing injectors, distinct base per hub
    dyn_vertices = [
        pu_dqhub(name=:hub1), pu_dqhub(name=:hub2),
        pu_dqswing_injector(name=:inj1), pu_dqswing_injector(name=:inj2),
    ]
    set_default!(dyn_vertices[1], :busbar₊S_b, 100.0)
    set_default!(dyn_vertices[2], :busbar₊S_b, 200.0)
    nw = Network(dyn_vertices, edges; warn_order=false)
    s0 = initialize_componentwise!(nw; subverbose=false, verbose=false, default_overrides=interface_values(pf))
    nw, s0
end

end #module
