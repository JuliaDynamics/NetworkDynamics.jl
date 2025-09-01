module Lib
using NetworkDynamics
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as Dt

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

@mtkmodel PQLoad begin
    @extend BusBase()
    @parameters begin
        Pset, [description = "Active Power setpoint", guess=0]
        Qset, [description = "Reactive Power setpoint", guess=0]
    end
    @equations begin
        0 ~ _no_simplify(u_r*i_r + u_i*i_i + Pset)
        0 ~ _no_simplify(u_r*i_i - u_i*i_r + Qset)
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
        0 ~ _no_simplify(u_r*i_r + u_i*i_i + Pfun(t))
        0 ~ _no_simplify(u_r*i_i - u_i*i_r + Qset)
    end
end

_no_simplify(x) = x
@register_symbolic _no_simplify(x)
@mtkmodel PVBus begin
    @extend BusBase()
    @parameters begin
        Pset, [description = "Active Power setpoint"]
        Vset = 1.0, [description = "Voltage setpoint"]
    end
    @equations begin
        0 ~ _no_simplify(u_r*i_r + u_i*i_i + Pset)
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
function dqbus_pq(; kwargs...)
    @named pq = PQLoad(; kwargs...)
    VertexModel(pq, [:i_r, :i_i], [:u_r, :u_i])
end
function dqbus_timedeppq(; kwargs...)
    @named pq_timedep = TimeDependentPQLoad(; kwargs...)
    VertexModel(pq_timedep, [:i_r, :i_i], [:u_r, :u_i])
end
function dqbus_pv(; kwargs...)
    @named pv = PVBus(; kwargs...)
    VertexModel(pv, [:i_r, :i_i], [:u_r, :u_i])
end
function dqbus_slack(; kwargs...)
    @named pv = SlackBus(; kwargs...)
    VertexModel(pv, [:i_r, :i_i], [:u_r, :u_i])
end
function dqbus_swing_and_load(; kwargs...)
    @named swing_and_load = SwingAndLoadDQ(; kwargs...)
    VertexModel(swing_and_load, [:i_r, :i_i], [:u_r, :u_i])
end
function dqline(; kwargs...)
    @named line = StaticPowerLineDQ(; kwargs...)
    EdgeModel(line, [:src_u_r, :src_u_i], [:dst_u_r, :dst_u_i], [:src_i_r, :src_i_i], [:dst_i_r, :dst_i_i])
end

end #module
