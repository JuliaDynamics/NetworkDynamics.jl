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

end #module
