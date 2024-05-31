module Lib
using NetworkDynamics

Base.@propagate_inbounds function diffusionedge!(e, v_s, v_d, (p,), _)
    e .= p * (v_s[1] .- v_d[1])
end
function diffusion_edge()
    StaticEdge(f=diffusionedge!, dim=1, pdim=1, pdef=[1], coupling=AntiSymmetric(),
               name=:diff_edge)
end
function diffusion_edge_closure()
    force_closure = rand()
    StaticEdge(dim=1,
               pdim=1, pdef=[1],
               coupling=AntiSymmetric(), name=:diff_edge_closure) do e, v_s, v_d, (p,), _
        e .= p * (v_s[1] .- v_d[1]) .+ 0 * force_closure
    end
end

Base.@propagate_inbounds function diffusionedge_fid!(e, v_s, v_d, (p,), _)
    e[1] = p * (v_s[1] .- v_d[1])
    e[2] = -e[1]
end
function diffusion_edge_fid()
    StaticEdge(f=diffusionedge_fid!,
               dim=2,
               pdim=1, pdef=[1],
               coupling=Fiducial(), name=:diff_edge_fid)
end

Base.@propagate_inbounds function diffusion_dedge!(de, e, v_s, v_d, (τ,), _)
    de[1] = 1 / τ * (sin(v_s[1] - v_d[1]) - e[1])
    de[2] = 1 / τ * (sin(v_d[1] - v_s[1]) - e[2])
    nothing
end
function diffusion_odeedge()
    ODEEdge(f=diffusion_dedge!,
            dim=2, sym=[:e_dst, :e_src],
            pdim=1, pdef=[100], psym=[:τ],
            coupling=Fiducial(), name=:diff_edge_ode)
end

Base.@propagate_inbounds function diffusionvertex!(dv, _, acc, _, _)
    dv[1] = acc[1]
    nothing
end
diffusion_vertex() = ODEVertex(f=diffusionvertex!, dim=1, pdim=0)

function diffusion_vertex_constraint()
    ODEVertex(f=diffusionvertex!, dim=1, pdim=0, mass_matrix=[0])
end

####
#### inhomogenious kuramoto system
####
Base.@propagate_inbounds function kuramoto_edge!(e, θ_s, θ_d, (K,), t)
    e .= K .* sin(θ_s[1] - θ_d[1])
end
function kuramoto_edge()
    StaticEdge(f=kuramoto_edge!,
               dim=1, sym=[:P],
               pdim=1, psym=[:K],
               coupling=AntiSymmetric())
end

Base.@propagate_inbounds function kuramoto_inertia!(dv, v, acc, p, t)
    M, D, Pm = p
    dv[1] = v[2]
    dv[2] = 1 / M * (Pm - D * v[2] + acc[1])
end
function kuramoto_second()
    ODEVertex(f=kuramoto_inertia!,
              dim=2, sym=[:δ, :ω], def=[0, 0],
              pdim=3, psym=[:M, :D, :Pm], pdef=[1, 0.1, 1])
end

end #module
