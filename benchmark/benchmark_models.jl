using NetworkDynamics
####
#### Diffusion system
####
Base.@propagate_inbounds function diffusionedge!(e, _, v_s, v_d, _, _)
    e[1] = v_s[1] - v_d[1]
    nothing
end
diffusion_edge() = UnifiedEdge(; g=AntiSymmetric(diffusionedge!), outdim=1, pdim=0)

Base.@propagate_inbounds function diffusion_dedge!(de, e, v_s, v_d, _, _)
    de[1] = 100.0 * (sin(v_s[1] - v_d[1]) - e[1])
    de[2] = 100.0 * (sin(v_d[1] - v_s[1]) - e[2])
    nothing
end
diffusion_dedge() = UnifiedEdge(; f=diffusion_dedge!, dim=2, pdim=0, g=Fiducial(dst=1:1, src=2:2))

Base.@propagate_inbounds function diffusionvertex!(dv, _, esum, _, _)
    dv[1] = esum[1]
    nothing
end
diffusion_vertex() = UnifiedVertex(; f=diffusionvertex!, dim=1, pdim=0, g=StateMask(1:1))

####
#### inhomogenious kuramoto system
####
Base.@propagate_inbounds function kuramoto_edge!(e, θ_s, θ_d, (K,), t)
    e[1] = K * sin(θ_s[1] - θ_d[1])
end
static_kuramoto_edge() = UnifiedEdge(; g=AntiSymmetric(kuramoto_edge!), outdim=1, pdim=1)

Base.@propagate_inbounds function kuramoto_vertex!(dθ, θ, esum, (ω,), t)
    dθ[1] = ω + esum[1]
end
kuramoto_vertex_1d() = UnifiedVertex(; f=kuramoto_vertex!, pdim=1, sym=[:θ], osym=[:θ])

Base.@propagate_inbounds function kuramoto_inertia!(dv, v, esum, (P,), t)
    dv[1] = v[2]
    dv[2] = P - 1.0 * v[2]
    dv[2] += esum[1]
end
kuramoto_vertex_2d() = UnifiedVertex(; f=kuramoto_inertia!, dim=2, pdim=1, sym=[:θ, :ω]);
