####
#### Diffusion system
####
Base.@propagate_inbounds function diffusionedge!(v_s, v_d, _, _)
    return v_s[1] - v_d[1]
end
diffusion_edge() = StaticEdge(f=diffusionedge!, dim=1, pdim=0, coupling=AntiSymmetric())

Base.@propagate_inbounds function diffusion_dedge!(de, e, v_s, v_d, _, _)
    de[1] = 100. * (sin(v_s[1] - v_d[1]) - e[1])
    de[2] = 100. * (sin(v_d[1] - v_s[1]) - e[2])
    nothing
end
diffusion_dedge() = ODEEdge(f=diffusion_dedge!, dim=2, pdim=0, coupling=Fiducial())

Base.@propagate_inbounds function diffusionvertex!(dv, _, acc, _, _)
    dv[1] = acc[1]
    nothing
end
diffusion_vertex() = ODEVertex(f=diffusionvertex!, dim=1, pdim=0)

####
#### inhomogenious kuramoto system
####
Base.@propagate_inbounds function kuramoto_edge!(e, θ_s, θ_d, (K,), t)
    e[1] = K * sin(θ_s[1] - θ_d[1])
end
static_kuramoto_edge() = StaticEdge(f=kuramoto_edge!, dim=1, pdim=1, coupling=AntiSymmetric())

Base.@propagate_inbounds function kuramoto_vertex!(dθ, θ, acc, (ω,), t)
    dθ[1] = ω
    dθ[1] += acc[1]
end
kuramoto_vertex_1d() = ODEVertex(f=kuramoto_vertex!, dim=1, pdim=1)

Base.@propagate_inbounds function kuramoto_inertia!(dv, v, acc, (P,), t)
    dv[1] = v[2]
    dv[2] = P - 1. * v[2]
    dv[2] += acc[1]
end
kuramoto_vertex_2d() = ODEVertex(f=kuramoto_inertia!, dim=2, pdim=1);
