using NDPrototype
import NetworkDynamics as OldND
####
#### Diffusion system
####
Base.@propagate_inbounds function diffusionedge!(e, v_s, v_d, _, _)
    e .= v_s[1] .- v_d[1]
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
    e .= K .* sin(θ_s[1] - θ_d[1])
end
static_kuramoto_edge() = StaticEdge(f=kuramoto_edge!, dim=1, pdim=1, coupling=AntiSymmetric())
static_kuramoto_edge_old() = OldND.StaticEdge(f=kuramoto_edge!, dim=1, coupling=:antisymmetric)

Base.@propagate_inbounds function kuramoto_vertex!(dθ, θ, acc, (ω,), t)
    dθ[1] = ω
    dθ[1] += acc[1]
end
kuramoto_vertex_1d() = ODEVertex(f=kuramoto_vertex!, dim=1, pdim=1)

Base.@propagate_inbounds function kuramoto_vertex_old!(dθ, θ, edges, (ω,), t)
    dθ[1] = ω
    for e in edges
        dθ[1] += e[1]
    end
end
kuramoto_vertex_1d_old() = OldND.ODEVertex(f=kuramoto_vertex_old!, dim=1)

Base.@propagate_inbounds function kuramoto_inertia!(dv, v, acc, (P,), t)
    dv[1] = v[2]
    dv[2] = P - 1. * v[2]
    dv[2] += acc[1]
end
kuramoto_vertex_2d() = ODEVertex(f=kuramoto_inertia!, dim=2, pdim=1);
Base.@propagate_inbounds function kuramoto_inertia_old!(dv, v, edges, (P,), t)
    dv[1] = v[2]
    dv[2] = P - 1. * v[2]
    for e in edges
        dv[2] += e[1]
    end
end
kuramoto_vertex_2d_old() = OldND.ODEVertex(f=kuramoto_inertia_old!, dim=2);


function homogeneous(N)
    rng = MersenneTwister(1)
    g = watts_strogatz(N, 3, 0.8, seed=1)
    edge = static_kuramoto_edge()
    vertex = kuramoto_vertex_2d()
    p = vcat(randn(rng, nv(g)), randn(rng, ne(g)))
    (p, vertex, edge, g)
end

function heterogeneous(N)
    rng = MersenneTwister(1)
    g = watts_strogatz(N, 3, 0.8, seed=1)
    edge = static_kuramoto_edge()
    vertex = [kuramoto_vertex_1d(), kuramoto_vertex_2d()]
    vertices = vertex[shuffle(rng, vcat([1 for _ in 1:N÷2], [2 for _ in 1:N÷2]))]
    p = vcat(ones(nv(g)), ones(ne(g)))
    # p = vcat(randn(rng, nv(g)), randn(rng, ne(g)))
    (p, vertices, edge, g)
end

function heterogeneous_old(N)
    rng = MersenneTwister(1)
    g = watts_strogatz(N, 3, 0.8, seed=1)
    edge = static_kuramoto_edge_old()
    vertex = [kuramoto_vertex_1d_old(), kuramoto_vertex_2d_old()]
    vertices = vertex[shuffle(rng, vcat([1 for _ in 1:N÷2], [2 for _ in 1:N÷2]))]
    # p = (randn(rng, nv(g)), randn(rng, ne(g)))
    p = (ones(nv(g)), ones(ne(g)))
    (p, vertices, edge, g)
end
