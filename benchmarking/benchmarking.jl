# Diffusion on a star graph

using NetworkDynamics
using LightGraphs
using OrdinaryDiffEq
using BenchmarkTools


### Functions for edges and vertices

@inline Base.@propagate_inbounds function diffusionedge!(e, v_s, v_d, p, t)
    e[1] = v_s[1] - v_d[1]
    nothing
end

@inline Base.@propagate_inbounds function diffusionvertex!(dv, v, e_s, e_d, p, t)
    dv[1] = 0.
    for e in e_s
        dv[1] -= e[1]
    end
    for e in e_d
        dv[1] += e[1]
    end
    nothing
end

nd_diffusion_vertex = ODEVertex(f! = diffusionvertex!, dim = 1)
nd_diffusion_edge = StaticEdge(f! = diffusionedge!, dim = 1)


for N =[10, 100, 1_000, 10_000, 100_000]
    local g = star_graph(N)
    println("Time to assemble star diffusion with $N nodes: ")
    @time local nd = network_dynamics(nd_diffusion_vertex, nd_diffusion_edge, g)
    local x0 = randn(N)
    println("Time to call star diffusion with $N nodes: ")
    @btime $nd($x0, $x0, nothing, 0.)
end
