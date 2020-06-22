using ForwardDiff
using ReverseDiff
import Zygote
using NetworkDynamics
using LightGraphs


N = 100 # number of nodes
k = 10  # average degree
g = barabasi_albert(N, k) # a little more exciting than a bare random graph


### Functions for edges and vertices


@inline Base.@propagate_inbounds function diffusionedge!(e, v_s, v_d, p, t)
    # usually e, v_s, v_d are arrays, hence we use the broadcasting operator .
    e[1] = p * (v_s[1] - v_d[1])
    nothing
end

@inline Base.@propagate_inbounds function diffusionvertex!(dv, v, e_s, e_d, p, t)
    # usually v, e_s, e_d are arrays, hence we use the broadcasting operator .
    dv[1] = p
    # edges for which v is the source
    for e in e_s
        dv[1] -= e[1]
    end
    # edges for which v is the destination
    for e in e_d
        dv[1] += e[1]
    end
    nothing
end
### Constructing the network dynamics

nd_diffusion_vertex = ODEVertex(f! = diffusionvertex!, dim = 1)
nd_diffusion_edge = StaticEdge(f! = diffusionedge!, dim = 1)

p_v = collect(1:nv(g))./nv(g) .- 0.5
p_e = .5 .* ones(ne(g))
p  = (p_v, p_e)

nd = network_dynamics(nd_diffusion_vertex, nd_diffusion_edge, g, parallel=false)
x0 = randn(N)


## Gradient

function gradnd(x)
    nd.f(x, zeros(N), p, 0.)
    sum(x)
end

ForwardDiff.gradient(gradnd, ones(N))
ReverseDiff.gradient(gradnd, ones(N)) # mutation not supported
Zygote.gradient(gradnd, ones(N))  # mutation not supported


## Gradient wrt. vertex parameters

function gradpnd(p)
    x0int = ones(N)
    nd.f(x0int, zeros(N), (p, p_e), 0.)
    sum(x0int)
end

ForwardDiff.gradient(gradpnd, ones(N)) # type problem
ReverseDiff.gradient(gradpnd, ones(N))
Zygote.gradient(gradpnd, ones(N))  # mutation not supported


## parameter array instead of tuple

@inline Base.@propagate_inbounds function diffusionedge2!(e, v_s, v_d, p, t)

    e[1] = p[2] * (v_s[1] - v_d[1])
    nothing
end

@inline Base.@propagate_inbounds function diffusionvertex2!(dv, v, e_s, e_d, p, t)

    dv[1] = p[1]
    # edges for which v is the source
    for e in e_s
        dv[1] -= e[1]
    end
    # edges for which v is the destination
    for e in e_d
        dv[1] += e[1]
    end
    nothing
end
nd_diffusion_vertex2 = ODEVertex(f! = diffusionvertex2!, dim = 1)
nd_diffusion_edge2 = StaticEdge(f! = diffusionedge2!, dim = 1)


p = [1., 2.]

nd2 = network_dynamics(nd_diffusion_vertex2, nd_diffusion_edge2, g, parallel=false)
nd2.f(x0, zeros(N), p, 0.)

function gradpnd2(p)
    x0 = ones(N)
    nd2.f(x0, zeros(N), p, 0.)
    sum(x0)
end

ForwardDiff.gradient(gradpnd2, p) # type problem
ReverseDiff.gradient(gradpnd2, p)
Zygote.gradient(gradpnd2, p)
