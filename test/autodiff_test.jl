using ForwardDiff: jacobian!
using ReverseDiff
using NetworkDynamics
using LightGraphs


N = 100 # number of nodes
k = 10  # average degree
g = barabasi_albert(N, k) # a little more exciting than a bare random graph


### Functions for edges and vertices

@inline Base.@propagate_inbounds function diffusionedge!(e, v_s, v_d, p, t)
    e[1] = p * (v_s[1] - v_d[1])
    nothing
end

@inline Base.@propagate_inbounds function diffusionvertex!(dv, v, p, t, i)
    dv[1] = i[1]
    nothing
end

### Constructing the network dynamics

nd_diffusion_vertex = ODEVertex(f! = diffusionvertex!, dim = 1)
nd_diffusion_edge = StaticEdge(f! = diffusionedge!, dim = 1)


const para  = (nothing, randn(ne(g)))
[1 randn(ne(g)-1)...][1,:]
nd = network_dynamics(nd_diffusion_vertex, nd_diffusion_edge, oriented_edge_sum, g, parallel=false)
x0 = randn(N)

nd.f(x0,x0,p,0)

function adndyn(x0)

    nd.f(x0, x0, p , 0.)
    x0
end

ReverseDiff.jacobian(adndyn, x0)
