using NDPrototype
using LightGraphs

@inline function diffusion_edge!(e,v_s,v_d,p,t)
    @inbounds for i in 1:length(e)
        e[i] = v_s[i] - v_d[i]
    end
    nothing
end

@inline function diffusion_vertex!(dv, v, agg, p, t)
    @inbounds for i in 1:length(dv)
        dv[i] = agg[i]
    end
    nothing
end

N = 10
g = SimpleGraph(watts_strogatz(N, 4, 0.5, seed=1))

# g = complete_graph(2)
vdim = 1
edim = 1
odevertex = ODEVertex(f = diffusion_vertex!, dim = vdim, pdim = 0)
staticedge = StaticEdge(f = diffusion_edge!, dim = edim, pdim = 0, coupling=AntiSymmetric())

nw = Network(g, odevertex, staticedge, verbose=true);
u = rand(dim(nw))
du = zeros(size(u)...)
p = Float64[]

@btime $nw($du, $u, $p, 0.0)

#############

odevertexv = [odevertex for i in 1:nv(g)];
staticedgev = [staticedge for i in 1:ne(g)];

nw = Network(g, odevertexv, staticedgev);

u = rand(dim(nw))
# u = [0.0, 1.0]
du = zeros(size(u)...)
p = Float64[]

nw(du, u, Float64[], 0.0)

#############

odevertex = [ODEVertex(f = (du,u,a,p,t)->diffusion_vertex!(du,u,a,p,i), dim = vdim, pdim = 0) for i in 1:nv(g)]
staticedge = [StaticEdge(f = (e,vs,vd,p,t)->diffusion_edge!(e,vs,vd,p,i), dim = edim, pdim = 0, coupling=AntiSymmetric()) for i in 1:ne(g)]

nw = Network(g, odevertex, staticedge);

u = rand(dim(nw))
# u = [0.0, 1.0]
du = zeros(size(u)...)
p = Float64[]

@btime $nw($du, $u, $p, 0.0)

#############
include("../benchmark/benchmark_utils.jl")

N = 10000
vertex = diffusion_vertex()
edge = diffusion_edge()
g = watts_strogatz(N, 3, 0.8, seed=1)
nd = Network(g, vertex, edge; verbose=true);
x0 = randn(dim(nd));
dx = similar(x0);
@btime $nd($dx, $x0, nothing, 0.0)

keys(nd.nl.cachepool.caches)

N = 100
vertex = diffusion_vertex()
edge = diffusion_edge()
@benchmark begin
    nd(dx, x0, nothing, 0.0)
end setup = begin
    g = watts_strogatz($N, 3, 0.8, seed=1)
    nd = Network(g, $vertex, $edge)
    x0 = randn(dim(nd))
    dx = similar(x0)
    nd(dx, x0, nothing, 0.0)
end
