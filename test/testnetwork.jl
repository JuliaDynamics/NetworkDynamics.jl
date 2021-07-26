using NDPrototype
using LightGraphs
using Random

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
cb = nw.nl.colorbatches[1]

#############
include("../benchmark/benchmark_utils.jl")

# diffusion test
N = 10
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
# parallel
@benchmark begin
    nd(dx, x0, nothing, 0.0)
end setup = begin
    g = watts_strogatz($N, 3, 0.8, seed=1)
    nd = Network(g, $vertex, $edge, parallel=true)
    x0 = randn(dim(nd))
    dx = similar(x0)
    nd(dx, x0, nothing, 0.0)
end

# odeedge
@benchmark begin
    nd(dx, x0, nothing, 0.0)
end setup = begin
    g = watts_strogatz($N, 3, 0.8, seed=1)
    edge = diffusion_dedge()
    nd = Network(g, $vertex, edge, accdim=1, parallel=true)
    x0 = randn(dim(nd))
    dx = similar(x0)
    nd(dx, x0, nothing, 0.0)
end

N = 10
g = watts_strogatz(N, 3, 0.8, seed=1)
edge = diffusion_dedge()
vertex = diffusion_vertex()
nd1 = Network(g, vertex, edge, accdim=1, parallel=true);
nd2 = Network(g, vertex, edge, accdim=1, parallel=false);
x0 = randn(dim(nd1));
dx1 = zero(x0);
dx2 = zero(x0);
nd1(dx1, x0, nothing, 0.0)
nd2(dx2, x0, nothing, 0.0)
@test dx1 ≈ dx2

# inhomgeneous network
function heterogeneous(N)
    rng = MersenneTwister(1)
    g = watts_strogatz(N, 3, 0.8, seed=1)
    edge = static_kuramoto_edge()
    vertex = [kuramoto_vertex_1d(), kuramoto_vertex_2d()]
    vertices = vertex[shuffle(rng, vcat([1 for _ in 1:N÷2], [2 for _ in 1:N÷2]))]
    p = vcat(randn(rng, nv(g)), randn(rng, ne(g)))
    (p, vertices, edge, g)
end

N = 10
@benchmark begin
    nd(dx, x0, p, 0.0)
end setup = begin
    (p, v, e, g) = heterogeneous($N)
    nd = Network(g, v, e)
    @assert length(p) == pdim(nd)
    x0 = randn(dim(nd))
    dx = similar(x0)
    nd(dx, x0, p, 0.0) # call to init caches, we don't want to benchmark this
end
@benchmark begin
    nd(dx, x0, p, 0.0)
end setup = begin
    (p, v, e, g) = heterogeneous($N)
    nd = Network(g, v, e, parallel=true)
    @assert length(p) == pdim(nd)
    x0 = randn(dim(nd))
    dx = similar(x0)
    nd(dx, x0, p, 0.0) # call to init caches, we don't want to benchmark this
end

begin
    N = 1000
    (p, v, e, g) = heterogeneous(N)
    nd = Network(g, v, e, parallel=true)
    x0 = randn(dim(nd))
    dx = similar(x0)
    nd(dx, x0, p, 0.0) # call to init caches, we don't want to benchmark this
    @btime $nd($dx, $x0, $p, 0.0) # call to init caches, we don't want to benchmark this
end
