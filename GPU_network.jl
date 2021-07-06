using CUDA
using LightGraphs
using NetworkDynamics
using BenchmarkTools

@inline function diffusion_edge!(e,v_s,v_d,p,t)
    e .= v_s .- v_d
    nothing
end
staticedge = StaticEdge(f! = diffusion_edge!, dim = 1, coupling=:directed)

@inline function diffusion_vertex_nd!(dv, v, edges, p, t)
    dv .= 0.
    sum_coupling!(dv, edges)
    nothing
end
odevertex = ODEVertex(f! = diffusion_vertex_nd!,dim = 1)

@inline function diffusion_vertex!(dv, v, agg::AbstractVector{<:AbstractFloat}, p, t)
    dv .= agg
    nothing
end

struct Network{FloatMat, IdxMat, Vertex, Edge, Aggr}
    "vertex function (ODE)"
    vfun::Vertex
    "dimensions of edge function"
    vdim::Int
    "edge function (static)"
    efun::Edge
    "dimensions of vertex function"
    edim::Int
    "Aggregation function"
    aggfun::Aggr
    "Matrix(2, #edges), [srcidx, dstidx] per column"
    vidx_per_e::IdxMat
    "Matrix(?, #vertices), [e1idx, e2idx,...] of *incoming* edges per column"
    eidx_per_v::IdxMat
    "cache matrix for the aggregation"
    _aggregation::FloatMat
    "cache matrix for the static edges"
    _edge_u::FloatMat
end

function Network(; g::SimpleDiGraph, vf, vdim, ef, edim, aggf, AT=Array, FT=Float64)
    edge_u = Matrix{FT}(undef, edim, ne(g))
    aggregation = Matrix{FT}(undef, edim, nv(g))

    vidx_per_e = Matrix{Int32}(undef, 2, ne(g))
    for (i, e) in enumerate(edges(g))
        vidx_per_e[:, i] .= [e.src, e.dst]
    end

    # indices of the outgoing edges
    # i.e. find idx for occurences of v in the dst row of vidx_per_edge
    eidx_per_v = zeros(Int, maximum(indegree(g)), nv(g))
    for v in 1:nv(g)
        eidx = findall(isequal(v), vidx_per_e[2, :])
        eidx_per_v[1:length(eidx), v] .= eidx
    end

    Network(vf, vdim, ef, edim, aggf,
            AT{Int32,2}(vidx_per_e),
            AT{Int32,2}(eidx_per_v),
            AT{FT,2}(aggregation),
            AT{FT,2}(edge_u))
end

LightGraphs.ne(nw::Network) = size(nw.vidx_per_e, 2)
LightGraphs.nv(nw::Network) = size(nw.eidx_per_v, 2)

function (nw::Network{<:Matrix})(du, u, p, t)
    # calculate edges
    @inbounds for eidx in 1:ne(nw)
        vu = view(nw._edge_u, :, eidx)
        vs = view(u, :, nw.vidx_per_e[1, eidx])
        vd = view(u, :, nw.vidx_per_e[2, eidx])
        nw.efun(vu, vs, vd, p, t)
    end

    # aggregate edge values
    fill!(nw._aggregation, 0)
    @inbounds for vidx in 1:nv(nw)
        agg = view(nw._aggregation, :, vidx)
        for eidx in view(nw.eidx_per_v, :, vidx)
            eidx == 0 && break # no more edges left
            agg .+= view(nw._edge_u, :, eidx)
        end
    end

    # calculate vertex functions
    fill!(du, 0)
    @inbounds for vidx in 1:nv(nw)
        vdu = view(du, :, vidx)
        vu  = view(du, :, vidx)
        agg = view(nw._aggregation, :, vidx)
        nw.vfun(vdu, vu, agg, p, t)
    end
end

# function (nw::Network{<:CuMatrix})(du, u, p, t)
#     println("GPU")
# end

# g = SimpleDiGraph(smallgraph(:house))
g = SimpleDiGraph(watts_strogatz(1_000, 4, 0.5))

nd = network_dynamics(odevertex, staticedge, g)
nw = Network(; g, vf=diffusion_vertex!, vdim=1, ef=diffusion_edge!, edim=1, aggf=+)

u = rand(1, nv(g))
du = zeros(size(u)...)
u_nd  = u[:] # copy and bring in right shape
du_nd = du[:]

nd(du_nd, u_nd, nothing, 0.0)
nw(du, u, nothing, 0.0)
# du_nd .- du[:]
@assert du_nd ≈ du[:]

@benchmark $nd($du_nd, $u_nd, nothing, 0.0)
@benchmark $nw($du, $u, nothing, 0.0)

# net_gpu = Network(; g, vf=nothing, vdim=1, ef=nothing, edim=1, AT=CuArray)
