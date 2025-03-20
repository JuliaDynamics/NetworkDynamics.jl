using NetworkDynamics
using Graphs: complete_graph

Base.@propagate_inbounds function diffusionedge!(e, v_s, v_d, (p,), _)
    e .= p * (v_s[1] .- v_d[1])
end

Base.@propagate_inbounds function diffusionvertex!(dv, _, acc, _, _)
    dv[1] = acc[1]
    nothing
end

g = complete_graph(3)
edge = EdgeModel(g=AntiSymmetric(diffusionedge!), outdim=1, pdim=1, pdef=[1])
vert = VertexModel(f=diffusionvertex!, g=1, dim=1, pdim=0)
nd = Network(g, vert, edge)
u = rand(dim(nd))
du = similar(u)
p = rand(pdim(nd))
nd(du,u,p,0.0)
