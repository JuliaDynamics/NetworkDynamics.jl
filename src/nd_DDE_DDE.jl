module nd_DDE_DDE_mod

using Parameters
using LightGraphs
using LinearAlgebra

export nd_DDE_DDE


#= This is under construction. It already works, but shows weird behaviour, for example
if one tests a diffusive network, the variables keep oscillating when they should converge. =#

# See the other nd_ Constructors for details on the fields.

@with_kw struct nd_DDE_DDE
    vertices!
    edges!
    s_e
    d_e
    num_v
    num_e
    e_int
    e_idx
    e_x_idx
    s_idx
    d_idx
    v_idx
    e_s
    e_d
    dim_v
    dim_e
    tau_s # Array of Arrays of Vertex Delays for the outgoing edges.
    tau_d # Array of Arrays of Vertex Delays for the incoming edges.
end

#= This function is needed for convenience in the function call, i.e. that I can
give a function like vertex! = -h(p, t - tau), one does not have to give the delays
in an awkward fashion. =#

struct indexed_h
    h
    idxs
end

function (ih::indexed_h)(args...)
    ih.h(args...; idxs = ih.idxs)
end

function e_s_delayed(h, p, t, tau_s, s_e, i, dim_e, num_e, e_x_idx)
    [[h(p,t-tau_s[i][j]; idxs = e_x_idx[k][j])[1] for j in 1:dim_e[k]] for k in 1:num_e if s_e[k] == i]
end

function e_d_delayed(h, p, t, tau_d, d_e, i, dim_e, num_e, e_x_idx)
    [[h(p,t-tau_d[i][j]; idxs = e_x_idx[k][j])[1] for j in 1:dim_e[k]] for k in 1:num_e if d_e[k] == i]
end

function (d::nd_DDE_DDE)(dx, x, h, p, t)
    @views begin
    for i in 1:d.num_e
        d.e_int[d.e_idx] .= x[d.e_x_idx]
        d.edges![i](dx[d.e_x_idx[i]], x[d.e_x_idx[i]], indexed_h(h,d.e_x_idx[i]), x[d.s_idx[i]], x[d.d_idx[i]], indexed_h(h,d.s_idx[i]), indexed_h(h,d.d_idx[i]), p[d.num_v + i], t)
    end
    for i in 1:d.num_v
        d.vertices![i](dx[d.v_idx[i]], x[d.v_idx[i]], indexed_h(h,d.v_idx[i]), d.e_s[i], d.e_d[i], e_s_delayed(h, p, t, d.tau_s, d.s_e, i, d.dim_e, d.num_e, d.e_x_idx), e_d_delayed(h, p, t, d.tau_d, d.t_e, i, d.dim_e, d.num_e, d.e_x_idx), p[i], t)
    end
    end
    nothing
end

function (d::nd_DDE_DDE)(dx, x, h, p::Nothing, t)
    @views begin
    for i in 1:d.num_e
        d.e_int[d.e_idx[i]] .= x[d.e_x_idx[i]]
        d.edges![i](view(dx,d.e_x_idx[i]),x[d.e_x_idx[i]], indexed_h(h,d.e_x_idx[i]), x[d.s_idx[i]], x[d.d_idx[i]], indexed_h(h,d.s_idx[i]), indexed_h(h,d.d_idx[i]), p, t)
    end
    for i in 1:d.num_v
        d.vertices![i](view(dx,d.v_idx[i]), x[d.v_idx[i]], indexed_h(h,d.v_idx[i]), d.e_s[i], d.e_d[i], e_s_delayed(h, p, t, d.tau_s, d.s_e, i, d.dim_e, d.num_e, d.e_x_idx), e_d_delayed(h, p, t, d.tau_d, d.d_e, i, d.dim_e, d.num_e, d.e_x_idx), p, t)
    end
    end
    nothing
end

function nd_DDE_DDE(vertices!, edges!, s_e, d_e, dim_v, dim_e, tau_s, tau_d)
    num_v = length(dim_v)
    num_e = length(dim_e)

    e_int = zeros(sum(dim_e))

    # x has length sum(dim_v)
    # v_idx is an Array of Array of indices of variables in x belonging to vertex

    counter = 1
    v_idx = [zeros(Int32, dim) for dim in dim_v]
    for i in 1:num_v
        v_idx[i] .= collect(counter:counter + dim_v[i] - 1)
        counter += dim_v[i]
    end

    counter = 1
    e_idx = [zeros(Int32, dim) for dim in dim_e]
    for i in 1:num_e
        e_idx[i] .= collect(counter:counter + dim_e[i] - 1)
        counter += dim_e[i]
    end

    e_x_idx = [e_idx[i] .+ sum(dim_v) for i in 1:num_e]
    # For every vertex, and for every edge, if the source of the edge is that vertex, take the view on the variables of that edge.
    # Thus e_s[i] is an array of views onto the variables of the edges for which i is the source.
    e_s = [[view(e_int, e_idx[i_e]) for i_e in 1:num_e if i_v == s_e[i_e]] for i_v in 1:num_v]
    e_d = [[view(e_int, e_idx[i_e]) for i_e in 1:num_e if i_v == d_e[i_e]] for i_v in 1:num_v]

    s_idx = [v_idx[s_e[i_e]] for i_e in 1:num_e]
    d_idx = [v_idx[d_e[i_e]] for i_e in 1:num_e]

    nd_DDE_DDE(
    vertices!,
    edges!,
    s_e,
    d_e,
    num_v, # Number of vertices
    num_e, # Number of edges
    e_int, # Variables living on edges
    e_idx, # Array of Array of indices of variables in e_int belonging to edges
    e_x_idx, #Array of Array of indices of variables in x belonging to edges
    s_idx, # Array of Array of indices of variables in x belonging to source vertex of edge
    d_idx, # Array of Array of indices of variables in x belonging to destination vertex of edge
    v_idx, # Array of Array of indices of variables in x belonging to vertex
    e_s, # Array of Array of views on the variables in e_int of the edges that are source of a vertex
    e_d, # Array of Array of views on the variables in e_int of the edges that are destination of a vertex
    dim_v,
    dim_e,
    tau_s,
    tau_d
    )
end

"""
    When called with a graph, we construct the source and target vectors."""
function nd_DDE_DDE(vertices!, edges!, g::AbstractGraph, dim_v, dim_e, tau_s, tau_d)
    s_e = [src(e) for e in edges(g)]
    d_e = [dst(e) for e in edges(g)]
    nd_DDE_DDE(vertices!, edges!, s_e, d_e, dim_v, dim_e, tau_s, tau_d)
end

end # module
