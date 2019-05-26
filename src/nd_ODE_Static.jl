module nd_ODE_Static_mod

include("NetworkStructures.jl")
using .NetworkStructures

using Parameters
using LightGraphs
using LinearAlgebra

export nd_ODE_Static
export StaticEdgeFunction

#= nd_ODE_Static constructs a (dx,x,p,t)-function from an Array of functions for the vertices,
 edges as well as a graph.
The arguments of the vertex functions must be of the form (dv,v,e_s,e_d,p,t),
where dv is the vertex variable derivative, v the vertex variable and e_s and e_d Arrays of edge variables that
have the vertex as source and destination respectively. p and t are as usual.
The arguments of the edge functions must be of the form (e,v_s,v_d,p,t),
where e is the edge variable and v_s and v_d the vertex variables of the vertices
the edge has as source and destination respectively.
This works for multi-dimensional variables as well. =#


@with_kw struct nd_ODE_Static
    vertices!
    edges!
    graph
    graph_stucture
end

function (d::nd_ODE_Static)(dx, x, p, t)
    gs = d.graph_stucture
    @views begin
    for i in 1:gs.num_e
        d.edges![i].f!(gs.e_int[gs.e_idx[i]], x[gs.s_idx[i]], x[gs.d_idx[i]], p, t)
    end
    for i in 1:gs.num_v
        d.vertices![i].f!(dx[gs.v_idx[i]], x[gs.v_idx[i]], gs.e_s[i], gs.e_d[i], p, t)
    end
    end
    nothing
end

function (d::nd_ODE_Static)(dx, x, p::T, t) where T <: AbstractArray
    gs = d.graph_stucture
    @views begin
    for i in 1:gs.num_e
        d.edges![i].f!(gs.e_int[gs.e_idx[i]], x[gs.s_idx[i]], x[gs.d_idx[i]], p[i + gs.num_v], t)
    end
    for i in 1:gs.num_v
        d.vertices![i].f!(dx[gs.v_idx[i]], x[gs.v_idx[i]], gs.e_s[i], gs.e_d[i], p[i], t)
    end
    end
    nothing
end



#= The struct comes in a version where edges and vertices are not arrays but
homogenous across the network.
=#

struct StaticEdgeFunction
    nd_ODE_Static::nd_ODE_Static
end

function StaticEdgeFunction(vertices!, edges!, graph::G) where G <: AbstractGraph
    dim_v = [v.dim for v in vertices!]
    dim_e = [e.dim for e in edges!]
    dim_nd = sum(dim_v)

    e_int = zeros(sum(dim_e))

    graph_stucture = create_graph_structure(graph, dim_v, dim_e, e_int)

    StaticEdgeFunction(nd_ODE_Static(vertices!, edges!, graph, graph_stucture))
end

function (sef::StaticEdgeFunction)(x, p, t)
    d = sef.nd_ODE_Static
    gs = d.graph_stucture
    @views begin
        for i in 1:gs.num_e
            d.edges![i].f!(gs.e_int[gs.e_idx[i]], x[gs.s_idx[i]], x[gs.d_idx[i]], p, t)
        end
    end
    (gs.e_s, gs.e_d)
end

function (sef::StaticEdgeFunction)(x, p::T, t) where T <: AbstractArray
    d = sef.nd_ODE_Static
    gs = d.graph_stucture
    @views begin
        for i in 1:d.num_e
            d.edges![i].f!(gs.e_int[gs.e_idx[i]], x[gs.s_idx[i]], x[gs.d_idx[i]], p[i + gs.num_v], t)
        end
    end
    (gs.e_s, gs.e_d)
end

end #module
