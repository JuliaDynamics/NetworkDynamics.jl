module nd_ODE_Static_mod

include("NetworkStructures.jl")
using .NetworkStructures

using ..NDFunctions

using Parameters
using LightGraphs
using LinearAlgebra

export nd_ODE_Static
export StaticEdgeFunction

#=  =#

@inline Base.@propagate_inbounds function maybe_idx(p::T, i) where T <: AbstractArray
    p[i]
end

@inline function maybe_idx(p, i)
    p
end

@inline function prep_gd(x::T, gd::GraphData{T}, gs) where T
    gd.v_array = x
    gd
end

@inline function prep_gd(x, gd, gs)
    e_array = similar(x, d.graph_structure.num_e)
    GraphData(x, e_array, d.graph_structure)
end


@Base.kwdef struct nd_ODE_Static{G, T, T1, T2}
    vertices!::T1
    edges!::T2
    graph::G
    graph_structure::GraphStruct
    graph_data::GraphData{T}
end


function (d::nd_ODE_Static)(dx, x, p, t)
    gd = prep_gd(x, d.graph_data, d.graph_structure)

    @inbounds begin

    for i in 1:d.graph_structure.num_e
        d.edges![i].f!(gd.e[i], gd.v_s_e[i], gd.v_d_e[i], maybe_idx(p, i+d.graph_structure.num_v), t)
    end

    for i in 1:d.graph_structure.num_v
        d.vertices![i].f!(view(dx,d.graph_structure.v_idx[i]), gd.v[i], gd.e_s_v[i], gd.e_d_v[i], maybe_idx(p, i), t)
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
        for i in 1:d.num_e
            d.edges![i].f!(gs.e_int[gs.e_idx[i]], x[gs.s_idx[i]], x[gs.d_idx[i]], maybe_idx(p, i + gs.num_v), t)
        end
    end
    (gs.e_s, gs.e_d)
end

end #module
