module nd_ODE_Static_mod


using ..NetworkStructures

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
    e_array = similar(x, gs.num_e)
    GraphData(x, e_array, gs)
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


function (d::nd_ODE_Static)(x, p, t, ::Type{GetGD})
    gd = prep_gd(x, d.graph_data, d.graph_structure)

    @inbounds begin

    for i in 1:d.graph_structure.num_e
        d.edges![i].f!(gd.e[i], gd.v_s_e[i], gd.v_d_e[i], maybe_idx(p, i+d.graph_structure.num_v), t)
    end

    end

    gd
end


# For compatibility with PowerDynamics

struct StaticEdgeFunction
    nd_ODE_Static::nd_ODE_Static
end


function (sef::StaticEdgeFunction)(x, p, t)
    gd = sef.nd_ODE_Static(x, p, t, GetGD)

    (gd.e_s_v, gd.e_d_v)
end

end #module
