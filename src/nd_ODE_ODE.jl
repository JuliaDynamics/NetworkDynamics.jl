module nd_ODE_ODE_mod

using ..NetworkStructures
using ..ComponentFunctions
using ..Utilities

export nd_ODE_ODE

#= The signature  of the vertex functions is expected to be (dv,v,e_s,e_d,p,t),
The signature of the edge functions is expected to be (de,e,v_s,v_d,p,t). =#


# In order to match the type, we need to pass both, a view that matches the type
# to be constructed, and the original array we want to construct a GD on top of.
@inline function prep_gd(dx::AbstractArray{T}, x::AbstractArray{T}, gd::GraphData{GDB, T, T}, gs) where {GDB, T}
    # println("Type match")
    # We don't need to check the size of x here since view() does that by default.
    swap_v_array!(gd, view(x, 1:gs.dim_v))
    swap_e_array!(gd, view(x, gs.dim_v+1:gs.dim_v+gs.dim_e))
    gd
end

@inline function prep_gd(dx, x, gd, gs)
    # println("Type mismatch")
    v_array = view(x, 1:gs.dim_v)
    e_array = view(x, gs.dim_v+1:gs.dim_v+gs.dim_e)
    GraphData(v_array, e_array, gs)
end

function component_loop!(dx, p, t, gd, gs,
                              unique_components, unique_c_indices, parallel)
    for j in 1:length(unique_components)
        # Function barrier
        _inner_loop!(dx, p, t, gd, gs,
                           unique_components[j], unique_c_indices[j], parallel)
    end
    return nothing
end

function _inner_loop!(dx, p, t, gd, gs,
                         component::T, indices, parallel) where T <: ODEVertex
    @nd_threads parallel for i in indices
        component.f!(view(dx,gs.v_idx[i]),
                  get_vertex(gd, i),
                  get_dst_edges(gd, i),
                  p_v_idx(p, i),
                  t)
    end
    return nothing
end

function _inner_loop!(dx, p, t, gd, gs,
                         component::T, indices, parallel) where T <: ODEEdge
    @nd_threads parallel for i in indices
        component.f!(view(dx, gs.e_idx[i] .+ gs.dim_v),
                     get_edge(gd, i),
                     get_src_vertex(gd, i),
                     get_dst_vertex(gd, i),
                     p_e_idx(p, i),
                     t)
    end
    return nothing
end


@Base.kwdef struct nd_ODE_ODE{G, GD, TUV, TUE}
    unique_vertices!::TUV
    unique_v_indices::Vector{Vector{Int}}
    unique_edges!::TUE
    unique_e_indices::Vector{Vector{Int}}
    graph::G
    graph_structure::GraphStruct
    graph_data::GD
    parallel::Bool
end

function (d::nd_ODE_ODE)(dx, x, p, t)
    gs = d.graph_structure
    checkbounds_p(p, gs.num_v, gs.num_e)
    gd = prep_gd(dx, x, d.graph_data, d.graph_structure)

    @assert size(dx) == size(x) "Sizes of dx and x do not match"

    # Pass nothing, because here we have only Static Edges
    component_loop!(dx, p, t, gd, gs,
                 d.unique_edges!, d.unique_e_indices, d.parallel)

    component_loop!(dx, p, t, gd, gs,
                 d.unique_vertices!, d.unique_v_indices, d.parallel)

    nothing
end


function (d::nd_ODE_ODE)(x, p, t, ::Type{GetGD})
    prep_gd(nothing, x, d.graph_data, d.graph_structure)
end

function (d::nd_ODE_ODE)(::Type{GetGS})
    d.graph_structure
end


end #module
