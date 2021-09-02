module nd_ODE_DDE_combined_mod

using ..NetworkStructures
using ..ComponentFunctions
using ..Utilities

export nd_ODE_DDE_combined

@inline function prep_gd(dx::AbstractArray{T}, x::AbstractArray{T}, gd::GraphData{GDB, T, T}, gs) where {GDB, T}
    # Type matching
    if size(x) == (gs.dim_v,)
         swap_v_array!(gd, x)
         return gd
    else
         error("Size of x does not match the dimension of the system.")
    end
end

@inline function prep_gd(dx, x, gd, gs)
    # Type mismatch
    if size(x) == (gs.dim_v,)
        e_array = similar(dx, gs.dim_e)
        return GraphData(x, e_array, gs)
    else
        error("Size of x does not match the dimension of the system.")
   end
end

# component_loop for ODE and DDE vertices
function component_loop!(dx, p, t, gd, gs,
                              unique_components, unique_c_indices, history, parallel)
    for j in 1:length(unique_components)
        # Function barrier
        _inner_loop!(dx, p, t, gd, gs,
                           unique_components[j], unique_c_indices[j], history, parallel)
    end
    return nothing
end

# inner loops for ODE + Static case

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
                         component::T, indices, parallel) where T <: StaticEdge
    @nd_threads parallel for i in indices
        component.f!(get_edge(gd, i),
                     get_src_vertex(gd, i),
                     get_dst_vertex(gd, i),
                     p_e_idx(p, i),
                     t)
    end
    return nothing
end

# inner loops for DDE + Static Delay

function _inner_loop!(dx, p, t, gd, gs,
                         component::T, indices, history, parallel) where T <: DDEVertex
    @nd_threads parallel for i in indices
        component.f!(view(dx, gs.v_idx[i]),
                  get_vertex(gd, i),
                  get_dst_edges(gd, i),
                  view(history, gs.v_idx[i]),
                  p_v_idx(p, i),
                  t)
    end
    return nothing
end

function _inner_loop!(dx, p, t, gd, gs,
                         component::T, indices, history, parallel) where T <: StaticDelayEdge
    @nd_threads parallel for i in indices
        component.f!(get_edge(gd, i),
                     get_src_vertex(gd, i),
                     get_dst_vertex(gd, i),
                     view(history, gs.s_e_idx[i]),
                     view(history, gs.d_e_idx[i]),
                     p_e_idx(p, i),
                     t)
    end
    return nothing
end

# struct for both cases

@Base.kwdef struct nd_ODE_DDE_combined{G, GDB, elV, elE, TUV, TUE, Th<:AbstractArray{elV}}
#@Base.kwdef struct nd_ODE_DDE_combined{G, GDB, elV, elE, TUV, TUE, Th<:Union{AbstractArray{elV}, nothing}}
    unique_vertices!::TUV
    unique_v_indices::Vector{Vector{Int}}
    unique_edges!::TUE
    unique_e_indices::Vector{Vector{Int}}
    graph::G #redundant?
    graph_structure::GraphStruct
    graph_data::GraphData{GDB, elV, elE}
#    history::Union{Th, nothing} # for ODE + Static case: nothing, for DDE + Static Delay case: Th
    history::Th # for ODE + Static case: nothing, for DDE + Static Delay case: Th
    parallel::Bool # enables multithreading for the core loop
end

function (d::nd_ODE_DDE_combined)(dx, x, h!, p, t)
    gs = d.graph_structure
    checkbounds_p(p, gs.num_v, gs.num_e)
    gd = prep_gd(dx, x, d.graph_data, d.graph_structure)
    h!(d.history, p, t - p[end])

    @assert size(dx) == size(x) "Sizes of dx and x do not match"

    # Pass nothing, because here we have only Static Edges or Static Delay Edges (if we
    # include the ODE ODE case than we have to pass dx)
    component_loop!(nothing, p, t, gd, gs,
                 d.unique_edges!, d.unique_e_indices, d.history, d.parallel)

    component_loop!(dx, p, t, gd, gs,
                 d.unique_vertices!, d.unique_v_indices, d.history, d.parallel)
    nothing
end

function (d::nd_ODE_DDE_combined)(x, p, t, ::Type{GetGD})
    gd = prep_gd(x, x, d.graph_data, d.graph_structure)
    gs = d.graph_structure
    checkbounds_p(p, gs.num_v, gs.num_e)

    component_loop!(nothing, p, t, gd, gs,
                    d.unique_edges!, d.unique_e_indices, d.history, d.parallel)

    gd
end

function (d::nd_ODE_DDE_combined)(::Type{GetGS})
    d.graph_structure
end

end # module
