module nd_DDE_Static_mod

using ..NetworkStructures
using ..ComponentFunctions
using ..Utilities


export nd_DDE_Static


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

# component_loop for DDE vertices
function component_loop!(dx, p, t, gd, gs,
                              unique_components, unique_c_indices, history, parallel)
    for j in 1:length(unique_components)
        # Function barrier
        _inner_loop!(dx, p, t, gd, gs,
                           unique_components[j], unique_c_indices[j], history, parallel)
    end
    return nothing
end


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

"""
nd_DDE_Static

vertices! has signature (dv, v, e_s, e_d, h_v, p ,t)
edges! has signature (e, v_s, v_d, h_v_s, h_v_d, p, t)
"""
#@Base.kwdef struct nd_DDE_Static{G, T1, T2, GDB, elV, elE, Th<:AbstractArray{elV}}
@Base.kwdef struct nd_DDE_Static{G, GDB, elV, elE, TUV, TUE, Th<:AbstractArray{elV}}
#@Base.kwdef struct nd_DDE_Static{G, GD, TUV, TUE, Th<:AbstractArray{elV}}
    unique_vertices!::TUV
    unique_v_indices::Vector{Vector{Int}}
    unique_edges!::TUE
    unique_e_indices::Vector{Vector{Int}}
    graph::G #redundant?
    graph_structure::GraphStruct
    graph_data::GraphData{GDB, elV, elE}
    history::Th # fixing the type to Th is my first guess, maybe the history should be stored
    # in the GraphData eventually
    parallel::Bool # enables multithreading for the core loop
end


function (d::nd_DDE_Static)(dx, x, h!, p, t)
    gs = d.graph_structure
    checkbounds_p(p, gs.num_v, gs.num_e)
    gd = prep_gd(dx, x, d.graph_data, d.graph_structure)
    h!(d.history, p, t - p[end])

    @assert size(dx) == size(x) "Sizes of dx and x do not match"

    # Pass nothing, because here we have only Static Edges
    component_loop!(nothing, p, t, gd, gs,
                 d.unique_edges!, d.unique_e_indices, d.history, d.parallel)

    component_loop!(dx, p, t, gd, gs,
                 d.unique_vertices!, d.unique_v_indices, d.history, d.parallel)
    nothing
end

function (d::nd_DDE_Static)(x, p, t, ::Type{GetGD})
    gd = prep_gd(x, x, d.graph_data, d.graph_structure)
    gs = d.graph_structure
    checkbounds_p(p, gs.num_v, gs.num_e)

    component_loop!(nothing, p, t, gd, gs,
                    d.unique_edges!, d.unique_e_indices, d.history, d.parallel)

    gd
end

function (d::nd_DDE_Static)(::Type{GetGS})
    d.graph_structure
end

end # module
