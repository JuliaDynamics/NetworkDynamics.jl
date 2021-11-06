# In order to match the type, we need to pass both, a view that matches the type
# to be constructed, and the original array we want to construct a GD on top of.
@inline function prep_gd(dx::AbstractArray{T}, x::AbstractArray{T}, gd::GraphData{GDB,T,T}, gs) where {GDB,T}
    # Type matching
    if size(x) == (gs.dim_v,)
        swap_v_array!(gd, x)
        return gd
    elseif size(x) == (gs.dim_v + gs.dim_e,)
        swap_v_array!(gd, view(x, 1:gs.dim_v))
        swap_e_array!(gd, view(x, gs.dim_v+1:gs.dim_v+gs.dim_e))
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
    elseif size(x) == (gs.dim_v + gs.dim_e,)
        v_array = view(x, 1:gs.dim_v)
        e_array = view(x, gs.dim_v+1:gs.dim_v+gs.dim_e)
        return GraphData(v_array, e_array, gs)
    else
        error("Size of x does not match the dimension of the system.")
    end
end


@inline function component_loop!(unique_components, unique_c_indices,
                                 dx, p, t, gd, gs, history, parallel)
    for j in 1:length(unique_components)
        # Function barrier
        _inner_loop!(unique_components[j], unique_c_indices[j],
                     dx, p, t, gd, gs, history, parallel)
    end
    return nothing
end


# inner loops for ODE + Static case

function _inner_loop!(component::ODEVertex, indices,
                      dx, p, t, gd, gs, history, parallel)
    @nd_threads parallel for i in indices
        component.f(view(dx, gs.v_idx[i]),
                    get_vertex(gd, i),
                    get_dst_edges(gd, i),
                    p_v_idx(p, i),
                    t)
    end
    return nothing
end

function _inner_loop!(component::StaticEdge, indices,
                      dx, p, t, gd, gs, history, parallel)
    @nd_threads parallel for i in indices
        component.f(get_edge(gd, i),
                    get_src_vertex(gd, i),
                    get_dst_vertex(gd, i),
                    p_e_idx(p, i),
                    t)
    end
    return nothing
end

# inner loops for DDE + Static Delay

function _inner_loop!(component::DDEVertex, indices,
                      dx, p, t, gd, gs, history, parallel)
    @nd_threads parallel for i in indices
        component.f(view(dx, gs.v_idx[i]),
                    get_vertex(gd, i),
                    get_dst_edges(gd, i),
                    view(history, gs.v_idx[i]),
                    p_v_idx(p, i),
                    t)
    end
    return nothing
end

function _inner_loop!(component::StaticDelayEdge, indices,
                      dx, p, t, gd, gs, history, parallel)
    @nd_threads parallel for i in indices
        component.f(get_edge(gd, i),
                    get_src_vertex(gd, i),
                    get_dst_vertex(gd, i),
                    view(history, gs.s_e_idx[i]),
                    view(history, gs.d_e_idx[i]),
                    p_e_idx(p, i),
                    t)
    end
    return nothing
end

function _inner_loop!(component::ODEEdge, indices,
                      dx, p, t, gd, gs, history, parallel)
    @nd_threads parallel for i in indices
        component.f(view(dx, gs.e_idx[i] .+ gs.dim_v),
                    get_edge(gd, i),
                    get_src_vertex(gd, i),
                    get_dst_vertex(gd, i),
                    p_e_idx(p, i),
                    t)
    end
    return nothing
end

# struct for both cases

#@Base.kwdef struct NetworkDE{G, GDB, elV, elE, TUV, TUE, Th<:AbstractArray{elV}}
Base.@kwdef struct NetworkDE{G,GDB,elV,elE,TUV,TUE,Th<:Union{AbstractArray,Nothing}}
    unique_vertices!::TUV
    unique_v_indices::Vector{Vector{Int}}
    unique_edges!::TUE
    unique_e_indices::Vector{Vector{Int}}
    graph::G #redundant?
    graph_structure::GraphStruct
    graph_data::GraphData{GDB,elV,elE}
    history::Th # for ODE + Static case: nothing, for DDE + Static Delay case: Th
    parallel::Bool # enables multithreading for the core loop
end

# for ODE case

function (d::NetworkDE)(dx, x, p, t)
    gs = d.graph_structure
    checkbounds_p(p, gs.num_v, gs.num_e)
    gd = prep_gd(dx, x, d.graph_data, gs)

    @assert size(dx) == size(x) "Sizes of dx and x do not match"

    # Pass nothing, because here we have only Static Edges or Static Delay Edges (if we
    # include the ODE ODE case than we have to pass dx)
    component_loop!(d.unique_edges!, d.unique_e_indices,
                    dx, p, t, gd, gs, d.history, d.parallel)

    component_loop!(d.unique_vertices!, d.unique_v_indices,
                    dx, p, t, gd, gs, d.history, d.parallel)
    return nothing
end
# for DDE case
function (d::NetworkDE)(dx, x, h!, p, t)
    # History computation happens beforehand and is cached in d.history
    h!(d.history, p, t - p[end])
    d(dx, x, p, t)
    return nothing
end

function (d::NetworkDE)(x, p, t, ::Type{GetGD})
    gs = d.graph_structure
    gd = prep_gd(x, x, d.graph_data, gs)
    # For networks with ODE edges all edge data x
    # Such network have size(x) == (gs.dim_v + gs.dim_e)
    if size(x) == (gs.dim_v,)
        checkbounds_p(p, gs.num_v, gs.num_e)
        component_loop!(d.unique_edges!, d.unique_e_indices,
                        nothing, p, t, gd, gs, d.history, d.parallel)
    end
    gd
end

function (d::NetworkDE)(::Type{GetGS})
    d.graph_structure
end
