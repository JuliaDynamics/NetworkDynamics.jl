module nd_ODE_ODE_mod

using ..NetworkStructures
using ..ComponentFunctions
using ..Utilities

export nd_ODE_ODE

#= The signature  of the vertex functions is expected to be (dv,v,e_s,e_d,p,t),
The signature of the edge functions is expected to be (de,e,v_s,v_d,p,t). =#


# In order to match the type, we need to pass both, a view that matches the type
# to be constructed, and the original array we want to construct a GD on top of.
@inline function prep_gd(dy::T, y::T, x, gd::GraphData{T, T}, gs) where T
    # println("Type match")
    gd.v_array = view(x, 1:gs.dim_v)
    gd.e_array = view(x, gs.dim_v+1:gs.dim_v+gs.dim_e)
    gd
end

@inline function prep_gd(dy, y, x, gd, gs)
    # println("Type mismatch")
    v_array = view(x, 1:gs.dim_v)
    e_array = view(x, gs.dim_v+1:gs.dim_v+gs.dim_e)
    GraphData(v_array, e_array, gs)
end


@Base.kwdef struct nd_ODE_ODE{G, T, T1, T2}
    vertices!::T1
    edges!::T2
    graph::G
    graph_structure::GraphStruct
    graph_data::GraphData{T, T}
    parallel::Bool
end

function (d::nd_ODE_ODE)(dx, x, p, t)
    gd = prep_gd(view(dx, 1:2), view(x, 1:2), x, d.graph_data, d.graph_structure)

    @nd_threads d.parallel for i in 1:d.graph_structure.num_e
        maybe_idx(d.edges!, i).f!(
            view(dx,d.graph_structure.e_idx[i] .+ d.graph_structure.dim_v),
            get_edge(gd, i),
            get_src_vertex(gd, i),
            get_dst_vertex(gd, i),
            p_e_idx(p, i),
            t)
    end

    @nd_threads d.parallel for i in 1:d.graph_structure.num_v
        maybe_idx(d.vertices!,i).f!(
            view(dx,d.graph_structure.v_idx[i]),
            get_vertex(gd, i),
            get_out_edges(gd, i),
            get_in_edges(gd, i),
            p_v_idx(p, i),
            t)
    end

    nothing
end


function (d::nd_ODE_ODE)(x, p, t, ::Type{GetGD})
    prep_gd(view(x, 1:2), view(x, 1:2), x, d.graph_data, d.graph_structure)
end

function (d::nd_ODE_ODE)(::Type{GetGS})
    d.graph_structure
end


end #module
