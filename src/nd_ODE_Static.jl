module nd_ODE_Static_mod

using ..NetworkStructures
using ..ComponentFunctions
using ..Utilities


export nd_ODE_Static
export StaticEdgeFunction


@inline function prep_gd(dx::T, x::T, gd::GraphData{T, T}, gs) where T
    # Type matching
    gd.v_array = x # Does not allocate
    gd
end

@inline function prep_gd(dx, x, gd, gs)
    # Type mismatch
    e_array = similar(dx, gs.dim_e)
    GraphData(x, e_array, gs)
end



@Base.kwdef struct nd_ODE_Static{G, Tv, Te, T1, T2}
    vertices!::T1
    edges!::T2
    graph::G #redundant?
    graph_structure::GraphStruct
    graph_data::GraphData{Tv, Te}
    parallel::Bool # enables multithreading for the core loop
end


function (d::nd_ODE_Static)(dx, x, p, t)
    gd = prep_gd(dx, x, d.graph_data, d.graph_structure)

    @nd_threads d.parallel for i in 1:d.graph_structure.num_e
        maybe_idx(d.edges!, i).f!(
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


function (d::nd_ODE_Static)(x, p, t, ::Type{GetGD})
    gd = prep_gd(x, x, d.graph_data, d.graph_structure)


    @nd_threads d.parallel for i in 1:d.graph_structure.num_e
        maybe_idx(d.edges!,i).f!(
            get_edge(gd, i),
            get_src_vertex(gd, i),
            get_dst_vertex(gd, i),
            p_e_idx(p, i),
            t)
    end

    gd
end


function (d::nd_ODE_Static)(::Type{GetGS})
    d.graph_structure
end


# For compatibility with PowerDynamics

struct StaticEdgeFunction
    nd_ODE_Static
end

function (sef::StaticEdgeFunction)(x, p, t)
    gd = sef.nd_ODE_Static(x, p, t, GetGD)

    # Allocating but non-breaking. No special knowlage about the fields of
    # GraphData needed anymore. Is there a real need for this in PD.jl?
    num_v = sef.nd_ODE_Static(GetGS).num_v
    ([get_out_edges(gd, i) for i in 1:gs.num_v],
     [get_in_edges(gd, i) for i in 1:gs.num_v])
end

end #module
