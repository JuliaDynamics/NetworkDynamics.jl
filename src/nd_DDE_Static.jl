module nd_DDE_Static_mod

using ..NetworkStructures
using ..NDFunctions
using ..Utilities


export nd_DDE_Static


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



@Base.kwdef struct nd_DDE_Static{G, Tv, Te, T1, T2}
    vertices!::T1
    edges!::T2
    graph::G #redundant?
    graph_structure::GraphStruct
    graph_data::GraphData{Tv, Te}
    history::Tv # fixing the type to Tv is my first guess, maybe the histry should be stored
    # in the GraphData eventually
    parallel::Bool # enables multithreading for the core loop
end


function (d::nd_DDE_Static)(dx, x, h!, p, t)
    gd = prep_gd(dx, x, d.graph_data, d.graph_structure)
    h!(d.history, p, t)

    @nd_threads d.parallel for i in 1:d.graph_structure.num_e
        maybe_idx(d.edges!, i).f!(gd.e[i], gd.v_s_e[i], gd.v_d_e[i], p_e_idx(p, i), t)
    end

    @nd_threads d.parallel for i in 1:d.graph_structure.num_v
        maybe_idx(d.vertices!,i).f!(view(dx,d.graph_structure.v_idx[i]), gd.v[i], gd.e_s_v[i], gd.e_d_v[i], view(d.history, d.graph_structure.v_idx[i]), p_v_idx(p, i), t)
    end

    nothing
end

end
