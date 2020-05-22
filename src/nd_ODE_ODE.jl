module nd_ODE_ODE_mod

using ..NetworkStructures
using ..NDFunctions
using ..Utilities

export nd_ODE_ODE

#= The arguments of the vertex functions must be of the form (dv,v,e_s,e_d,p,t),
where dv is the vertex variable derivative, v the vertex variable and e_s and e_d Arrays of edge variables that
have the vertex as source and destination respectively. p and t are as usual.
The arguments of the edge functions must be of the form (de,e,v_s,v_d,p,t),
where de is the derivative of the edge variable, e the edge variable, v_s and v_d the vertex variables of the vertices
the edge has as source and destination respectively. This works for arbitrary dimensional vertex and edge functions, they only
need to fit, i.e. don't do something like edges! = v_s - v_d when v_s and v_d have not the same dimension. =#


# In order to match the type, we need to pass both, a view that matches the type
# to be constructed, and the original array we want to construct a GD on top of.
@inline function prep_gd(dy::T, y::T, x, gd::GraphData{GDB, T, T}, gs) where {GDB, T}
    # println("Type match")
    gd.gdb.v_array = view(x, 1:gs.dim_v)
    gd.gdb.e_array = view(x, gs.dim_v+1:gs.dim_v+gs.dim_e)
    gd
end

@inline function prep_gd(dy, y, x, gd, gs)
    # println("Type mismatch")
    v_array = view(x, 1:gs.dim_v)
    e_array = view(x, gs.dim_v+1:gs.dim_v+gs.dim_e)
    GraphData(v_array, e_array, gs)
end


@Base.kwdef struct nd_ODE_ODE{T1, GL}
    vertices!::T1
    graph_layer::GL
    parallel::Bool
end

function (d::nd_ODE_ODE)(dx, x, p, t)

    gs = d.graph_layer.graph_structure
    es! = d.graph_layer.edges!
    gd = prep_gd(view(dx, 1:2), view(x, 1:2), x, d.graph_layer.graph_data, gs) # prep_gd will have to prep the whol tuple of layers at once
    aggregator = d.graph_layer.aggregator

    @nd_threads d.parallel for i in 1:gs.num_e
            maybe_idx(es!, i).f!(view(dx,gs.e_idx[i] .+ gs.dim_v), gd.e[i], gd.v_s_e[i], gd.v_d_e[i], p_e_idx(p, i), t)
    end

    @nd_threads d.parallel for i in 1:gs.num_v
            maybe_idx(d.vertices!,i).f!(view(dx,gs.v_idx[i]), gd.v[i], p_v_idx(p, i), t, aggregator(gd.e_s_v[i], gd.e_d_v[i]))
    end

    nothing
end


function (d::nd_ODE_ODE)(x, p, t, ::Type{GetGD})
    prep_gd(view(x, 1:2), view(x, 1:2), x, d.graph_layer.graph_data, d.graph_layer.graph_structure)
end

function (d::nd_ODE_ODE)(::Type{GetGS})
    d.graph_layer.graph_structure
end


end #module
