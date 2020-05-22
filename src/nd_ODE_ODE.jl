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
@inline function prep_gls(view_dx::T, view_x::T, x, gls, e_array::T) where {T}
    # println("Type match")
    gs = gls[1].graph_structure
    total_dim_e = sum((gl.graph_structure.dim_e for gl in gls)...)
    gd = gls[1].graph_data
    gd.gdb.v_array = view(x, 1:gs.dim_v)
    gd.gdb.e_array = view(x, gs.dim_v+1:gs.dim_v+total_dim_e)
    # These dimension calculations need to be changed to support multiple layers.
    gls
end

@inline function prep_gls(view_dx, view_x, x, gls, e_array)
    gs = gls[1].graph_structure
    # println("Type mismatch")
    v_array = view(x, 1:gs.dim_v)
    e_array = view(x, gs.dim_v+1:gs.dim_v+gs.dim_e)
    (GraphLayer(gl.edges!,
        gl.aggregator,
        gl.graph,
        gl.graph_structure,
        GraphData(v_array, e_array, gl.graph_structure)) for gl in gls)
end


@Base.kwdef struct nd_ODE_ODE{T1, GL}
    vertices!::T1
    graph_layers::GL
    parallel::Bool
end

function (d::nd_ODE_ODE)(dx, x, p, t)

    gls = prep_gls(view(dx, 1:2), view(x, 1:2), x, d.graph_layers, d.graph_layers[1].graph_data.gdb.e_array)

    p_offset = 0

    for gl in gls
        gs = gl.graph_structure
        es! = gl.edges!
        gd = gl.graph_data

        @nd_threads d.parallel for i in 1:gs.num_e
            maybe_idx(es!, i).f!(view(dx,gs.e_idx[i] .+ gs.dim_v .+ gd.global_offset), gd.e[i], gd.v_s_e[i], gd.v_d_e[i], p_e_idx(p, i + p_offset), t)
        end

        p_offset += gs.num_e
    end

    gs = iterate(gls)[1].graph_structure # The tuple comprehension in prep_gls provides a generator, thus we have to use iterate to get the first element.
    gd = iterate(gls)[1].graph_data

    @nd_threads d.parallel for i in 1:gs.num_v
        maybe_idx(d.vertices!,i).f!(view(dx,gs.v_idx[i]), gd.v[i], p_v_idx(p, i), t, (gl.aggregator(gl.graph_data.e_s_v[i], gl.graph_data.e_d_v[i]) for gl in  gls)...)
    end

    nothing
end

function (d::nd_ODE_ODE)(x, p, t, ::Type{GetGD})
    gls = prep_gls(view(x, 1:2), view(x, 1:2), x, d.graph_layers, d.graph_layers[1].graph_data.gdb.e_array)

    println("Warning, only returning first layer data. API Change needed.")
    iterate(gls)[1].graph_data
end

function (d::nd_ODE_ODE)(::Type{GetGS})
    d.graph_layers[1].graph_structure
end


end #module
