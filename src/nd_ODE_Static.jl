module nd_ODE_Static_mod

using ..NetworkStructures
using ..NDFunctions
using ..Utilities


export nd_ODE_Static
export StaticEdgeFunction

#=  =#

# If the type of dx and x match, we swap out v_array for x
@inline function prep_gd(dx::T, x::T, gd::GraphData{GDB, T, T}, gs) where {GDB, T}
    gd.gdb.v_array = x # Does not allocate
    gd
end

# If the type of dx and x do not match, we swap initialize a new GraphData object
# that is based on the type of dx for the edge buffer.
# Some solvers take the derivative with respect to time, thus x will not be dual
# but dx will be, leading to errors otherwise
@inline function prep_gd(dx, x, gd, gs)
    e_array = similar(dx, gs.dim_e)
    GraphData(x, e_array, gs)
end



@Base.kwdef struct nd_ODE_Static{T1, GL}
    vertices!::T1
    graph_layer::GL
    parallel::Bool # enables multithreading for the core loop
end


function (d::nd_ODE_Static)(dx, x, p, t)

    gs = d.graph_layer.graph_structure
    es! = d.graph_layer.edges!
    gd = prep_gd(dx, x, d.graph_layer.graph_data, gs) # prep_gd will have to prep the whol tuple of layers at once
    aggregator = d.graph_layer.aggregator

    @nd_threads d.parallel for i in 1:gs.num_e
        maybe_idx(es!, i).f!(gd.e[i], gd.v_s_e[i], gd.v_d_e[i], p_e_idx(p, i), t)
    end

    @nd_threads d.parallel for i in 1:gs.num_v
        maybe_idx(d.vertices!,i).f!(view(dx,gs.v_idx[i]), gd.v[i], p_v_idx(p, i), t, aggregator(gd.e_s_v[i], gd.e_d_v[i]))
    end

    nothing
end


function (d::nd_ODE_Static)(x, p, t, ::Type{GetGD})

    gs = d.graph_layer.graph_structure
    es! = d.graph_layer.edges!
    gd = prep_gd(dx, x, d.graph_layer.graph_data, gs) # prep_gd will have to prep the whol tuple of layers at once
    aggregator = d.graph_layer.aggregator

    @nd_threads d.parallel for i in 1:d.graph_structure.num_v
        maybe_idx(d.vertices!,i).f!(view(dx,gs.v_idx[i]), gd.v[i], p_v_idx(p, i), t, aggregator(gd.e_s_v[i], gd.e_d_v[i]))
    end

    gd
end


function (d::nd_ODE_Static)(::Type{GetGS})
    d.graph_layer.graph_structure
end


# For compatibility with PowerDynamics

struct StaticEdgeFunction
    nd_ODE_Static
end

function (sef::StaticEdgeFunction)(x, p, t)
    gd = sef.nd_ODE_Static(x, p, t, GetGD)

    (gd.e_s_v, gd.e_d_v)
end

end #module
