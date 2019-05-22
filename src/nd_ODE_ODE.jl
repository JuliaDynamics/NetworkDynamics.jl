module nd_ODE_ODE_mod

using LightGraphs
using LinearAlgebra
using Parameters

export nd_ODE_ODE

#= nd_ODE_ODE_scalar constructs a (dx,x,p,t)-function from an Array of functions for the vertices,
 edges as well as a graph.
The arguments of the vertex functions must be of the form (dv,v,e_s,e_d,p,t),
where dv is the vertex variable derivative, v the vertex variable and e_s and e_d Arrays of edge variables that
have the vertex as source and destination respectively. p and t are as usual.
The arguments of the edge functions must be of the form (de,e,v_s,v_d,p,t),
where de is the derivative of the edge variable, e the edge variable, v_s and v_d the vertex variables of the vertices
the edge has as source and destination respectively. This works for arbitrary dimensional vertex and edge functions, they only
need to fit, i.e. don't do something like edges! = v_s - v_d when v_s and v_d have not the same dimension. =#

@with_kw struct nd_ODE_ODE
    vertices!
    edges!
    graph
    graph_stucture
end

function (d::nd_ODE_ODE)(dx, x, p, t)
    gs = d.graph_stucture
    @views begin
    for i in 1:gs.num_e
        gs.e_int[gs.e_idx[i]] .= x[gs.e_x_idx[i]]
        d.edges![i].f!(dx[gs.e_x_idx[i]], gs.e_int[gs.e_idx[i]], x[gs.s_idx[i]], x[gs.d_idx[i]], p, t)
    end
    for i in 1:gs.num_v
        d.vertices![i].f!(dx[gs.v_idx[i]], x[gs.v_idx[i]], gs.e_s[i], gs.e_d[i], p, t)
    end
    end
    nothing
end

function (d::nd_ODE_ODE)(dx, x, p::T, t) where T <: AbstractArray
    gs = d.graph_stucture
    @views begin
    for i in 1:gs.num_e
        gs.e_int[gs.e_idx[i]] .= x[gs.e_x_idx[i]]
        d.edges![i].f!(dx[gs.e_x_idx[i]], gs.e_int[gs.e_idx[i]], x[gs.s_idx[i]], x[gs.d_idx[i]], p[i + gs.num_v], t)
    end
    for i in 1:gs.num_v
        d.vertices![i].f!(dx[gs.v_idx[i]], x[gs.v_idx[i]], gs.e_s[i], gs.e_d[i], p[i], t)
    end
    end
    nothing
end

end #module
