module NetworkDynamics

using Parameters
using LightGraphs
using LinearAlgebra

export diffusive_network_dynamics

"""
    diffusive_network_dynamics(L, nodes)
    L: Matrix
    nodes: scalar function ``x \\rightarrow nodes(x)``"""
@with_kw struct diffusive_network_dynamics{T}
    L::AbstractArray{T,2}
    nodes::Function
end

"""
    (dnd::diffusive_network_dynamics)(dx, x, p, t)

    Calling a struct of type diffusive network dynamics implements the ODE:

        ``\\frac{dx_i}{dt} = nodes(x_i) + \\sum_j L_{ij} x_j``

    dnd = diffusive_network_dynamics(L, nodes)
    dnd(dx, x, p, t)"""

function (dnd::diffusive_network_dynamics)(dx, x, p, t)
    mul!(dx, dnd.L, x) # dx .= L * x
    dx_temp = 0
    for i in 1:length(dx)
        dnd.nodes(dx_temp, x[i], p, t)
        dx[i] = dx_temp - x[i]
    end
    nothing
end

"""
    When called with a graph, the dynamics defaults to using the laplacian.
"""
function diffusive_network_dynamics(g::AbstractGraph, nodes)
    diffusive_network_dynamics(laplacian_matrix(g), nodes)
end

# The main, fully flexible network dynamics implementation.

"""
    The key functions or function arrays are:

    nodes: ``nodes!(dx, x, [l]_s,[l]_t, p, t)``

    lines: ``lines!(dl, l, x_s, x_t, p, t)``

    Given edges ``e``, ans nodes ``n``, as well as an orientation encoded by
    the source function ``s(e)`` and the target function ``t(e)``
    this implements the system of ODEs:

    ``\\frac{dx_n}{dt} = dx_n``

    ``\\frac{dl_e}{dt} = dl_e``

    with ``dx`` and ``dl`` calculated by

    ``[l]_s = [l_e \\text{ if } s(e) = n]``

    ``[l]_t = [l_e \\text{ if } t(e) = n]``

    ``nodes![n](dx_n, x_n, [l]_s, [l]_t, p_n, t)``

    ``lines![e](dl_e, l_e, x_{s(e)}, x_{t(e)}, p_e, t)``

    Alternative design:

    Something that relaxes to a diffusive network would for example be
    implemented by

        lines = (dl, l, x_1, x_2) -> dl .= 1000. * ((x_1 - x_2) - l)
        agg = (list_1, list_2) -> sum(list_1) - sum(list_2)"""


export scalar_static_lines

"""
Documentation!!
"""
@with_kw struct scalar_static_lines
    vertices!
    edges!
    s_e
    t_e
    e
    e_int
    e_s
    e_t
    num_e
    num_v
end


function (d::scalar_static_lines)(dx, x, p, t)
    for i in 1:d.num_e
        d.edges![i](d.e[i], x[d.s_e[i]], x[d.t_e[i]], p, t)
    end
    for i in 1:d.num_v
        d.vertices![i](view(dx,i), x[i], sum.(d.e_s[i]), sum.(d.e_t[i]), p, t)
    end
    nothing
end

function scalar_static_lines(vertices!, edges!, s_e, t_e)
    num_e = length(edges!)
    num_v = length(vertices!)

    # Will get longer once we have more variables per edge
    e_int = rand(num_e)

    # This will be views to more than one variable eventually
    e = [
        view(e_int, i)
        for i in 1:num_e ]

    # Create an array of views into the edges that selects only the sources
    # for each vertex
    e_s = [
         [e[j] for j in 1:num_e if s_e[j] == i]
        for i in 1:num_v ]
    # Create an array of views into the lines that selects only the targets
    # for each node
    e_t = [
        [e[j] for j in 1:num_e if t_e[j] == i]
        for i in 1:num_v ]

    scalar_static_lines(vertices!, edges!, s_e, t_e, e, e_int, e_s, e_t, num_e, num_v)
end

"""
When called with a graph, we construct the source and target vectors.
"""
function scalar_static_lines(vertices!, edges!, g::AbstractGraph)
    s_e = [src(e) for e in edges(g)]
    t_e = [dst(e) for e in edges(g)]
    scalar_static_lines(vertices!, edges!, s_e, t_e)
end


export scalar_dynamic_lines

@with_kw struct scalar_dynamic_lines
    vertices!
    edges!
    s_e
    t_e
    e
    e_int
    e_s
    e_t
    num_e
    num_v
end


function (d::scalar_dynamic_lines)(dx, x, p, t)
    for i in 1:d.num_e
        d.e[i] .= x[d.num_v+i]
        d.edges![i](view(dx,d.num_v+i),x[d.num_v+i], x[d.s_e[i]], x[d.t_e[i]], p, t)
    end
    for i in 1:d.num_v
        d.vertices![i](view(dx,i), x[i], sum.(d.e_s[i]), sum.(d.e_t[i]), p, t)
    end
    nothing
end

function scalar_dynamic_lines(vertices!, edges!, s_e, t_e)
    num_e = length(edges!)
    num_v = length(vertices!)

    # Will get longer once we have more variables per edge
    e_int = rand(num_e)

    # This will be views to more than one variable eventually
    e = [
        view(e_int, i)
        for i in 1:num_e ]

    # Create an array of views into the edges that selects only the sources
    # for each vertex
    e_s = [
         [e[j] for j in 1:num_e if s_e[j] == i]
        for i in 1:num_v ]
    # Create an array of views into the edges that selects only the targets
    # for each vertex
    e_t = [
        [e[j] for j in 1:num_e if t_e[j] == i]
        for i in 1:num_v ]

    scalar_dynamic_lines(vertices!, edges!, s_e, t_e, e, e_int, e_s, e_t, num_e, num_v)
end

"""
    When called with a graph, we construct the source and target vectors."""

function scalar_dynamic_lines(vertices!, edges!, g::AbstractGraph)
    s_e = [src(e) for e in edges(g)]
    t_e = [dst(e) for e in edges(g)]
    scalar_dynamic_lines(vertices!, edges!, s_e, t_e)
end


export static_lines

@with_kw struct static_lines
    edges!
    vertices!
    num_v # Number of vertices
    num_e # Number of edges
    e_int # Variables living on edges
    e_idx # Array of Array of indices of variables in e_int belonging to edges
    s_idx # Array of Array of indices of variables in x belonging to source vertex of edge
    d_idx # Array of Array of indices of variables in x belonging to destination vertex of edge
    v_idx # Array of Array of indices of variables in x belonging to vertex
    e_s # Array of Array of views on the variables in e_int of the edges that are source of a vertex
    e_d # Array of Array of views on the variables in e_int of the edges that are destination of a vertex
end

function (d::static_lines)(dx, x, p, t)
    @views begin

    for i in 1:d.num_e
        # d.l_e[i] = d.l_e_int[l_idx[i]] as view
        d.edges![i](d.e_int[d.e_idx[i]], x[d.s_idx[i]], x[d.d_idx[i]], p, t)
    end
    for i in 1:d.num_v
        d.vertices![i](dx[d.v_idx[i]], x[d.v_idx[i]], d.e_s[i], d.e_d[i], p, t)
    end

    end
    nothing
end

"""
    dim_v is an array of the number of variables per vertex
    dim_e is an array of the number of variables per edge"""

function static_lines(vertices!, edges!, s_e, d_e, dim_v, dim_e)
    num_v = length(dim_v)
    num_e = length(dim_e)

    e_int = zeros(sum(dim_e))

    # x has length sum(dim_v)
    # v_idx is an Array of Array of indices of variables in x belonging to vertex

    counter = 1
    v_idx = [zeros(Int32, dim) for dim in dim_v]
    for i in 1:num_v
        v_idx[i] .= collect(counter:counter + dim_v[i] - 1)
        counter += dim_v[i]
    end

    counter = 1
    e_idx = [zeros(Int32, dim) for dim in dim_e]
    for i in 1:num_e
        e_idx[i] .= collect(counter:counter + dim_e[i] - 1)
        counter += dim_e[i]
    end

    # For every vertex, and for every edge, if the source of the edge is that vertex, take the view on the variables of that edge.
    # Thus e_s[i] is an array of views onto the variables of the edges for which i is the source.
    e_s = [[view(e_int, e_idx[i_e]) for i_e in 1:num_e if i_v == s_e[i_e]] for i_v in 1:num_v]
    e_d = [[view(e_int, e_idx[i_e]) for i_e in 1:num_e if i_v == d_e[i_e]] for i_v in 1:num_v]

    s_idx = [v_idx[s_e[i_e]] for i_e in 1:num_e]
    d_idx = [v_idx[d_e[i_e]] for i_e in 1:num_e]


    static_lines(
    edges!,
    vertices!,
    num_v, # Number of vertices
    num_e, # Number of edges
    e_int, # Variables living on edges
    e_idx, # Array of Array of indices of variables in e_int belonging to edges
    s_idx, # Array of Array of indices of variables in x belonging to source vertex of edge
    d_idx, # Array of Array of indices of variables in x belonging to destination vertex of edge
    v_idx, # Array of Array of indices of variables in x belonging to vertex
    e_s, # Array of Array of views on the variables in e_int of the edges that are source of a vertex
    e_d) # Array of Array of views on the variables in e_int of the edges that are destination of a vertex
end

function static_lines(vertices!, edges!, g::AbstractGraph, dim_v, dim_e)
    s_e = [src(e) for e in edges(g)]
    d_e = [dst(e) for e in edges(g)]
    static_lines(vertices!, edges!, s_e, d_e, dim_v, dim_e)
end


export dynamic_lines

@with_kw struct dynamic_lines
    edges!
    vertices!
    num_v # Number of vertices
    num_e # Number of edges
    e_int # Variables living on edges
    e_idx # Array of Array of indices of variables in e_int belonging to edges
    e_x_idx # Array of Array of indices of variables in x belonging to edges
    s_idx # Array of Array of indices of variables in x belonging to source vertex of edge
    d_idx # Array of Array of indices of variables in x belonging to destination vertex of edge
    v_idx # Array of Array of indices of variables in x belonging to vertex
    e_s # Array of Array of views on the variables in e_int of the edges that are source of a vertex
    e_d # Array of Array of views on the variables in e_int of the edges that are destination of a vertex
end

function (d::dynamic_lines)(dx, x, p, t)
    @views begin

    for i in 1:d.num_e
        d.e_int[d.e_idx[i]] .= x[d.e_x_idx[i]]
        d.edges![i](dx[d.e_x_idx[i]],d.e_int[d.e_idx[i]], x[d.s_idx[i]], x[d.d_idx[i]], p, t)
    end
    for i in 1:d.num_v
        d.vertices![i](dx[d.v_idx[i]], x[d.v_idx[i]], d.e_s[i], d.e_d[i], p, t)
    end

    end
    nothing
end

"""
dim_v is an array of the number of variables per vertex
dim_e is an array of the number of variables per edge
"""
function dynamic_lines(vertices!, edges!, s_e, d_e, dim_v, dim_e)
    num_v = length(dim_v)
    num_e = length(dim_e)

    e_int = zeros(sum(dim_e))

    # x has length sum(dim_v)
    # v_idx is an Array of Array of indices of variables in x belonging to vertex

    counter = 1
    v_idx = [zeros(Int32, dim) for dim in dim_v]
    for i in 1:num_v
        v_idx[i] .= collect(counter:counter + dim_v[i] - 1)
        counter += dim_v[i]
    end

    counter = 1
    e_idx = [zeros(Int32, dim) for dim in dim_e]
    for i in 1:num_e
        e_idx[i] .= collect(counter:counter + dim_e[i] - 1)
        counter += dim_e[i]
    end

    e_x_idx = [e_idx[i] .+ sum(dim_v) for i in 1:num_e]
    # For every vertex, and for every edge, if the source of the edge is that vertex, take the view on the variables of that edge.
    # Thus e_s[i] is an array of views onto the variables of the edges for which i is the source.
    e_s = [[view(e_int, e_idx[i_e]) for i_e in 1:num_e if i_v == s_e[i_e]] for i_v in 1:num_v]
    e_d = [[view(e_int, e_idx[i_e]) for i_e in 1:num_e if i_v == d_e[i_e]] for i_v in 1:num_v]

    s_idx = [v_idx[s_e[i_e]] for i_e in 1:num_e]
    d_idx = [v_idx[d_e[i_e]] for i_e in 1:num_e]

    dynamic_lines(
    edges!,
    vertices!,
    num_v, # Number of vertices
    num_e, # Number of edges
    e_int, # Variables living on edges
    e_idx, # Array of Array of indices of variables in e_int belonging to edges
    e_x_idx, #Array of Array of indices of variables in x belonging to edges
    s_idx, # Array of Array of indices of variables in x belonging to source vertex of edge
    d_idx, # Array of Array of indices of variables in x belonging to destination vertex of edge
    v_idx, # Array of Array of indices of variables in x belonging to vertex
    e_s, # Array of Array of views on the variables in e_int of the edges that are source of a vertex
    e_d) # Array of Array of views on the variables in e_int of the edges that are destination of a vertex
end

function dynamic_lines(vertices!, edges!, g::AbstractGraph, dim_v, dim_e)
    s_e = [src(e) for e in edges(g)]
    d_e = [dst(e) for e in edges(g)]
    dynamic_lines(vertices!, edges!, s_e, d_e, dim_v, dim_e)
end


end # module
