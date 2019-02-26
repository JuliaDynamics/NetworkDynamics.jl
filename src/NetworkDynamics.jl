module NetworkDynamics

using Parameters
using LightGraphs
using LinearAlgebra

export diffusive_network_dynamics

"""
    diffusive_network_dynamics(L, nodes)
L: Matrix
nodes: scalar function ``x \\rightarrow nodes(x)``
"""
@with_kw struct diffusive_network_dynamics{T}
    L::AbstractArray{T,2}
    nodes::Function
end

"""
    (dnd::diffusive_network_dynamics)(dx, x, p, t)

Calling a struct of type diffusive network dynamics implements the ODE:

``\\frac{dx_i}{dt} = nodes(x_i) + \\sum_j L_{ij} x_j``

    dnd = diffusive_network_dynamics(L, nodes)
    dnd(dx, x, p, t)
"""

function (dnd::diffusive_network_dynamics)(dx, x, p, t)
    mul!(dx, dnd.L, x) # dx .= L * x
    dx_temp = 0
    for i in 1:length(dx)
        dnd.nodes(dx_temp, x[i], p, t)
        dx[i] = dx_temp - dx[i]
    end
    nothing
end

"""
When called with a graph, the dynamics defaults to using the laplacian.
"""
function diffusive_network_dynamics(g::AbstractGraph, nodes)
    diffusive_network_dynamics(laplacian_matrix(g), nodes)
end

# The main, fully flexible network dynamics implementation
export network_dynamics

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
    agg = (list_1, list_2) -> sum(list_1) - sum(list_2)
"""

@with_kw struct network_dynamics{T}
    B::AbstractArray{T,2}
    nodes::Function
    lines::Function
end

function (nd::network_dynamics)(dv, v, p, t)
    dx_temp=0
    dl_temp=0
    dx=zeros(size(nd.B,1))
    dl=zeros(size(nd.B,2))
    for i in 1:length(dx)
        dv[i]= nd.nodes(dx_temp,v[i],p,t) - sum(Float64[v[k+size(nd.B,1)] for k in 1:size(nd.B,2) if nd.B[i,k]==1])
        + sum(Float64[v[k+size(nd.B,1)] for k in 1:size(nd.B,2) if nd.B[i,k]==-1])
    end
    for j in (length(dx)+1):(length(dx)+length(dl))
        dv[j]=nd.lines(dl_temp,v[j],p,t)
    end
    nothing
end

function network_dynamics(g::AbstractGraph, nodes, lines)
    network_dynamics(incidence_matrix(g, oriented=true), nodes, lines)
end



export network_dynamics_on_the_line

"""
Documentation!!
"""
@with_kw struct network_dynamics_on_the_line{T}
    B::AbstractArray{T,2}
    nodes::Function
    lines::Function
end


function (onl::network_dynamics_on_the_line)(dx, x, p, t)
    dx_temp=0
    l=1
    for i in 1:length(dx)
        onl.nodes(dx_temp,x[i],p,t)
        dx[i]= dx_temp - sum(Float64[onl.lines(l,x,p,t)[k] for k in 1:size(onl.B,2) if onl.B[i,k]==1])
        + sum(Float64[onl.lines(l,x,p,t)[k] for k in 1:size(onl.B,2) if onl.B[i,k]==-1])
    end
    nothing
end

"""
When called with a graph, the dynamics defaults to using the laplacian.
"""
function network_dynamics_on_the_line(g::AbstractGraph, nodes, lines)
    network_dynamics_on_the_line(incidence_matrix(g, oriented=true), nodes, lines)
end

export static_line_network_dynamics

"""
Documentation!!
"""
@with_kw struct static_line_network_dynamics
    nodes!
    lines!
    s_e
    t_e
    l_e
    l_e_int
    l_s
    l_t
    len_l
    len_n
end


function (d::static_line_network_dynamics)(dx, x, p, t)
    for i in 1:d.len_l
        d.lines![i](d.l_e[i], x[d.s_e[i]], x[d.t_e[i]], p, t)
    end
    for i in 1:d.len_n
        d.nodes![i](view(dx,i), view(x,i), d.l_s[i], d.l_t[i], p, t)
    end
    nothing
end

function static_line_network_dynamics(nodes!, lines!, s_e, t_e)
    len_l = length(lines!)
    len_n = length(nodes!)

    # Will get longer once we have more variables per line
    l_e_int = rand(len_l)

    # This will be views to more than one variable eventually
    l_e = [
        view(l_e_int, i)
        for i in 1:len_l ]

    # This list contains the indices of the variables belonging to the nodes
    #vorher x_idxs=[i for i in 1:len_n] gab einen 2d array raus

    # Create an array of views into the lines that selects only the sources
    # for each node
    l_s = [
         [l_e[j] for j in 1:len_l if s_e[j] == i]
        for i in 1:len_n ]
    # Create an array of views into the lines that selects only the targets
    # for each node
    l_t = [
        [l_e[j] for j in 1:len_l if t_e[j] == i]
        for i in 1:len_n ]

    static_line_network_dynamics(nodes!, lines!, s_e, t_e, l_e, l_e_int, l_s, l_t, len_l, len_n)
end

"""
When called with a graph, we construct the source and target vectors.
"""
function static_line_network_dynamics(nodes!, lines!, g::AbstractGraph)
    s_e = [src(e) for e in edges(g)]
    t_e = [dst(e) for e in edges(g)]
    static_line_network_dynamics(nodes!, lines!, s_e, t_e)
end


end # module
