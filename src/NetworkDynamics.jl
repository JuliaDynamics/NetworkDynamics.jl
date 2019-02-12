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
    @. dx = dx + dnd.nodes(dx, x, p, t) # ToDo test that this does the right thing
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
    Placeholder
end
function (dnd::network_dynamics)(dx, x, p, t)
    nothing
end



export network_dynamics_on_the_line

"""
edge_list: list of edges
nodes: scalar function ``x \\rightarrow nodes(x)``
lines: scalar function ``x \\rightarrow lines(x)``

Implements the ODE
``\\frac{dx_i}{dt} = nodes(x_i) + \\sum_j A_{ij} lines(x_i - x_j)``

where ``A`` is the adjacency matrix of ``B``
"""
@with_kw struct network_dynamics_on_the_line{T}
    B::AbstractArray{T,2}
    nodes::Function
    lines::Function
end
function (dnd::network_dynamics_on_the_line)(dx, x, p, t)
    @. dx = dnd.nodes(x)
    nothing
end

"""
When called with a graph, the dynamics defaults to using the laplacian.
"""
function network_dynamics_on_lines(g::AbstractGraph, nodes)
    diffusive_network_dynamics(incidence_matrix(g, oriented=true), nodes)
end



end # module
