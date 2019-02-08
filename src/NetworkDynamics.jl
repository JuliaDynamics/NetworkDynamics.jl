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
    mul!(dx, dnd.L, x)
    @. dx = dx + dnd.nodes(dx, x, p, t)
    nothing
end

"""
When called with a graph, the dynamics defaults to using the laplacian.
"""
function diffusive_network_dynamics(g::AbstractGraph, nodes)
    diffusive_network_dynamics(laplacian_matrix(g), nodes)
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
