module NetworkDynamics

struct diffusive_network_dynamics{T}
    """
    Implements the ODE
    ``\\frac{dx_i}{dt} = nodes(x) + A \\cdot x``
    """
    L::AbstractArray{T,2}
    nodes::Function
end
function (dnd::diffusive_network_dynamics)(dx, x, p, t)
    mul!(dx, dnd.L, x)
    @. dx = dx + dnd.nodes(x)
    nothing
end

end # module
