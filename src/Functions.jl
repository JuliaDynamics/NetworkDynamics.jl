module NDFunctions
#=
This module needs to be documented. Together with Constructors
it forms the backbone of the core API, users should provide us
with Arrays of VertexFunction and EdgeFunction as well as a graph and that's it.
=#

struct StaticVertex #???
    f! # ToDo
    dim # number of dimensions of x
    sym # Symbols for the dimensions
end

struct StaticEdge
    f! # (l, v_s, v_t, p, t) -> nothing
    dim # number of dimensions of x
    sym # Symbols for the dimensions
end

struct ODEVertex
    f! # The function with signature (dx, x, e_s, e_t, p, t) -> nothing
    dim # number of dimensions of x
    massmatrix # Mass matrix for the equation
    sym # Symbols for the dimensions
end

function ODEVertex(f!, dim)
    ODEVertex(f!, dim, massmatrix=sparse(1.0I, dim, dim), nothing)
end

struct ODEEdge
    f! # The function with signature (dx, x, e_s, e_t, p, t) -> nothing
    dim # number of dimensions of x
    massmatrix # Mass matrix for the equation
    sym # Symbols for the dimensions
end

function ODEEdge(f!, dim)
    ODEEdge(f!, dim, massmatrix=sparse(1.0I, dim, dim), [:e for i in 1:dim])
end

const VertexFunction = Union{ODEVertex, StaticVertex}
const EdgeFunction = Union{ODEEdge, StaticEdge}

function edge_constraint!(f!, de, e, v_s, v_t, p, t)
    # If mass matrix = 0 the differential equation sets de = 0.
    # To set e to the value calculated by f! we first write the value calculated
    # by f! into de, the subtract e. This leads to the  constraint
    # 0 = - e + f(v_s, v_p, p, t)
    # where f(...) denotes the value that f!(a, ...) writes into a.
    f!(de, v_s, v_t, p, t)
    de .-= e
    nothing
end

function ODEEdge(se::StaticEdge)
    ODEEdge(
    edge_constraint!(se.f!, de, e, v_s, v_t, p, t),
    massmatrix = 0., # should be zero(T)
    se.dim,
    se.sym)
end

# investigate whether convert(::Type(ODEEdge), se::StaticEdge) = ODEEdge(se)
# is useful, check out PowerDynamics for what it does...
end
