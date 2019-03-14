module NDFunctions

const VertexFunction = Union{ODEVertex, StaticVertex, DDEVertex}

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
    massmatrix # Mass matrix for the equation
    dim # number of dimensions of x
    sym # Symbols for the dimensions
end

struct ODEEdge
    f! # The function with signature (dx, x, e_s, e_t, p, t) -> nothing
    massmatrix # Mass matrix for the equation
    dim # number of dimensions of x
    sym # Symbols for the dimensions
end

struct DDEVertex
    f! # The function with signature (dx, x, h, e_s, h_e_s, e_t, h_e_t, p, t) -> nothing
    massmatrix # Mass matrix for the equation
    dim # number of dimensions of x
    sym # Symbols for the dimensions
    taus # List of delays
end

struct DDEEdge{T}
    f! # The function with signature (dx, x, h, e_s, h_e_s, e_t, h_e_t, p, t) -> nothing
    massmatrix # Mass matrix for the equation
    dim # number of dimensions of x
    sym # Symbols for the dimensions
    taus # list of delays
end

# SDEEdge, SDEVertex

function DDEVertex(ov::ODEVertex)
    DDEVertex(
    (dx, x, h, e_s, h_e_s, e_t, h_e_t, p, t) -> ov.f!(dx, x, e_s, e_t, p, t),
    ov.massmatrix # Mass matrix for the equation
    ov.dim,
    ov.sym,
    [])
end

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
    massmatrix = 0. # should be zero(T)
    se.dim,
    se.sym)
end

# investigate whether convert(::Type(ODEEdge), se::StaticEdge) = ODEEdge(se)
# is useful, check out PowerDynamics for what it does...
end
