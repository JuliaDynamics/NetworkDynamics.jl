# Together with Constructors this file forms the backbone of the core API.
# It provide the basic types to construct Arrays of VertexFunction and
# EdgeFunction which can be handled by network_dynamics.


using LinearAlgebra

import Base.convert
import Base.promote_rule


export VertexFunction
export EdgeFunction
export StaticVertex
export StaticEdge
export ODEVertex
export ODEEdge
export DDEVertex
export StaticDelayEdge

"""
Abstract supertype for all vertex functions.
"""
abstract type VertexFunction end
"""
Abstract supertype for all edge functions.
"""
abstract type EdgeFunction end


"""
    Helper function to check correct input args for ComponentFunctions
"""
@inline function _dimcheck(dim, sym)
    dim ≤ 0 && error("dim has to be a positive number.")
    dim != length(sym) && error("Please specify a symbol for every dimension.")
    return nothing
end

"""
    Helper function to check correct input args for ComponentFunctions
"""
@inline function _argcheck(f,n)
    !((n+1) in getproperty.(methods(f), :nargs)) && error("f does not take the right number of arguments ($n args)")
    return nothing
end


"""
    StaticEdge(; f, dim, coupling, sym)

Wrapper that ensures compatibility of a **mutating** function `f` with
the key constructor `network_dynamics`.

`f`  describes the local behaviour at a static edge and has to respect
the following calling syntax

```julia
f(e, v_s, v_d, p, t) -> nothing
```

Here  `e`, `p` and `t` are the usual arguments, while
`v_s` and `v_d` are arrays containing the vertices which are
the source and destination of the described edge.

  - `dim` is the number of independent variables in the edge equations and
  - `sym` is an array of symbols for these variables.
  - `coupling` is a Symbol describing if the EdgeFunction is intended for a directed graph (`:directed`) or for an undirected graph (`{:undirected, :symmetric, :antisymmetric, :fiducial}`). `:directed` is intended for directed graphs. `:undirected` is the default option and is only compatible with SimpleGraph.  In this case f should specify the coupling from a source vertex to a destination vertex. `:symmetric` and `:antisymmetric` trigger performance optimizations, if `f` has that symmetry property. `:fiducial` lets the user specify both the coupling from src to dst, as well as the coupling from dst to src and is intended for advanced users, i.e. the edge gets passed a vector of `2dim` edge states, where the first `dim` elements will be presented to the dst and the second `dim` elements will be presented to src,

Optional metadata keywords:
  - `name` string representation of model
  - `psym` vector of symbols to name parameters
  - `obsf` function `g(e, v_s, v_d, p, t) -> obs` which returns a vector of virtual/observed states
  - `obssym` vector of symbols to describe order of observed states

For more details see the documentation.
"""
Base.@kwdef struct StaticEdge{T,G} <: EdgeFunction
    f::T # (e, v_s, v_t, p, t) -> nothing
    dim::Int # number of dimensions of e
    coupling = :undefined # :directed, :symmetric, :antisymmetric, :fiducial, :undirected
    sym = [:e for i in 1:dim] # Symbols for the dimensions
    name::String = "StaticEdge"          # optional: name of edge
    psym::Vector{Symbol} = Symbol[]   # optional: symbols of parameters in order
    obsf::G = nothing                 # optional: observables
    obssym::Vector{Symbol} = Symbol[] # optional: symbols of observables

    function StaticEdge(user_f::T,
                        dim::Int,
                        coupling::Symbol,
                        sym::Vector{Symbol},
                        name="StaticEdge", psym=Symbol[], obsf::G=nothing, obssym=Symbol[]) where {T,G}

        coupling_types = (:undefined, :directed, :fiducial, :undirected, :symmetric,
                          :antisymmetric)

        coupling ∉ coupling_types && error("Coupling type not recognized. Choose from $coupling_types.")

        _dimcheck(dim, sym)
        _argcheck(user_f, 5)

        if coupling ∈ [:undefined, :directed]
            return new{T,G}(user_f, dim, coupling, sym, name, psym, obsf, obssym)

        elseif coupling == :fiducial
            dim % 2 != 0 && error("Fiducial edges are required to have even dim. ",
                                  "The first dim args are used for src -> dst ",
                                  "the second for dst -> src coupling.")
            return new{T,G}(user_f, dim, coupling, sym, name, psym, obsf, obssym)

        elseif coupling == :undirected
            # This might cause unexpected behaviour if source and destination vertex don't
            # have the same internal arguments.
            # Make sure to explicitly define the edge is :fiducial in that case.
            f = @inline (e, v_s, v_d, p, t) -> begin
                @inbounds user_f(view(e, 1:dim), v_s, v_d, p, t)
                @inbounds user_f(view(e, dim+1:2dim), v_d, v_s, p, t)
                nothing
            end
        elseif coupling == :antisymmetric
            f = @inline (e, v_s, v_d, p, t) -> begin
                @inbounds user_f(view(e, 1:dim), v_s, v_d, p, t)
                @inbounds for i in 1:dim
                    e[dim+i] = -1.0 * e[i]
                end
                nothing
            end
        elseif coupling == :symmetric
            f = @inline (e, v_s, v_d, p, t) -> begin
                @inbounds user_f(view(e, 1:dim), v_s, v_d, p, t)
                @inbounds for i in 1:dim
                    e[dim+i] = e[i]
                end
                nothing
            end
        end
        # For edges with mass matrix this will be a little more complicated
        return new{typeof(f),G}(f, 2dim, coupling, _double_syms(sym), name, psym, obsf, obssym)
    end
end




"""
    ODEVertex(; f, dim, mass_matrix, sym)

Wrapper that ensures compatibility of a **mutating** function **`f`** with
the key constructor `network_dynamics`.

**`f`**  describes the local behaviour at a dynamic node and has to respect
the following calling syntax

```julia
f(dv, v, edges, p, t) -> nothing
```

Here `dv`, `v`, `p` and `t` are the usual ODE arguments, while
`edges` is an Array containing the edges for which the vertex is the destination (in-edges for directed graphs).

**`dim`** is the number of independent variables in the vertex equations and
**`sym`** is an array of symbols for these variables.
**`mass_matrix`** is an optional argument that defaults to the identity
matrix `I`. If a mass matrix M is given the system `M * dv = f` will be
solved.

Optional metadata keywords:
  - `name` string representation of model
  - `psym` vector of symbols to name parameters
  - `obsf` function `g(v, edges, p, t) -> obs` which returns a vector of virtual/observed states
  - `obssym` vector of symbols to describe order of observed states

For more details see the documentation.
"""
Base.@kwdef struct ODEVertex{T,G} <: VertexFunction
    f::T # signature (dx, x, edges, p, t) -> nothing
    dim::Int
    mass_matrix = I
    sym = [:v for i in 1:dim]
    name::String = "ODEVertex"        # optional: name of vertex
    psym::Vector{Symbol} = Symbol[]   # optional: symbols of parameters in order
    obsf::G = nothing                 # optional: observable
    obssym::Vector{Symbol} = Symbol[] # optional: symbols of observables

    function ODEVertex(f, dim, mass_matrix, sym,
                       name="ODEVertex", psym=Symbol[], obsf=nothing, obssym=Symbol[])
        _dimcheck(dim, sym)
        _argcheck(f, 5)
        return new{typeof(f), typeof(obsf)}(f, dim, mass_matrix, sym, name, psym, obsf, obssym)
    end
end

"""
    StaticVertex(; f, dim, sym)

Wrapper for ODEVertex with 0 mass matrix, i.e. static behaviour / algebraic constraint in mass matrix form.

**`f`**  describes the local behaviour at a static node and has to respect
the following calling syntax

```julia
f(v, edges, p, t) -> nothing
```

Here  `v`, `p` and `t` are the usual arguments, while
`edges` is an arrays containing the edges for which the
described vertex is the destination (in-edges for directed graphs).

`dim` is the number of independent variables in the vertex equations and
`sym` is an array of symbols for these variables.

Optional metadata keywords:
  - `name` string representation of model
  - `psym` vector of symbols to name parameters
  - `obsf` function `g(v, edges, p, t) -> obs` which returns a vector of virtual/observed states
  - `obssym` vector of symbols to describe order of observed states

For more details see the documentation.
"""
function StaticVertex(; f, dim::Int, kwargs...)
    # If mass matrix = 0 the differential equation sets dx = 0.
    # To set x to the value calculated by f we first write the value calculated
    # by f into dx, then subtract x. This leads to the  constraint
    # 0 = - x + f(...)
    # where f(...) denotes the value that f(a, ...) writes into a.
    _argcheck(f, 4)
    _f = (dx, x, edges, p, t) -> begin
        f(dx, edges, p, t)
        @inbounds for i in eachindex(dx)
            dx[i] = dx[i] - x[i]
        end
        return nothing
    end
    return ODEVertex(; f=_f, dim, mass_matrix=0.0, kwargs...)
end


"""
    ODEEdge(; f, dim, mass_matrix, sym)

Wrapper that ensures compatibility of a **mutating** function **`f`** with
the key constructor `network_dynamics`.

**`f`**  describes the local behaviour at a dynamic edge and has to respect
the following calling syntax

```julia
f(de, e, v_s, v_d, p, t) -> nothing
```

Here  `de`, `e`, `p` and `t` are the usual arguments, while
`v_s` and `v_d` are arrays containing the vertices which are
the source and destination of the described edge.

**`dim`** is the number of independent variables in the edge equations and
**`sym`** is an array of symbols for these variables. For more details see
the documentation.
**`mass_matrix`** is an optional argument that defaults to the identity
matrix `I`. If a mass matrix M is given the system `M * de = f` will be
solved.

Optional metadata keywords:
  - `name` string representation of model
  - `psym` vector of symbols to name parameters
  - `obsf` function `g(e, v_s, v_d, p, t) -> obs` which returns a vector of virtual/observed states
  - `obssym` vector of symbols to describe order of observed states

For more details see the documentation.
"""
Base.@kwdef struct ODEEdge{T,G} <: EdgeFunction
    f::T # (de, e, v_s, v_t, p, t) -> nothing
    dim::Int # number of dimensions of e
    coupling = :undefined # :directed, :symmetric, :antisymmetric, :fiducial, :undirected
    mass_matrix = I # Mass matrix for the equation
    sym = [:e for i in 1:dim] # Symbols for the dimensions
    name::String = "ODEEdge"          # optional: name of edge
    psym::Vector{Symbol} = Symbol[]   # optional: symbols of parameters in order
    obsf::G = nothing                 # optional: observables
    obssym::Vector{Symbol} = Symbol[] # optional: symbols of observables

    function ODEEdge(user_f::T,
                     dim::Int,
                     coupling::Symbol,
                     mass_matrix,
                     sym::Vector{Symbol},
                     name="ODEEdge", psym=Symbol[], obsf::G=nothing, obssym=Symbol[]) where {T,G}

        coupling_types = (:directed, :fiducial, :undirected)

        coupling == :undefined && error("ODEEdges with undefined coupling type are not implemented at the " *
                                        "moment. Choose `coupling` from $coupling_types.")
        coupling ∈ (:symmetric, :antisymmetric) && error("Coupling type $coupling is not available for ODEEdges.")
        coupling ∉ coupling_types && error("Coupling type not recognized. Choose from $coupling_types.")

        _argcheck(user_f, 6)
        _dimcheck(dim, sym)

        if coupling == :directed
            return new{T,G}(user_f, dim, coupling, mass_matrix, sym, name, psym, obsf, obssym)
        elseif coupling == :fiducial
            dim % 2 != 0 && error("Fiducial edges are required to have even dim. ",
                                  "The first dim args are used for src -> dst ",
                                  "the second for dst -> src coupling.")
            return new{T,G}(user_f, dim, coupling, mass_matrix, sym, name, psym, obsf, obssym)

        elseif coupling == :undirected
            # This might cause unexpected behaviour if source and destination vertex don't
            # have the same internal arguments.
            # Make sure to explicitly define the edge is :fiducial in that case.
            f = @inline (de, e, v_s, v_d, p, t) -> begin
                @inbounds user_f(view(de, 1:dim), view(e, 1:dim), v_s, v_d, p, t)
                @inbounds user_f(view(de, dim+1:2dim), view(e, dim+1:2dim), v_d, v_s, p, t)
                nothing
            end
            let M = mass_matrix
                if M === I
                    newM = M
                elseif M isa Number
                    newM = M
                elseif M isa Vector
                    newM = repeat(M, 2)
                elseif M isa Matrix
                    newM = [M zeros(size(M)); zeros(size(M)) M]
                end
                return new{typeof(f),G}(f, 2dim, coupling, newM, _double_syms(sym), name, psym, obsf, obssym)
            end
        end
    end
end


"""
    DDEVertex(; f, dim, mass_matrix, sym)

Wrapper that ensures compatibility of a **mutating** function **`f`** with
the key constructor `network_dynamics`.

**`f`**  describes the local behaviour at a dynamic node and has to respect
the following calling syntax

```julia
f(dv, v, edges, h, p, t) -> nothing
```

Here `dv`, `v`, `p` and `t` are the usual ODE arguments, while
`edges` is an Arry of incoming edges. `h` is the history array for `v`.

  - `dim` is the number of independent variables in the edge equations and
  - `sym` is an array of symbols for these variables.
  - `mass_matrix` is an optional argument that defaults to the identity
    matrix `I`. If a mass matrix M is given the system `M * de = f` will be
    solved.

For more details see the documentation.
"""
Base.@kwdef struct DDEVertex{T} <: VertexFunction
    f::T # The function with signature (dx, x, edges, h, p, t) -> nothing
    dim::Int # number of dimensions of x
    mass_matrix = I # Mass matrix for the equation
    sym = [:v for i in 1:dim] # Symbols for the dimensions
    function DDEVertex(f::T, dim::Int, mass_matrix=I, sym=[:v for i in 1:dim]) where {T}
        _dimcheck(dim, sym)
        _argcheck(f, 6)
        return new{typeof(f)}(f, dim, mass_matrix, sym)
    end
end



"""
    StaticDelayEdge(; f, dim, coupling, sym)

Like a static edge but with extra arguments for the history of the source and destination vertices. This is NOT a DDEEdge.
"""
Base.@kwdef struct StaticDelayEdge{T} <: EdgeFunction
    f::T # (e, v_s, v_t, h_v_s, h_v_d, p, t) -> nothing
    dim::Int # number of dimensions of e
    coupling = :undefined # :directed, :symmetric, :antisymmetric, :fiducial, :undirected
    sym = [:e for i in 1:dim] # Symbols for the dimensions


    function StaticDelayEdge(user_f::T,
                             dim::Int,
                             coupling::Symbol,
                             sym::Vector{Symbol}) where {T}

        coupling_types = (:undefined, :directed, :fiducial, :undirected, :symmetric,
                          :antisymmetric)

        coupling ∉ coupling_types && error("Coupling type not recognized. Choose from $coupling_types.")
        _argcheck(user_f, 7)
        _dimcheck(dim, sym)

        if coupling ∈ [:undefined, :directed]
            return new{T}(user_f, dim, coupling, sym)

        elseif coupling == :fiducial
            dim % 2 != 0 && error("Fiducial edges are required to have even dim. ",
                                  "The first dim args are used for src -> dst ",
                                  "the second for dst -> src coupling.")
            return new{T}(user_f, dim, coupling, sym)
        elseif coupling ∈ (:antisymmetric, :symmetric)
            error("$coupling coupling not implemented for edges with delay. If you need it please open an issue on GitHub.")
        elseif coupling == :undirected
        # This might cause unexpected behaviour if source and destination vertex don't
        # have the same internal arguments.
        # Make sure to explicitly define the edge is :fiducial in that case.
            f = @inline (e, v_s, v_d, h_v_s, h_v_d, p, t) -> begin
                @inbounds user_f(view(e,1:dim), v_s, v_d, h_v_s, h_v_d, p, t)
                @inbounds user_f(view(e,dim+1:2dim), v_d, v_s, h_v_d, h_v_s, p, t)
                nothing
            end
        end
        return new{typeof(f)}(f, 2dim, coupling, _double_syms(sym))
    end
end

# Promotion rules



convert(::Type{ODEEdge}, x::StaticEdge) = ODEEdge(x)
promote_rule(::Type{ODEEdge}, ::Type{StaticEdge}) = ODEEdge


"""
    ODEEdge(se::StaticEdge)

Promotes a StaticEdge to an ODEEdge with zero mass matrix.
"""
function ODEEdge(se::StaticEdge)
    # If mass matrix = 0 the differential equation sets dx = 0.
    # To set x to the value calculated by f we first write the value calculated
    # by f into dx, then subtract x. This leads to the  constraint
    # 0 = - x + f(...)
    # where f(...) denotes the value that f(a, ...) writes into a.
    let _f=se.f, dim=se.dim, coupling=se.coupling, sym=se.sym, name=se.name, psym=se.psym, obsf=se.obsf, obssym=se.obssym

        f = (dx, x, v_s, v_d, p, t) -> begin
            _f(dx, v_s, v_d, p, t)
            @inbounds for i in eachindex(dx)
                dx[i] = dx[i] - x[i]
            end
            return nothing
        end

        if coupling == :undirected
            # undirected means the reconstruction has already happend and the edge may now
            # be considered fiducial
            return ODEEdge(f, dim, :fiducial, 0.0, sym)
        else
            return ODEEdge(f, dim, coupling, 0.0, sym)
        end
    end
end

syms(v::VertexFunction) = v.sym
syms(e::EdgeFunction)   = e.sym

observed_syms(v::VertexFunction) = v.obssym
observed_f(v::VertexFunction) = v.obsf

function _double_syms(syms)
    strs = string.(syms)
    Symbol.(vcat("dst_" .* strs, "src_" .* strs))
end
