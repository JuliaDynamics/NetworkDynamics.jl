

@Base.kwdef struct StaticVertex{T} <: VertexFunction
    f!::T # The function with signature (dx, x, edges, h, p, t) -> nothing
    dim::Int # number of dimensions of x
    mass_matrix=I # Mass matrix for the equation
    sym=[:v for i in 1:dim] # Symbols for the dimensions
    jac::bool
    VertexJacobian!::F1
    z::AbstractArray

    function JacVertex(f!::T, # The function with signature (dx, x, edges, h, p, t) -> nothingdim::Int, # number of dimensions of x
                       mass_matrix=I, # Mass matrix for the equation
                       sym=[:v for i in 1:dim], # Symbols for the dimensions
                       VertexJacobian!::F1.
                       z:AbstractArray)
        #J = internal_jacobian ...
end



@Base.kwdef struct StaticEdge{T} <: EdgeFunction
    f!::T # (e, v_s, v_t, p, t) -> nothing
    dim::Int # number of dimensions of e
    coupling = :undefined # :directed, :symmetric, :antisymmetric, :fiducial, :undirected
    sym=[:e for i in 1:dim] # Symbols for the dimensions
    jac::bool
    EdgeJacobian!::F


    function JacEdge(user_f!::T,
                           dim::Int,
                           coupling::Symbol,
                           sym::Vector{Symbol},
                           EdgeJacobian!::F) where T,

        coupling_types = (:undefined, :directed, :fiducial, :undirected, :symmetric,
                          :antisymmetric)

        coupling ∈ coupling_types ? nothing :
            error("Coupling type not recognized. Choose from $coupling_types.")

        dim > 0 ? nothing : error("dim has to be a positive number.")

        dim == length(sym) ? nothing : error("Please specify a symbol for every dimension.")

        if coupling ∈ [:undefined, :directed]
            return new{T}(user_f!, dim, coupling, sym)

        elseif coupling == :fiducial
            dim % 2 == 0 ? nothing : error("Fiducial edges are required to have even dim.
                                            The first dim args are used for src -> dst,
                                            the second for dst -> src coupling.")
            return new{T}(user_f!, dim, coupling, sym)

        elseif coupling == :undirected
            # This might cause unexpected behaviour if source and destination vertex don't
            # have the same internal arguments.
            # Make sure to explicitly define the edge is :fiducial in that case.
            f! = @inline (e, v_s, v_d, p, t) -> begin
                @inbounds user_f!(view(e,1:dim), v_s, v_d, p, t)
                @inbounds user_f!(view(e,dim+1:2dim), v_d, v_s, p, t)
                nothing
            end
        elseif coupling == :antisymmetric
            f! = @inline (e, v_s, v_d, p, t) -> begin
                @inbounds user_f!(view(e,1:dim), v_s, v_d, p, t)
                @inbounds for i in 1:dim
                    e[dim + i] = -1.0 * e[i]
                end
                nothing
            end
        elseif coupling == :symmetric
            f! = @inline (e, v_s, v_d, p, t) -> begin
                @inbounds user_f!(view(e,1:dim), v_s, v_d, p, t)
                @inbounds for i in 1:dim
                    e[dim + i] = e[i]
                end
                nothing
            end
        end
        J_s = source_jacobian(v_s, v_d, p, t)
        J_d = edge_jacobian(v_s, v_d, p, t)
        # For edges with mass matrix this will be a little more complicated
        return new{typeof(f!)}(f!, 2dim, coupling, repeat(sym, 2))
    end
end
