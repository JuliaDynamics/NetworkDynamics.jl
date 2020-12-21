module Utilities

using NLsolve
using LinearAlgebra
using SparseArrays

export maybe_idx, p_v_idx, p_e_idx, @nd_threads, construct_mass_matrix, warn_parallel, checkbounds_p


function warn_parallel(b::Bool)
    if b
        haskey(ENV, "JULIA_NUM_THREADS") &&
        parse(Int, ENV["JULIA_NUM_THREADS"]) > 1 ? nothing :
        print("Warning: You are using multi-threading with only one thread ",
        "available to Julia. Consider re-starting Julia with the environment ",
        "variable JULIA_NUM_THREADS set to the number of physical cores of your CPU.")
    end
end

"""
nd_threads: Allows control over threading by the 1st argument,
a boolean (likely from the network object).
This allows you to control for multithreading at runtime without
code duplication.
"""

macro nd_threads(trigger,args...)
    na = length(args)
    if na != 1
        throw(ArgumentError("wrong number of arguments in @nd_threads"))
    end
    ex = args[1]
    if ex.head == :for
        # Need to escape the bounds in the loop, but not the whole expression,
        # in order to work with the base Threads.@thread macro.
        # That macro does it's own escaping, but also it's own expression
        # checking (the head of the first expression must be :for. So we will
        # escape the bounds here, keeping the expression formed in a suitable
        # way for @threads

        # make sure we have the form :for, iterator, block
        @assert ex.args[1].head === :(=)
        @assert ex.args[2].head === :block
        @assert length(ex.args) == 2

        # extract the loop iterator
        loop_iter = Expr(:(=), esc.(ex.args[1].args)...)
        # extract the loop body and annotate @inbounds
        loop_body_ib = :(@inbounds $(esc(ex.args[2])))
        # build new loop
        loop = Expr(:for, loop_iter, loop_body_ib)

        return quote
            if $(esc(trigger))
                Threads.@threads $loop
            else
                @inbounds $(esc(ex))
            end
        end
    else
        throw(ArgumentError("unrecognized argument to @nd_threads"))
    end
end


"""
Utility function that drops the indexing operation when the argument is not
a subtype of AbstractArray. Used in the inner loop of Network Dynamics. Should
eventually be replaced by a macro that writes out the dispatches.
"""
@inline Base.@propagate_inbounds function maybe_idx(p::T, i) where T <: AbstractArray
    p[i]
end

@inline function maybe_idx(p, i)
    p
end

## non-allocating but code duplication

@inline Base.@propagate_inbounds function p_v_idx(p::Tuple{T1,T2,T3}, i) where {T1 <: AbstractArray, T2, T3}
    @view p[1][:, i]
end

@inline Base.@propagate_inbounds function p_v_idx(p::Tuple{T1,T2,T3}, i) where {T1 <: AbstractArray{T4, 1} where T4, T3, T2}
    p[1][i]
end

@inline Base.@propagate_inbounds function p_v_idx(p::Tuple{T1,T2, T3}, i) where {T1, T2, T3}
    p[1]
end

@inline Base.@propagate_inbounds function p_v_idx(p::Tuple{T1,T2}, i) where {T1 <: AbstractArray, T2}
    @view p[1][:, i]
end

@inline Base.@propagate_inbounds function p_v_idx(p::Tuple{T1,T2}, i) where {T1 <: AbstractArray{T3, 1} where T3, T2}
    p[1][i]
end

@inline Base.@propagate_inbounds function p_v_idx(p::Tuple{T1,T2}, i) where {T1, T2}
    p[1]
end


@inline function p_v_idx(p, i)
    p
end



## non-allocating but code duplication

@inline Base.@propagate_inbounds function p_e_idx(p::Tuple{T1,T2,T3}, i) where {T1 <: AbstractArray, T2, T3}
    @view p[2][:, i]
end

@inline Base.@propagate_inbounds function p_e_idx(p::Tuple{T1,T2,T3}, i) where {T1 <: AbstractArray{T4, 1} where T4, T3, T2}
    p[2][i]
end

@inline Base.@propagate_inbounds function p_e_idx(p::Tuple{T1,T2, T3}, i) where {T1, T2, T3}
    p[2]
end

@inline Base.@propagate_inbounds function p_e_idx(p::Tuple{T1,T2}, i) where {T1, T2 <: AbstractArray}
    @view p[2][:, i]
end

@inline Base.@propagate_inbounds function p_e_idx(p::Tuple{T1,T2}, i) where {T1, T2 <: AbstractArray{T3, 1} where T3}
    p[2][i]
end

@inline Base.@propagate_inbounds function p_e_idx(p::Tuple{T1,T2}, i) where {T1, T2}
    p[2]
end


@inline function p_e_idx(p, i)
    p
end


@inline function checkbounds_p(p, nv, ne)
    nothing
end

@inline function checkbounds_p(p::T, nv, ne) where T <: Tuple
    if p[1] isa AbstractArray
        size(p[1])[end] == nv ? nothing : error("Error: The size of the parameter array does not match the number of nodes. Make sure the correct number of parameters is given when using the tuple syntax.")
    end
    if p[2] isa AbstractArray
        size(p[2])[end] == ne ? nothing : error("Error: The size of the parameter array does not match the number of edges. Make sure the correct number of parameters is given when using the tuple syntax.")
    end
    nothing
end



export sum_coupling!

"""
A small utility function for writing diffusion dynamics. It provides the
 sum of all incoming edges.
"""
@inline function sum_coupling!(e_sum, in_edges)
    @inbounds for e in in_edges
        e_sum .+= e
    end
    nothing
end

export find_fixpoint

"""
Utility function for finding fixpoints.
"""

function find_fixpoint(nd, p, initial_guess)
    nl_res = nlsolve((dx, x) -> nd(dx, x, p, 0.), initial_guess)
    if converged(nl_res) == true
        return nl_res.zero
    else
        println("Failed to find fixpoint. Algorithm did not converge.")
    end
end

export RootRhs

"""
A utility function that provides a root function, to find valid initial
conditions. The operator ``M^\\dagger M - 1`` projects onto the kernel of ``M``.
The differential equation is ``M dx/dt = f(x)``, and can thus only be satisfied
if ``f(x)`` is zero on the kernel of ``M``. We look for such initial conditions
by returning the projection of ``f(x)`` onto the kernel as residue.

The use of the pseudoinverse here means this does not play nice with sparse
matrices. Beware when using for large systems.
"""
struct RootRhs
    rhs
    mpm
end
function (rr::RootRhs)(x)
    f_x = similar(x)
    rr.rhs(f_x, x, nothing, 0.)
    rr.mpm * f_x .- f_x
end

function RootRhs(of)
    mm = of.mass_matrix
    @assert mm !== nothing
    mpm = pinv(mm) * mm
    RootRhs(of.f, mpm)
end


export find_valid_ic


"""
Try to find valid initial conditions for problems involving mass matrices.
Uses RootRhs as residual function.
"""
function find_valid_ic(of, ic_guess)
    rr = RootRhs(of)
    nl_res = nlsolve(rr, ic_guess)
    if converged(nl_res) == true
        return nl_res.zero
    else
        println("Failed to find initial conditions on the constraint manifold!")
        println("Try running nlsolve with other options.")
    end
end


export syms_containing, idx_containing

"""
Find all symbols present in a network_dynamics that contain the string or symbol
`str`
"""
function syms_containing(nd, str)
    [s for s in nd.syms if occursin(string(str), string(s))]
end


"""
Find all indices of variables with symbols containing the string or symbol
`str`
"""
function idx_containing(nd, str)
    [i for (i, s) in enumerate(nd.syms) if occursin(string(str), string(s))]
end


function construct_mass_matrix(mmv_array, gs)
    if all([mm == I for mm in mmv_array])
        mass_matrix = I
    else
        mass_matrix = sparse(1.0I,gs.dim_v, gs.dim_v)
        for (i, mm) in enumerate(mmv_array)
            ind = gs.v_idx[i]
            if ndims(mm) == 0
                copyto!(@view(mass_matrix[ind, ind]), mm*I)
            elseif ndims(mm) == 1
                copyto!(@view(mass_matrix[ind, ind]), Diagonal(mm))
            elseif ndims(mm) == 2 # ndims(I) = 2
                # `I` does not support broadcasting but copyto! combined with views
                copyto!(@view(mass_matrix[ind, ind]), mm)
            else
                error("The mass matrix needs to be interpretable as a 2D matrix.")
            end
        end
    end
    mass_matrix
end

function construct_mass_matrix(mmv_array, mme_array, gs)
    if all([mm == I for mm in mmv_array]) && all([mm == I for mm in mme_array])
        mass_matrix = I
    else
        dim_nd = gs.dim_v + gs.dim_e
        mass_matrix = sparse(1.0I,dim_nd,dim_nd)
        for (i, mm) in enumerate(mmv_array)
            ind = gs.v_idx[i]
            if ndims(mm) == 0
                copyto!(@view(mass_matrix[ind, ind]), mm*I)
            elseif ndims(mm) == 1
                copyto!(@view(mass_matrix[ind, ind]), Diagonal(mm))
            elseif ndims(mm) == 2 # ndims(I) = 2
                # `I` does not support broadcasting but copyto!
                copyto!(@view(mass_matrix[ind, ind]), mm)
            else
                error("The mass matrix needs to be interpretable as a 2D matrix.")
            end
        end
        for (i, mm) in enumerate(mme_array)
            ind = gs.dim_v .+ (gs.e_idx[i])
            if ndims(mm) == 0
                copyto!(@view(mass_matrix[ind, ind]), mm*I)
            elseif ndims(mm) == 1
                copyto!(@view(mass_matrix[ind, ind]), Diagonal(mm))
            elseif ndims(mm) == 2 # ndims(I) = 2
                # `I` does not support broadcasting but copyto!
                copyto!(@view(mass_matrix[ind, ind]), mm)
            else
                error("The mass matrix needs to be interpretable as a 2D matrix.")
            end
        end
    end
    mass_matrix
end

end #module
