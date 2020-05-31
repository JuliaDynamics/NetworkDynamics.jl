module Utilities

using NLsolve
using LinearAlgebra

export maybe_idx, p_v_idx, p_e_idx, @nd_threads


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
        # Need to escape the bounds of the loop, but not the whole expression,
        # in order to work with the base Threads.@thread macro.
        # That macro does it's own escaping, but also it's own expression
        # checking (the head of the first expression must be :for. So we will
        # escape the bounds here, keeping the expression formed in a suitable
        # way for @threads
        #
        # This is a brute-force drill down an AST of a for loop to the limits
        # of said loop. See: dump(quote for i in e.a:e.b print(i) end end)
        ex.args[1].args[1] = esc(ex.args[1].args[1])
        ex.args[1].args[2].args[2] = esc(ex.args[1].args[2].args[2])
        ex.args[1].args[2].args[3] = esc(ex.args[1].args[2].args[3])

        # Same for the body of the loop
        ex.args[2] = esc(ex.args[2])

        return quote
            @inbounds if $(esc(trigger))
                Threads.@threads $ex
            else
                $(ex)
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

## allocating - fix later
@inline Base.@propagate_inbounds function p_v_idx(p::Tuple{T1,T2,T3}, i) where {T1,T2,T3}
    p_v_idx(p[1:2],i)
end

@inline Base.@propagate_inbounds function p_v_idx(p::Tuple{T1,T2}, i) where {T1 <: AbstractArray, T2}
    p[1][:, i]
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



## allocating - fix later
@inline Base.@propagate_inbounds function p_e_idx(p::Tuple{T1,T2,T3}, i) where {T1,T2,T3}
    p_e_idx((p[1], p[2]), i)
end

@inline Base.@propagate_inbounds function p_e_idx(p::Tuple{T1,T2}, i) where {T1, T2 <: AbstractArray}
    p[2][:, i]
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


export oriented_edge_sum!

"""
A small utility function for writing diffusion dynamics. It provides the
oriented sum of all the incident edges.
"""
@inline function oriented_edge_sum!(e_sum, e_s, e_d)
    @inbounds for e in e_s
        e_sum .-= e
    end
    @inbounds for e in e_d
        e_sum .+= e
    end
    nothing
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
    @assert mm != nothing
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


end #module
