import NLsolve: nlsolve, converged
using LinearAlgebra
using SparseArrays

"""
    warn_parallel(b::Bool)

If `b=true` (run in parallel mode) give a warning if either

  - `ENV["JULIA_NUM_THREADS"]` not found in ENV
  - `ENV["JULIA_NUM_THREADS"] ≤ 1`, i.e. there is just a single thread
"""
function warn_parallel(b::Bool)
    if b
        if !haskey(ENV, "JULIA_NUM_THREADS") || parse(Int, ENV["JULIA_NUM_THREADS"]) ≤ 1
            print("Warning: You are using multi-threading with only one thread ",
                  "available to Julia. Consider re-starting Julia with the environment ",
                  "variable JULIA_NUM_THREADS set to the number of physical cores of your CPU.")
        end
    end
end

"""
    @nd_threads trigger expr

nd_threads: Allows control over threading by the 1st argument,
a boolean (likely from the network object).
This allows you to control for multithreading at runtime without
code duplication.
"""
macro nd_threads(trigger, args...)
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
    maybe_idx(p::T, i) where T <: AbstractArray

Utility function that drops the indexing operation when the argument is not
a subtype of AbstractArray. Used in the inner loop of Network Dynamics. Should
eventually be replaced by a macro that writes out the dispatches.
"""
Base.@propagate_inbounds function maybe_idx(p::T, i) where {T<:AbstractArray}
    p[i]
end

@inline function maybe_idx(p, i)
    p
end

## non-allocating but code duplication

Base.@propagate_inbounds function p_v_idx(p::Tuple{T1,T2}, i) where {T1<:AbstractArray,T2}
    @view p[1][:, i]
end

Base.@propagate_inbounds function p_v_idx(p::Tuple{T1,T2}, i) where {T1<:AbstractVector{T3} where {T3},T2}
    p[1][i]
end

Base.@propagate_inbounds function p_v_idx(p::Tuple{T1,T2}, i) where {T1,T2}
    p[1]
end


@inline function p_v_idx(p, i)
    p
end



## non-allocating but code duplication

Base.@propagate_inbounds function p_e_idx(p::Tuple{T1,T2}, i) where {T1,T2<:AbstractArray}
    @view p[2][:, i]
end

Base.@propagate_inbounds function p_e_idx(p::Tuple{T1,T2}, i) where {T1,T2<:AbstractVector{T3} where {T3}}
    p[2][i]
end

Base.@propagate_inbounds function p_e_idx(p::Tuple{T1,T2}, i) where {T1,T2}
    p[2]
end


@inline function p_e_idx(p, i)
    p
end


@inline function checkbounds_p(p, nv, ne)
    nothing
end

@inline function checkbounds_p(p::T, nv, ne) where {T<:Tuple}
    if p[1] isa AbstractArray
        size(p[1])[end] != nv &&
            error("Error: The size of the parameter array does not match the number of nodes. Make sure the correct number of parameters is given when using the tuple syntax.")
    end
    if p[2] isa AbstractArray
        size(p[2])[end] != ne &&
            error("Error: The size of the parameter array does not match the number of edges. Make sure the correct number of parameters is given when using the tuple syntax.")
    end
    nothing
end



export sum_coupling!

"""
    sum_coupling!(e_sum, dst_edges)

A small utility function for writing diffusion dynamics. It provides the
sum of all incoming edges.
"""
@inline function sum_coupling!(e_sum, dst_edges)
    @inbounds for e in dst_edges
        e_sum .+= e
    end
    nothing
end

export find_fixpoint

"""
    find_fixpoint(nd, p, initial_guess)

Utility function for finding fixpoints.
"""
function find_fixpoint(nd, p, initial_guess)
    nl_res = nlsolve((dx, x) -> nd(dx, x, p, 0.0), initial_guess)
    if converged(nl_res) == true
        return nl_res.zero
    else
        println("Failed to find fixpoint. Algorithm did not converge.")
        println("Try running nlsolve with other options.")
    end
end

export RootRhs

"""
    struct RootRhs
    RootRhs(nd)

A utility function that provides a root function, to find valid initial
conditions. The operator ``M^\\dagger M - 1`` projects onto the kernel of ``M``.
The differential equation is ``M dx/dt = f(x)``, and can thus only be satisfied
if ``f(x)`` is zero on the kernel of ``M``. We look for such initial conditions
by returning the projection of ``f(x)`` onto the kernel as residue.

The use of the pseudoinverse here means this does not play nice with sparse
matrices. Beware when using for large systems.
"""
struct RootRhs{T,S}
    rhs::T
    mpm::S
end
function (rr::RootRhs)(x, p=nothing, t=nothing)
    f_x = similar(x)
    rr.rhs(f_x, x, p, t)
    return rr.mpm * f_x .- f_x
end

function RootRhs(nd)
    mm = nd.mass_matrix
    @assert mm !== nothing
    mpm = pinv(mm) * mm
    return RootRhs(nd.f, mpm)
end


export find_valid_ic


"""
    find_valid_ic(nd, ic_guess)

Try to find valid initial conditions for problems involving mass matrices.
Uses RootRhs as residual function.
"""
function find_valid_ic(nd, initial_guess; p=nothing)
    rr = RootRhs(nd)
    nl_res = nlsolve(x -> rr(x, p), initial_guess)
    if converged(nl_res) == true
        return nl_res.zero
    else
        println("Failed to find initial conditions on the constraint manifold!")
        println("Try running nlsolve with other options.")
    end
end


export syms_containing, idx_containing

"""
    syms_containing(nd, expr)

Find all symbols present in a network_dynamics that contain the string, regex or symbol `expr`.
"""
function syms_containing(nd, expr)
    [s for s in nd.syms if occursin(expr, string(s))]
end
syms_containing(nd, expr::Symbol) = syms_containing(nd, string(expr))

"""
    idx_containing(nd, expr)

Find all indices of variables with symbols containing the string, regex or symbol `expr`
"""
function idx_containing(nd, expr)
    [i for (i, s) in enumerate(nd.syms) if occursin(expr, string(s))]
end
idx_containing(nd, expr::Symbol) = idx_containing(nd, string(expr))

function construct_mass_matrix(mmv_array, gs)
    if all([mm == I for mm in mmv_array])
        mass_matrix = I
    else
        mass_matrix = sparse(1.0I, gs.dim_v, gs.dim_v)
        for (i, mm) in enumerate(mmv_array)
            ind = gs.v_idx[i]
            if ndims(mm) == 0
                copyto!(@view(mass_matrix[ind, ind]), mm * I)
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
        mass_matrix = sparse(1.0I, dim_nd, dim_nd)
        for (i, mm) in enumerate(mmv_array)
            ind = gs.v_idx[i]
            if ndims(mm) == 0
                copyto!(@view(mass_matrix[ind, ind]), mm * I)
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
                copyto!(@view(mass_matrix[ind, ind]), mm * I)
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


"""
    allocation_report(nd, p; verbose=true)

Tests the given `network_dynamics` object with parameters `p` for allocations.
Returns the overall allocations of the call to the rhs as well as the individual
allocations for all node and edge functions.

Returns `(a_nd, [a_vert..], [a_edge..])` where `a_nd` are the total allocations of the
nd call, `a_vert` are the allocations of the unique vertex functions and `a_edge` are the
allocations of the unique edge functions.
"""
function allocation_report(nd, p; verbose=true)
    ndf = nd.f;

    u  = rand(length(nd.syms));
    du = zeros(length(u));
    nd(du, u, p, 0.0) # warmup

    allocations_nd = @allocated nd(du, u, p, 0.0)
    verbose && @info "Execution of Network" allocations_nd

    # get minimum dimension of edges and vertices to create "fake" data
    gd = nd(u,p,0.0,GetGD)
    min_e_dim = mapreduce(min, get_dst_edges(gd, vertices(ndf.graph))) do edgelist
        minimum(length.(edgelist))
    end
    min_v_dim = mapreduce(v->v.dim, min, ndf.unique_vertices!)

    allocations_unique_verts = Int[]
    for (ui, vert) in enumerate(ndf.unique_vertices!)
        vert isa ODEVertex || error("Allocation report not supported for $(typeof(vert))")
        # get p for specific indices
        pv = NetworkDynamics.p_v_idx(p, first(ndf.unique_v_indices[ui]))
        u  = rand(vert.dim)
        du = zeros(vert.dim)
        edges = [rand(min_e_dim) for i in 1:5]

        vert.f(du, u, edges, pv, 0.0) # warmup
        allocations = @allocated vert.f(du, u, edges, pv, 0.0)
        verbose && @info "Execution of vertex model for $(ndf.unique_v_indices[ui])" allocations
        push!(allocations_unique_verts, allocations)
    end

    allocations_unique_edges = Int[]
    for (ui, edg) in enumerate(ndf.unique_edges!)
        # get p for specific indices
        pe = NetworkDynamics.p_e_idx(p, first(ndf.unique_e_indices[ui]))
        u  = rand(edg.dim)
        src = rand(min_v_dim)
        dst = rand(min_v_dim)

        if edg isa ODEEdge
            du = zeros(edg.dim)
            edg.f(du, u, src, dst, pe, 0.0) # warmup
            allocations = @allocated edg.f(du, u, src, dst, pe, 0.0)
            verbose && @info "Execution of edge model for $(ndf.unique_e_indices[ui])" allocations
        elseif edg isa StaticEdge
            edg.f(u, src, dst, pe, 0.0) # warmup
            allocations = @allocated edg.f(u, src, dst, pe, 0.0)
            verbose && @info "Execution of edge model for $(ndf.unique_e_indices[ui])" allocations
        else
            error("Allocation report not supportet for $(typeof(edg)).")
        end
        push!(allocations_unique_edges, allocations)
    end
    (allocations_nd, allocations_unique_verts, allocations_unique_edges)
end

allocations(nd, p) = allocation_report(nd, p; verbose=false)[1]
