"""
    DependencyAwareObsf(f, columns, needs_input)

Observed function which can compute a subset of its observables.

The wrapped `f(out, u, ins..., p, t, mask)` evaluates its internal assignments, and `mask`
switches each of them on or off; `mask = nothing` computes everything. `columns[i]` marks the
assignments observable `i` needs, `needs_input[i]` says whether any of them reads an input.

Called like a plain obsf it computes everything. The keyword `required` restricts it to a set of
observables, the keyword `mask` takes a precomputed [`assignment_mask`](@ref) instead. Everything
that is not needed for the requested observables comes out as `NaN`.
"""
struct DependencyAwareObsf{F} <: Function
    f::F
    columns::Vector{BitVector}
    needs_input::BitVector
    function DependencyAwareObsf(f::F, columns, needs_input) where {F}
        isempty(columns) && throw(ArgumentError("DependencyAwareObsf needs at least one observable."))
        length(columns) == length(needs_input) || throw(ArgumentError("`columns` and `needs_input` differ in length."))
        allequal(length, columns) || throw(ArgumentError("All columns must cover the same assignments."))
        new{F}(f, columns, needs_input)
    end
end

function (o::DependencyAwareObsf)(out, args::Vararg{Any,N}; required=nothing, mask=nothing) where {N}
    if !isnothing(required)
        isnothing(mask) || throw(ArgumentError("Pass either `required` or `mask`, not both."))
        mask = assignment_mask(o, required)
    end
    o.f(out, args..., mask)
end

function Base.:(==)(a::DependencyAwareObsf, b::DependencyAwareObsf)
    a.f == b.f && a.columns == b.columns && a.needs_input == b.needs_input
end
function Base.hash(o::DependencyAwareObsf, h::UInt)
    hash(o.needs_input, hash(o.columns, hash(o.f, h)))
end

"""
    requires_input(obsf, i)
    requires_input(obsf, idxs)

Whether observable `i` (or any of `idxs`) of `obsf` depends on the component inputs. `idxs`
works like an index into the observables, so a Bool mask is fine as well. Plain obsf are assumed
to.
"""
requires_input(o::DependencyAwareObsf, i::Integer) = o.needs_input[i]
requires_input(o::DependencyAwareObsf, idxs) = any(view(o.needs_input, idxs))
requires_input(_, _) = true

"""
    assignment_mask(obsf, idxs)

Assignment mask which computes the observables `idxs`, as a fresh `BitVector`. `idxs` works
like an index into the observables, so a Bool mask is fine as well. Pass the result as
`obsf(out, u, ins..., p, t; mask)`. Returns `nothing` for a plain obsf, which can only be
evaluated fully.
"""
function assignment_mask(o::DependencyAwareObsf, idxs)
    mask = falses(length(o.columns[1]))
    for col in view(o.columns, idxs)
        _or!(mask, col)
    end
    mask
end
assignment_mask(_, _) = nothing

# word-level OR, much faster than broadcasting over the bits
function _or!(a::BitVector, b::BitVector)
    ac, bc = a.chunks, b.chunks
    @inbounds for k in eachindex(ac, bc)
        ac[k] |= bc[k]
    end
    a
end
