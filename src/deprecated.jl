# Deprecated non plural methods
@deprecate delete_initconstraint!(args...; kwargs...) delete_initconstraints!(args...; kwargs...)
@deprecate delete_initformula!(args...; kwargs...) delete_initformulas!(args...; kwargs...)

# Deprecated names with typos
const ContinousComponentCallback = ContinuousComponentCallback
const VectorContinousComponentCallback = VectorContinuousComponentCallback

# Deprecated symbolic indexing methods
vidxs(compfilter::Union{Int,AbstractVector{Int},NTuple{<:Any,Int},Colon}, varfilter) = _deprecated_idxs_gen(VIndex, compfilter, varfilter)
eidxs(compfilter::Union{Int,AbstractVector{Int},NTuple{<:Any,Int},Colon}, varfilter) = _deprecated_idxs_gen(EIndex, compfilter, varfilter)
vpidxs(compfilter::Union{Int,AbstractVector{Int},NTuple{<:Any,Int},Colon}, varfilter) = _deprecated_idxs_gen(VPIndex, compfilter, varfilter, VIndex)
epidxs(compfilter::Union{Int,AbstractVector{Int},NTuple{<:Any,Int},Colon}, varfilter) = _deprecated_idxs_gen(EPIndex, compfilter, varfilter, EIndex)

# Deprecated linear analysis methods with nw as first argument
function isfixpoint(nw::Network, s0::NWState; kwargs...)
    @warn "isfixpoint(nw, s0) is deprecated, use isfixpoint(s0) instead." maxlog=1
    extract_nw(s0) === nw || throw(ArgumentError("Network mismatch: s0 is not associated with the given nw"))
    isfixpoint(s0; kwargs...)
end

function jacobian_eigenvals(nw::Network, s0; kwargs...)
    @warn "jacobian_eigenvals(nw, s0) is deprecated, use jacobian_eigenvals(s0) instead." maxlog=1
    extract_nw(s0) === nw || throw(ArgumentError("Network mismatch: s0 is not associated with the given nw"))
    jacobian_eigenvals(s0; kwargs...)
end

function is_linear_stable(nw::Network, s0; kwargs...)
    @warn "is_linear_stable(nw, s0) is deprecated, use jacobian_eigenvals(s0) instead." maxlog=1
    extract_nw(s0) === nw || throw(ArgumentError("Network mismatch: s0 is not associated with the given nw"))
    is_linear_stable(s0; kwargs...)
end

function _deprecated_idxs_gen(constructor, compiter, variter, type=constructor)
    @warn "*idxs(compfilter, varfilter) methods are deprecated. Use *idxs(nw, compfilter, varfilter) instead."
    if compiter isa Colon
        throw(ArgumentError("compfilter cannot be `:`, use *idxs(nw_like_thing,  ...) for that!"))
    end
    if variter isa Colon
        throw(ArgumentError("varfilter cannot be `:`, use *idxs(nw_like_thing,  ...) for that!"))
    end
    if variter isa Symbol
        variter = (variter,)
    end
    indices = type[]
    for ci in compiter
        for si in variter
            push!(indices, constructor(ci, si))
        end
    end
    indices
end
