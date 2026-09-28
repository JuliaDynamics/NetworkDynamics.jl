"""
    ParameterSnapshot{T,L}

Lazy view of the flat parameter vector at save index `i` of a [`ParameterLog`](@ref).
Single entries are looked up in the history, `copy` or `copyto!` rebuild the full vector.
"""
struct ParameterSnapshot{T,L} <: AbstractVector{T}
    log::L
    i::Int
end

"""
    ParameterLog{T}

History of the flat parameter vector over all `save_parameters!` calls. It stores the
initial vector once and afterwards only the values which actually changed, sorted by
parameter. This keeps saving cheap when callbacks flip a few parameters of a large network
over and over again.

It is used as the `u` field of the parameter timeseries `DiffEqArray`. Indexing it returns
a lazy [`ParameterSnapshot`](@ref) instead of a full vector.
"""
mutable struct ParameterLog{T} <: AbstractVector{ParameterSnapshot{T,ParameterLog{T}}}
    base::Vector{T}
    head::Vector{T}
    nsaves::Int
    # per parameter the (save index, new value) of each change, unassigned if it never changed
    tracks::Vector{Vector{Tuple{Int,T}}}
end

function ParameterLog(p::AbstractVector{T}) where {T}
    log = ParameterLog{T}(T[], T[], 0, Vector{Tuple{Int,T}}[])
    push!(log, p)
end

function Base.copy(log::ParameterLog)
    tracks = similar(log.tracks)
    for j in eachindex(tracks)
        isassigned(log.tracks, j) && (tracks[j] = copy(log.tracks[j]))
    end
    ParameterLog(copy(log.base), copy(log.head), log.nsaves, tracks)
end

Base.size(log::ParameterLog) = (log.nsaves,)
Base.IndexStyle(::Type{<:ParameterLog}) = IndexLinear()
function Base.getindex(log::ParameterLog, i::Int)
    @boundscheck checkbounds(log, i)
    ParameterSnapshot{eltype(log.base),typeof(log)}(log, i)
end

function Base.push!(log::ParameterLog{T}, p::AbstractVector) where {T}
    p isa Vector || (p = Array(p))
    log.nsaves += 1
    if log.nsaves == 1
        log.base = copy(p)
        log.head = copy(p)
        log.tracks = Vector{Vector{Tuple{Int,T}}}(undef, length(p))
        return log
    end
    length(p) == length(log.head) || throw(DimensionMismatch("Parameter vector changed size."))
    for j in eachindex(p)
        isequal(p[j], log.head[j]) && continue
        log.head[j] = p[j]
        isassigned(log.tracks, j) || (log.tracks[j] = Tuple{Int,T}[])
        push!(log.tracks[j], (log.nsaves, p[j]))
    end
    log
end

# value of a track at save `i`, or `default` if it didn't change up to then
function _track_value(track, i, default)
    k = searchsortedlast(track, (i,); by=first)
    iszero(k) ? default : track[k][2]
end

Base.size(s::ParameterSnapshot) = size(s.log.base)
Base.IndexStyle(::Type{<:ParameterSnapshot}) = IndexLinear()
function Base.getindex(s::ParameterSnapshot, j::Int)
    @boundscheck checkbounds(s, j)
    log = s.log
    isassigned(log.tracks, j) || return log.base[j]
    _track_value(log.tracks[j], s.i, log.base[j])
end

function Base.copyto!(dst::AbstractVector, s::ParameterSnapshot)
    log = s.log
    copyto!(dst, log.base)
    for j in eachindex(log.tracks)
        isassigned(log.tracks, j) || continue
        dst[j] = _track_value(log.tracks[j], s.i, log.base[j])
    end
    dst
end
Base.copy(s::ParameterSnapshot) = copyto!(similar(s.log.base), s)
Base.Vector(s::ParameterSnapshot) = copy(s)
