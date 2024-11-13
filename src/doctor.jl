abstract type Access end
struct Read <: Access
    idx::Int
end
struct Write <: Access
    idx::Int
end
struct DeriveSimilar <: Access
    dims::Dims{}
end

mutable struct AccessTracker{T} <: AbstractVector{T}
    actions::Vector{Access}
    data::Vector{T}
end
function AccessTracker(data)
    AccessTracker(Access[], data)
end

Base.IndexStyle(::Type{<:AccessTracker}) = IndexLinear()
Base.size(a::AccessTracker) = size(a.data)
function Base.getindex(a::AccessTracker, i::Int)
    push!(a.actions, Read(i))
    if i in eachindex(a.data)
        return a.data[i]
    else
        return NaN
    end
end
function Base.setindex!(a::AccessTracker, v, i::Int)
    push!(a.actions, Write(i))
    if i in eachindex(a.data)
        return a.data[i] = v
    else
        v
    end
end
function Base.similar(a::AccessTracker, ::Type{T}, dims::Dims) where {T}
    push!(a.actions, DeriveSimilar(dims))
    AccessTracker(similar(a.data, T, dims))
end

Base.iterate(a::AccessTracker, i::Int=1) = (a[i], i+1)

# see julia docs on interfaces, used to track similar calls
Base.BroadcastStyle(::Type{<:AccessTracker}) = Broadcast.ArrayStyle{AccessTracker}()
function Base.similar(bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{AccessTracker}}, ::Type{ElType}) where ElType
    a = find_tracker(bc)
    push!(a.actions, DeriveSimilar(length.(axes(bc))))
    AccessTracker(similar(Array{ElType}, axes(bc)))
end
find_tracker(bc::Base.Broadcast.Broadcasted) = find_tracker(bc.args)
find_tracker(args::Tuple) = find_tracker(find_tracker(args[1]), Base.tail(args))
find_tracker(x) = x
find_tracker(::Tuple{}) = nothing
find_tracker(a::AccessTracker, rest) = a
find_tracker(::Any, rest) = find_tracker(rest)

function Base.show(io::IO, m::MIME"text/plain", a::AccessTracker)
    println(io, "AccessTracker with $(length(reads(a))) reads and $(length(writes(a))) writes:")
    show(io, m, a.data)
end
Base.show(io::IO, a::AccessTracker) = print(io, "AccessTracker($(a.data))")

reads(a::AccessTracker) = Int[action.idx for action in a.actions if action isa Read]
writes(a::AccessTracker) = Int[action.idx for action in a.actions if action isa Write]
has_reads(a::AccessTracker) = any(a -> isa(a, Read), a.actions)
has_writes(a::AccessTracker) = any(a -> isa(a, Write), a.actions)
has_similars(a::AccessTracker) = any(a -> isa(a, DeriveSimilar), a.actions)
function has_uninit_reads(a::AccessTracker)
    writtenidx = Int[]
    for action in a.actions
        if action isa Write
            push!(writtenidx, action.idx)
        end
        if action isa Read && !(action.idx in writtenidx)
            return true
        end
    end
    return false
end
oob_reads(a::AccessTracker) = filter(i -> i ∉ eachindex(a.data), reads(a))
oob_writes(a::AccessTracker) = filter(i -> i ∉ eachindex(a.data), writes(a))
has_oob(a::AccessTracker) = !isempty(oob_reads(a)) || !isempty(oob_writes(a))


# chk_component(::VertexFunction) = @warn "Check on VertexFunction not implemented"
# chk_component(::EdgeFunction) = @warn "Check on EdgeFunction not implemented"
function chk_component(c::ComponentFunction)
    @nospecialize
    du = AccessTracker(rand(dim(c)))
    u = AccessTracker(rand(dim(c)))
    p = AccessTracker(rand(pdim(c)))
    ins = if hasindim(c)
        Tuple(AccessTracker(rand(l)) for l in values(indim(c)))
    else
        # we don't know the size of the input but this might be reasonable guess
        indim_guess = max(outdim_normalized(c)...)
        Tuple(AccessTracker(rand(indim_guess)) for _ in outdim_normalized(c))
    end
    outs = Tuple(AccessTracker(rand(l)) for l in values(outdim(c)))
    t = NaN

    try
        compfg(c)(outs..., du, u, ins..., p, t)
    catch e
        if e isa MethodError
            @warn "Encountered MethodError. All arguments are AbstractArrays, make sure to allways index into them: $e"
        elseif e isa BoundsError
            @warn "Call of component functions lead to out of bounds access! Maybe you're unpacking some function input?"
        elseif e isa DimensionMismatch
            # ignore, its probably because we don't know the sizes of esum, vsrc and vdst
            if hasindim(c)
                @warn "Call of component function lead to dimension mismatch: $e."
            end
        else
            @warn "Error while calling component function: $e"
        end
        return nothing
    end

    # check for out of bound read access
    has_oob(du) && @warn "There is out of bound acces to du: reads $(oob_reads(du)) and writes $(oob_writes(du))! Check dim/sym!"
    has_oob(u) && @warn "There is out of bound acces to u: reads $(oob_reads(u)) and writes $(oob_writes(u))! Check dim/sym!"
    has_oob(p) && @warn "There is out of bound acces to p: reads $(oob_reads(p)) and writes $(oob_writes(p))! Check pdim/psim!"
    for (j, o) in enumerate(outs)
        has_oob(o) && @warn "There is out of bound acces to output#$j: reads $(oob_reads(o)) and writes $(oob_writes(o))!"
    end
    if hasindim(c)
        for (j, i) in enumerate(ins)
            has_oob(i) &&  @warn "There is out of bound acces to input#$j: reads $(oob_reads(i)) and writes $(oob_writes(i))!"
        end
    end

    has_uninit_reads(du) && @warn "There is uninitialized read access to du: $(reads(du))!"
    has_writes(u) && @warn "There is write access to u: $(writes(u))!"
    written = unique!(sort!(writes(du)))
    written != 1:dim(c) && @warn "Not all state idx 1:$(dim(c)) are set, only $(written)!"

    for (j, o) in enumerate(outs)
        has_uninit_reads(o) && @warn "There is uninitialized read access to output#$j: $(reads(o))!"
        written = unique!(sort!(writes(o)))
        written != 1:length(o) && @warn "Not all idx of output#$j 1:$(length(o)) were set, only $(written)!"
    end

    for i in ins
        has_writes(i) && @warn "There is write access to input: $(writes(i))!"
    end

    has_writes(p) && @warn "There is write access to p: $(writes(p))!"

    similars = String[]
    has_similars(du) && push!(similars, "du")
    has_similars(u) && push!(similars, "u")
    has_similars(p) && push!(similars, "p")
    if any(has_similars, ins)
        push!(similars, "inputs")
    end
    if any(has_similars, outs)
        push!(similars, "outputs")
    end
    if !isempty(similars)
        @warn "Component function allocates similar arrays to $(join(similars, ", "))!"
    end
    nothing
end

_ninout(::EdgeFunction) = 2
_ninout(::VertexFunction) = 1
