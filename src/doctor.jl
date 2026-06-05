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


# chk_component(::VertexModel) = @warn "Check on VertexModel not implemented"
# chk_component(::EdgeModel) = @warn "Check on EdgeModel not implemented"
function chk_component(c::ComponentModel)
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
    ext = AccessTracker(rand(extdim(c)))
    if has_external_input(c)
        ins = (ins..., ext)
    end
    t = NaN

    try
        compfg(c)(outs, du, u, ins, p, t)
    catch e
        if e isa MethodError
            @warn "Encountered MethodError. All arguments are AbstractArrays, make sure to allways index into them: $e"
        elseif e isa BoundsError
            @warn "Call of component models lead to out of bounds access! Maybe you're unpacking some function input?"
        elseif e isa DimensionMismatch
            # ignore, its probably because we don't know the sizes of esum, vsrc and vdst
            if hasindim(c)
                @warn "Call of component model lead to dimension mismatch: $e."
            end
        else
            @warn "Error while calling component model: $e"
        end
        return nothing
    end

    # check for out of bound read access
    has_oob(du) && @warn "There is out of bound acces to du: reads $(oob_reads(du)) and writes $(oob_writes(du))! Check dim/sym!"
    has_oob(u) && @warn "There is out of bound acces to u: reads $(oob_reads(u)) and writes $(oob_writes(u))! Check dim/sym!"
    has_oob(p) && @warn "There is out of bound acces to p: reads $(oob_reads(p)) and writes $(oob_writes(p))! Check pdim/psim!"
    has_oob(ext) && @warn "There is out of bound acces to external input: reads $(oob_reads(ext)) and writes $(oob_writes(ext))!"
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
    has_writes(ext) && @warn "There is write access to external inputs: $(writes(ext))!"

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
        @warn "Component model allocates similar arrays to $(join(similars, ", "))!"
    end

    # smoketest for observed function
    if !isnothing(c.obsf)
        out = zeros(length(obssym(c)))
        _ins = map(t -> t.data, ins)
        _u = u.data
        _p = p.data
        try
            c.obsf(out, _u, _ins..., _p, t)
        catch e
            @warn "Component Check: Error while calling observable function!" exception=(e, catch_backtrace())
        end
    end
    nothing
end

_ninout(::EdgeModel) = 2
_ninout(::VertexModel) = 1

####
#### scope consistency checks
####
# basename used for grouping scoped parameters: the trailing part of the symbol
# name after the last namespace separator `₊`, i.e. matching `r"NAME$"`.
_scope_basename(sym) = Symbol(last(split(string(sym), '₊')))

"""
    chk_global_parameters(nw::Network; verbose=true)
    chk_global_parameters(s::NWState; verbose=true)
    chk_global_parameters(p::NWParameter; verbose=true)

Check the consistency of *scoped* parameters (see [`VariableScope`](@ref)).

Parameters are grouped by the trailing part of their symbol name (matching
`r"NAME\$"`):

- parameters with scope `:global` must hold the same value across the **whole
  network**,
- parameters with scope `:device` must hold the same value **within a single**
  `VertexModel`/`EdgeModel`,
- parameters with scope `:local` are ignored.

For a `Network` the metadata **defaults** are compared. For an `NWState`/
`NWParameter` the **current values** in the flat parameter array are compared.

Returns `true` if all scoped parameters are consistent, `false` otherwise. When
`verbose=true` (the default) a warning describing each inconsistency is emitted.

See also: [`VariableScope`](@ref), [`get_scope`](@ref).
"""
function chk_global_parameters end

chk_global_parameters(nw::Network; kw...) = _chk_global_parameters(nw, nothing; kw...)
chk_global_parameters(p::NWParameter; kw...) = _chk_global_parameters(p.nw, p; kw...)
chk_global_parameters(s::NWState; kw...) = _chk_global_parameters(s.nw, s.p; kw...)

function _chk_global_parameters(nw::Network, valueprovider; verbose=true)
    # value getter for a parameter `sym` of component `c` at index `ci` (idxtype VIndex/EIndex)
    getval = if isnothing(valueprovider)
        (idxtype, ci, c, sym) -> has_default(c, sym) ? get_default(c, sym) : nothing
    else
        (idxtype, ci, c, sym) -> valueprovider[idxtype(ci, sym)]
    end

    ok = true
    # global groups span the whole network: basename => Vector{Pair{value, location}}
    globalgroups = OrderedDict{Symbol, Vector{Pair{Any,String}}}()

    function _handle(c, ci, idxtype, prefix)
        # scope is only meaningful (and enforced) for parameters; warn if it was
        # accidentally attached to a state/observable/input/output symbol
        if verbose
            pset = Set(psym(c))
            for (s, smd) in symmetadata(c)
                if haskey(smd, :scope) && s ∉ pset
                    @warn "Symbol :$s in $prefix$ci carries `:scope` metadata but is not a parameter. \
                           Scope is only enforced for parameters and is ignored here."
                end
            end
        end
        # device groups are local to this single component
        devicegroups = OrderedDict{Symbol, Vector{Pair{Any,String}}}()
        for s in psym(c)
            has_scope(c, s) || continue
            scope = get_scope(c, s)
            scope == :local && continue
            val = getval(idxtype, ci, c, s)
            isnothing(val) && continue
            base = _scope_basename(s)
            loc = "$prefix$ci :$s"
            if scope == :global
                push!(get!(globalgroups, base, Pair{Any,String}[]), val => loc)
            elseif scope == :device
                push!(get!(devicegroups, base, Pair{Any,String}[]), val => loc)
            else
                throw(ArgumentError("Unknown scope $(repr(scope)) for $loc. Must be :local, :device or :global."))
            end
        end
        for (base, entries) in devicegroups
            ok &= _chk_consistent_scope("device parameter :$base ($prefix$ci)", entries, verbose)
        end
    end

    for (ci, c) in pairs(nw.im.vertexm)
        _handle(c, ci, VIndex, "VIndex ")
    end
    for (ci, c) in pairs(nw.im.edgem)
        _handle(c, ci, EIndex, "EIndex ")
    end
    for (base, entries) in globalgroups
        ok &= _chk_consistent_scope("global parameter :$base", entries, verbose)
    end
    ok
end

function _chk_consistent_scope(desc, entries, verbose)
    vals = unique(first.(entries))
    length(vals) <= 1 && return true
    if verbose
        msg = "Inconsistent $desc:\n" * join(("  $loc = $val" for (val, loc) in entries), "\n")
        @warn msg
    end
    false
end
