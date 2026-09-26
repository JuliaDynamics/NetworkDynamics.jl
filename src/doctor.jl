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


"""
    chk_component(c::ComponentModel; ad=true)

Call the component functions with test data and warn about wrong array access, allocations
and failed calls of the observable function. With `ad=true` f and g are also called with
ForwardDiff Duals, as implicit solvers do for the Jacobian. That catches models which can't
take Duals or which only allocate for them.

Runs on every component construction, there without the Dual pass, see
`NetworkDynamics.CHECK_COMPONENT`.

If f only allocates with Duals, it is most likely type unstable. `@code_warntype` marks the
culprits in red, e.g. for a vertex with 2 states and a 1-dim input:

```julia
using ForwardDiff: Dual
D = Dual{Nothing,Float64,1}
du, u, in = zeros(D, 2), D.(rand(2)), D.(rand(1))
p = rand(pdim(c))
@code_warntype NetworkDynamics.compf(c)(du, u, in, p, 0.0)
```

For MTK models f is a `RuntimeGeneratedFunction`, and the call above only shows a wrapper.
Its body shows up with

```julia
@code_warntype NetworkDynamics.RuntimeGeneratedFunctions.generated_callfunc(
    NetworkDynamics.compf(c), du, u, in, p, 0.0)
```
"""
function chk_component(c::ComponentModel; ad=true)
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

    D = ForwardDiff.Dual{ForwardDiff.Tag{typeof(chk_component),Float64},Float64,1}
    report = _allocation_report(c, length(du), length.(ins), length.(outs), length(p);
                                dualT = ad ? D : nothing)
    if !isempty(report)
        @warn "Component :$(c.name)\n" * join(("- $name $problem" for (name, problem, _) in report), "\n") *
              "\n" * join(unique(last.(report)), "\n")
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

"""
    chk_network(nw::Network; io=stderr)

Check the network rhs `nw(du, u, p, t)` for allocations and print the result to `io`. It is
called with Float64 and with the ForwardDiff Duals an implicit solver uses, Dual states for the
Jacobian and a Dual time for the time derivative. If any of these calls allocates, each batch
of components is checked once like in [`chk_component`](@ref) and the allocating ones are listed.

The check runs on a copy with sequential execution and aggregation, since the parallel ones
allocate by design.
"""
function chk_network(nw::Network; io=stderr)
    # threaded execution and aggregation allocate on their own, so check a sequential copy
    if !(executionstyle(nw) isa SequentialExecution && nw.layer.aggregator isa SequentialAggregator)
        nw = Network(nw; execution=SequentialExecution{usebuffer(executionstyle(nw))}(),
                     aggregator=SequentialAggregator(aggfun(nw.layer.aggregator)),
                     copy_components=false, sparse=false)
    end
    N = ad_chunksize(nw.im)
    D = ForwardDiff.Dual{ForwardDiff.Tag{typeof(chk_network),Float64},Float64,N}
    Dt = ForwardDiff.Dual{ForwardDiff.Tag{typeof(chk_network),Float64},Float64,1}
    n = dim(nw)
    u, p = rand(n), rand(pdim(nw))
    counts = try
        ("Float64"     => _allocations(nw, zeros(n), u, p, 0.0),
         "Dual states" => _allocations(nw, zeros(D, n), D.(u), p, 0.0),
         "Dual time"   => _allocations(nw, zeros(Dt, n), u, p, Dt(0.0)))
    catch e
        printstyled(io, "✗ "; color=:red)
        println(io, "Calling the network with random states failed: ", _errorline(e))
        return nothing
    end
    if all(iszero ∘ last, counts)
        printstyled(io, "✓ "; color=:green)
        println(io, "The network rhs doesn't allocate, neither with Float64 nor with Duals.")
        return nothing
    end

    printstyled(io, "The network rhs allocates on every call\n"; bold=true)
    for (label, n) in counts
        print(io, "  ", rpad(label, 12))
        printstyled(io, lpad(n, 4), " allocations\n"; color = iszero(n) ? :green : :red)
    end
    allreports = Tuple{String,String,String}[]
    for (kind, batches, comps) in (("vertices", nw.vertexbatches, nw.im.vertexm),
                                   ("edges", nw.layer.edgebatches, nw.im.edgem))
        for batch in batches
            c = comps[first(batch.indices)]
            report = _allocation_report(c, _batch_dims(c, batch)...; dualT=D)
            isempty(report) && continue
            print(io, "  ")
            printstyled(io, ":", c.name; bold=true)
            printstyled(io, " ($kind $(_shortlist(batch.indices)))\n"; color=:light_black)
            _print_report(io, report)
            append!(allreports, report)
        end
    end
    if isempty(allreports)
        println(io, "  No component function allocates on its own, so the allocations come from the \
                     network code around them.")
    end
    _print_hints(io, allreports)
    nothing
end
# the dims one batch member sees in the coreloop
function _batch_dims(c, batch)
    indims = if c isa VertexModel
        (length(in_range(batch, 1)),)
    else
        (length(in_range(batch, 1, :src)), length(in_range(batch, 1, :dst)))
    end
    has_external_input(c) && (indims = (indims..., extdim(c)))
    (dim(c), indims, Tuple(values(outdim(c))), pdim(c))
end
function _shortlist(idxs)
    length(idxs) <= 10 ? join(idxs, ", ") : join(first(idxs, 10), ", ") * ", …"
end

"""
    _allocation_report(c, dim, indims, outdims, pdim; dualT)

Call f and g with Float64 and, unless `dualT` is `nothing`, with Duals of that type. Returns a
list of `(function name, problem, hint)`, empty if nothing is wrong. An error of the Float64
call is not reported here, the access tracker check covers that.
"""
function _allocation_report(c, dim, indims, outdims, pdim; dualT)
    float = _component_allocations(c, dim, indims, outdims, pdim, Float64)
    dual = isnothing(dualT) ? nothing : _component_allocations(c, dim, indims, outdims, pdim, dualT)
    report = Tuple{String,String,String}[]
    for name in (:f, :g)
        fb = float[name]
        db = isnothing(dual) ? nothing : dual[name]
        (isnothing(fb) || fb isa Exception) && continue
        if fb > 0
            push!(report, (string(name), "allocates $(_times(fb)) per call",
                           "For MTK models, look for array expressions like `array_literal` in \
                            `NetworkDynamics.pretty_f(c)`."))
        end
        if db isa Exception
            push!(report, (string(name), "fails with ForwardDiff Duals: $(_errorline(db))",
                           "Implicit solvers call f and g with Duals for the Jacobian, so buffers \
                            inside the model must take Duals too."))
        elseif db isa Integer && db > 0 && fb == 0
            push!(report, (string(name), "allocates $(_times(db)) per call, but only with ForwardDiff Duals",
                           "This is usually a type instability, see `?chk_component` for how to \
                            find it."))
        end
    end
    report
end

function _print_report(io, report)
    for (name, problem, _) in report
        print(io, "    ")
        printstyled(io, "✗ "; color=:red)
        printstyled(io, name; bold=true)
        println(io, " ", problem)
    end
end
function _print_hints(io, report)
    for hint in unique(last.(report))
        printstyled(io, "  ", hint, "\n"; color=:light_black)
    end
end
_errorline(e) = first(split(sprint(showerror, e), '\n'))
_times(n) = n == 1 ? "once" : "$n times"

"""
    _component_allocations(c, dim, indims, outdims, pdim, T)

Count the allocations of `f` and `g` per call, with states, inputs and outputs of eltype `T`.
The arguments are views into plain vectors, the same types the coreloop passes, so this matches
what the network will see. Each entry of the returned `(; f, g)` is the count, the exception if
the call errors, or `nothing` if there is no f.
"""
function _component_allocations(c, dim, indims, outdims, pdim, ::Type{T}) where {T}
    _view(n) = view(T.(rand(n)), 1:n)
    du, u = _view(dim), _view(dim)
    # like the coreloop, pass nothing instead of an empty parameter view; p and t stay Float64
    p = iszero(pdim) ? nothing : view(rand(pdim), 1:pdim)
    ins = map(_view, indims)
    outs = map(_view, outdims)
    t = 0.0
    f = isnothing(compf(c)) ? nothing : _try_allocations(apply_compf, compf(c), du, u, ins, p, t)
    gp = fftype(c) == PureStateMap() ? nothing : p
    g = _try_allocations(apply_compg, fftype(c), compg(c), outs, u, ins, gp, t)
    (; f, g)
end
function _try_allocations(args...)
    try
        _allocations(args...)
    catch e
        e
    end
end

# function barrier, so dispatch on the component types doesn't count as allocation
@noinline function _allocations(apply::A, args...) where {A}
    apply(args...)
    @allocations apply(args...)
end

_ninout(::EdgeModel) = 2
_ninout(::VertexModel) = 1
