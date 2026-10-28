# The mode of both caches, `:off`, `:on` or `:check`, see `set_mtk_model_cache!`.
const MODEL_CACHE_MODE = Ref(:off)

# The properties which make up the structure of a system. Two systems are structurally equal
# when all of them agree, recursing through `systems`. Names, guesses, descriptions and most
# defaults are not on the list, those are read from each instance.
const STRUCTURAL_PROPERTIES = (
    eqs = ModelingToolkitBase.get_eqs,
    unknowns = ModelingToolkitBase.get_unknowns,
    ps = ModelingToolkitBase.get_ps,
    iv = ModelingToolkitBase.get_iv,
    observed = ModelingToolkitBase.get_observed,
    bindings = ModelingToolkitBase.get_bindings,
    connector_type = ModelingToolkitBase.get_connector_type,
    ignored_connections = ModelingToolkitBase.get_ignored_connections,
    initf = s -> SymbolicUtils.getmetadata(s, NetworkDynamics.SystemInitFormulas, nothing),
    guessf = s -> SymbolicUtils.getmetadata(s, NetworkDynamics.SystemGuessFormulas, nothing),
    postprocessing = s -> SymbolicUtils.getmetadata(s, ComponentPostprocessing, nothing),
    systems = ModelingToolkitBase.get_systems,
)

# Symbols hash without their metadata, so the key collects the variable metadata which the
# build reads on top.
const STRUCTURAL_VAR_METADATA = (
    VariableInitFormula,
    VariableInitFormulaWeak,
    VariableInitFormulaOptional,
    VariableGuessFormula,
    ModelingToolkitBase.VariableConnectType,
    ModelingToolkitBase.VariableInput,
    ModelingToolkitBase.VariableOutput
)

# The key is a hash of the structure of `sys` and of the other arguments of
# `generate_io_function`. With 64 bits, collisions between a few hundred entries are too
# unlikely to be worth comparing the systems themselves.
function _structural_key(sys, inputss, outputss, kwargs)
    args = (
        map(x -> getname.(x), inputss),
        map(x -> getname.(x), outputss),
        NamedTuple(kwargs)
    )
    _structural_hash(sys, hash(args))
end

# Hashes every structural property of `sys`, then the metadata of its variables. The name of
# `sys` itself does not count, children are hashed together with their name.
function _structural_hash(sys, h::UInt)
    # build incremental hash of all structural properties (recurses on system)
    for (name, prop) in pairs(STRUCTURAL_PROPERTIES)
        h = _hash_property(name, prop(sys), h)
    end
    # build incremental hash of the variable metadata the build reads
    for vars in (ModelingToolkitBase.get_unknowns(sys), ModelingToolkitBase.get_ps(sys))
        for v in vars
            var_metadata = map(K -> SymbolicUtils.getmetadata(v, K, nothing), STRUCTURAL_VAR_METADATA)
            h = hash(var_metadata, h)
        end
    end
    # parameters of a domain connector which have a default get bound across its connections
    if ModelingToolkitBase.is_domain_connector(sys)
        ics = ModelingToolkitBase.get_initial_conditions(sys)
        h = hash(map(p -> !isnothing(get(ics, p, nothing)), ModelingToolkitBase.get_ps(sys)), h)
    end
    h
end
function _hash_property(name, x, h)
    if name === :systems
        foldl((h, child) -> _hash_child(child, h), x; init=h)
    elseif name === :eqs || name === :observed
        foldl((h, eq) -> hash(_comparable_eq(eq), h), x; init=h)
    else
        hash(x, h)
    end
end
_hash_child(child, h) = _structural_hash(child, hash(nameof(child), h))

# The systems in a connection equation hash by identity. Their names are enough to describe
# the connection, since the systems themselves are hashed as children.
function _comparable_eq(eq::Equation)
    lhs, rhs = unwrap_const(eq.lhs), unwrap_const(eq.rhs)
    lhs isa ModelingToolkitBase.Connection && rhs isa ModelingToolkitBase.Connection || return eq
    conns = map(s -> s isa ModelingToolkitBase.AbstractSystem ? nameof(s) : s, rhs.systems)
    (lhs.systems, conns)
end

# `value` stays `nothing` until `owner`, the task which created the entry, has built it.
# The owner is dropped then, so that the entry doesn't keep the task alive.
mutable struct CacheEntry
    const ready::Base.Event
    owner::Union{Task, Nothing}
    value::Any
end

# A session-wide LRU cache which several tasks can fill at once. Beyond `maxsize` entries, the
# least recently used ones are dropped. The lock is never held while a value is built.
mutable struct ThreadsafeCache{K}
    const name::String
    const lock::ReentrantLock
    # insertion order doubles as recency order, the front is dropped first
    const entries::OrderedDict{K, CacheEntry}
    maxsize::Int
    hits::Int
    misses::Int
    waitings::Int
end
function ThreadsafeCache{K}(name; maxsize=500) where {K}
    ThreadsafeCache{K}(name, ReentrantLock(), OrderedDict{K, CacheEntry}(), maxsize, 0, 0, 0)
end

# the structure of a system to the output of `generate_io_function`
const COMPILE_CACHE = ThreadsafeCache{UInt}("MTK Model Cache")
# a `@cached_component` call to the system the user's function built
const COMPONENT_CACHE = ThreadsafeCache{Any}("Cached Component")
const ALL_CACHES = (COMPILE_CACHE, COMPONENT_CACHE)

# Returns the cached value for `key`, calling `f` to fill the entry on a miss. Other tasks
# asking for the same key wait for that build. If it fails, they call `f` themselves.
function Base.get!(f, cache::ThreadsafeCache, key)
    entry, isowner = @lock cache.lock begin
        entry = get(cache.entries, key, nothing)
        if isnothing(entry)
            cache.misses += 1
            entry = cache.entries[key] = CacheEntry(Base.Event(), current_task(), nothing)
            while length(cache.entries) > cache.maxsize
                delete!(cache.entries, first(keys(cache.entries)))
            end
            entry, true
        else
            cache.hits += 1
            if isnothing(entry.value)
                cache.waitings += 1
            end
            # move to the back
            delete!(cache.entries, key)
            cache.entries[key] = entry
            entry, false
        end
    end

    if isowner
        try
            value = f()
            @lock cache.lock begin
                entry.value = value
                entry.owner = nothing
            end
        catch
            @lock cache.lock begin
                get(cache.entries, key, nothing) === entry && delete!(cache.entries, key)
            end
            rethrow()
        finally
            notify(entry.ready)
        end
    else
        # the task building the entry would wait for itself forever
        if entry.owner === current_task() && isnothing(entry.value)
            error("A cache entry was requested again while it was being built. For a cached component \
                which calls itself, a `@structural` expression probably hides the choice which ends the recursion.")
        end
        wait(entry.ready)
    end
    isnothing(entry.value) ? f() : entry.value
end

function Base.empty!(cache::ThreadsafeCache)
    @lock cache.lock empty!(cache.entries)
    cache
end

function reset_counts!(cache::ThreadsafeCache)
    @lock cache.lock begin
        cache.hits = cache.misses = cache.waitings = 0
    end
    cache
end

function print_stats(cache::ThreadsafeCache)
    @lock cache.lock begin
        (; hits, misses, waitings) = cache
        println(cache.name, " Stats:")
        println("  Entries:  ", length(cache.entries), " (max ", cache.maxsize, ")")
        println("  Misses:   ", misses)
        println("  Hits:     ", hits, " (", waitings, " after waiting)")
        println("  Hit Rate: ", round(hits / max(hits + misses, 1) * 100, digits=2), "%")
    end
end

# `generate_io_function`, reused for systems with the same structure. Callers read all
# per-instance data (defaults, guesses, metadata) from `sys`.
function generate_io_function_cached(sys, inputss, outputss; verbose=false, mtkcompile=false, kwargs...)
    build() = generate_io_function(sys, inputss, outputss; verbose, mtkcompile, kwargs...)
    # Verbose builds run fresh so that all their messages show up. The key only knows what ND's
    # own simplification reads, so `mtkcompile` builds are not cached.
    mode = MODEL_CACHE_MODE[]
    if mode === :off || verbose || mtkcompile
        return build()
    end
    hit = true
    key = _structural_key(sys, inputss, outputss, kwargs)
    gen = get!(COMPILE_CACHE, key) do
        hit = false
        build()
    end
    if hit && mode === :check
        rebuilt = build()
        diff = [k for k in keys(rebuilt) if !_samebuild(gen[k], rebuilt[k])]
        isempty(diff) || error("MTK model cache returned a different build for system $(getname(sys)), \
            differing in $(join(diff, ", ")). Please report this as a bug. \
            `NetworkDynamics.set_mtk_model_cache!(:off)` switches the cache off.")
    end
    # the containers which end up in the component are not shared with the cache
    fresh = map(copy, (; gen.obsrules, gen.aliasmap, gen.observed, gen.equations,
                         gen.full_equations, gen.outputeqs, gen.full_outputeqs, gen.mass_matrix))
    (; gen..., fresh...)
end

# Compares a cache hit with a fresh build. Generated functions count as equal when their types
# agree, since the type of a RuntimeGeneratedFunction carries the hash of its body.
function _samebuild(a, b)
    if a isa AbstractVector && b isa AbstractVector
        length(a) == length(b) && all(i -> _samebuild(a[i], b[i]), eachindex(a, b))
    elseif a isa Union{ResolutionRule, InitFormula, GuessFormula, NetworkDynamics.DependencyAwareObsf}
        typeof(a) == typeof(b) && all(fn -> _samebuild(getfield(a, fn), getfield(b, fn)), fieldnames(typeof(a)))
    elseif a isa Function
        typeof(a) == typeof(b)
    else
        isequal(a, b)
    end
end

"""
    NetworkDynamics.set_mtk_model_cache!(mode; maxsize=nothing)

Set the mode of the MTK model cache:
- `:off` (default): every component is built from scratch. Switching off also empties the cache.
- `:on`: components with the same structure as an earlier one reuse its generated functions.
- `:check`: like `:on`, but every cache hit is also built from scratch and compared with the cached
  result. Errors if they differ. This is slow and meant for debugging.

`true` and `false` stand for `:on` and `:off`. `maxsize` sets how many distinct models the cache
keeps (500 by default).

The mode also applies to constructors defined with [`@cached_component`](@ref), which keep their
own cache of the same size.

See also [`NetworkDynamics.with_mtk_model_cache`](@ref).
"""
function NetworkDynamics.set_mtk_model_cache!(mode; maxsize=nothing)
    mode = cache_mode(mode)
    for cache in ALL_CACHES
        @lock cache.lock begin
            if !isnothing(maxsize)
                cache.maxsize = maxsize
            end
            mode === :off && empty!(cache.entries)
        end
    end
    MODEL_CACHE_MODE[] = mode
end
cache_mode(use::Bool) = use ? :on : :off
function cache_mode(mode::Symbol)
    mode in (:off, :on, :check) || throw(ArgumentError("Cache mode must be :off, :on or :check, got :$mode."))
    mode
end

"""
    NetworkDynamics.wipe_mtk_model_cache!()

Empty the MTK model cache and the cache of [`@cached_component`](@ref) constructors. Needed when
the code generating models changes within a session, for example when editing NetworkDynamics or
a cached component with Revise.
"""
function NetworkDynamics.wipe_mtk_model_cache!()
    foreach(empty!, ALL_CACHES)
    nothing
end

"""
    NetworkDynamics.with_mtk_model_cache(f, mode=:on)

Run `f` with the MTK model cache in the given mode (see [`NetworkDynamics.set_mtk_model_cache!`](@ref))
and restore the previous mode afterwards. Also resets the counters of
[`NetworkDynamics.mtk_cache_stats`](@ref).

```julia
NetworkDynamics.with_mtk_model_cache(false) do
    # code which builds mtk models without the cache
end
```

The cache is off by default. When it is on and a [`VertexModel`](@ref) or [`EdgeModel`](@ref) is built from
a `System` with the same equations, variables, bindings, init and guess formulas, subsystems
and inputs/outputs as an earlier one, the generated functions and equations are reused. Defaults,
guesses and other metadata are always read from the new system. Warnings of the build only show
up once, when the cache entry is filled. Builds with `verbose=true` or `mtkcompile=true` bypass
the cache. The cache is thread safe, so components can be built concurrently.

See also [`NetworkDynamics.set_mtk_model_cache!`](@ref), [`NetworkDynamics.wipe_mtk_model_cache!`](@ref).
"""
function NetworkDynamics.with_mtk_model_cache(f, mode=:on)
    oldmode = MODEL_CACHE_MODE[]
    foreach(reset_counts!, ALL_CACHES)
    NetworkDynamics.set_mtk_model_cache!(mode)
    try
        f()
    finally
        NetworkDynamics.set_mtk_model_cache!(oldmode)
    end
end

"""
    NetworkDynamics.mtk_cache_stats()

Print statistics about the MTK model cache and the cache of [`@cached_component`](@ref)
constructors since the session started or since the last
[`NetworkDynamics.with_mtk_model_cache`](@ref) block began.
"""
function NetworkDynamics.mtk_cache_stats()
    foreach(print_stats, ALL_CACHES)
end

function NetworkDynamics.mtk_model_cache_enabled()
    MODEL_CACHE_MODE[] !== :off
end

# the default of the `cache` kwarg of `@cached_component` constructors
NetworkDynamics.mtk_model_cache_mode() = MODEL_CACHE_MODE[]
