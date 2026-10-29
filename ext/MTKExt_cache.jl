# The mode of both caches, `:off`, `:on` or `:check`, see `set_mtk_cache!`.
const MTK_CACHE_MODE = Ref(:off)

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
# The owner is dropped then, so that the entry doesn't keep the task alive. The label, the
# timings and the hits are only there for `print_stats`. `keytime` adds up the time spent on
# computing the key over all calls which ended up at this entry.
mutable struct CacheEntry
    const ready::Base.Event
    const label::String
    owner::Union{Task, Nothing}
    value::Any
    buildtime::Float64
    compiletime::Float64
    keytime::Float64
    hits::Int
end

# A session-wide LRU cache which several tasks can fill at once. Beyond `maxsize` entries, the
# least recently used ones are dropped. The lock is never held while a value is built.
mutable struct ThreadsafeCache
    const name::String
    const lock::ReentrantLock
    # insertion order doubles as recency order, the front is dropped first
    const entries::OrderedDict{UInt, CacheEntry}
    maxsize::Int
    evictions::Int
end
function ThreadsafeCache(name; maxsize=500)
    ThreadsafeCache(name, ReentrantLock(), OrderedDict{UInt, CacheEntry}(), maxsize, 0)
end

# the structure of a system to the output of `generate_io_function`
const COMPILE_CACHE = ThreadsafeCache("Compile Cache")
# a `@cached_component` call to the system the user's function built
const COMPONENT_CACHE = ThreadsafeCache("Component Cache")
const ALL_CACHES = (COMPILE_CACHE, COMPONENT_CACHE)

# Returns the cached value for `key`, calling `f` to fill the entry on a miss. Other tasks
# asking for the same key wait for that build. If it fails, they call `f` themselves.
# For the statistics, `keytime` is the `@timed` computation of the key, and `label()` describes
# the entry. It is only called on a miss.
function Base.get!(f, cache::ThreadsafeCache, key::UInt; label=() -> "", keytime=nothing)
    keytime = isnothing(keytime) ? 0.0 : _runtime(keytime)
    entry, isowner = @lock cache.lock begin
        entry = get(cache.entries, key, nothing)
        if isnothing(entry)
            entry = CacheEntry(Base.Event(), label(), current_task(), nothing, NaN, NaN, keytime, 0)
            cache.entries[key] = entry
            while length(cache.entries) > cache.maxsize
                delete!(cache.entries, first(keys(cache.entries)))
                cache.evictions += 1
            end
            entry, true
        else
            entry.hits += 1
            entry.keytime += keytime
            # move to the back
            delete!(cache.entries, key)
            cache.entries[key] = entry
            entry, false
        end
    end

    if isowner
        try
            timed = @timed f()
            @lock cache.lock begin
                entry.buildtime = timed.time
                entry.compiletime = _compiletime(timed)
                entry.value = timed.value
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
        cache.evictions = 0
        foreach(values(cache.entries)) do e
            e.hits = 0
            e.keytime = 0.0
        end
    end
    cache
end

# All stats are read off the entries in the cache: every built entry was a miss, and the hits
# are summed over them. Entries which were dropped take their stats with them.
function print_stats(cache::ThreadsafeCache; n=10)
    @lock cache.lock begin
        built = filter(e -> !isnothing(e.value), collect(values(cache.entries)))
        misses = length(built)
        hits = sum(e -> e.hits, built; init=0)
        if isempty(built)
            println(cache.name, " Stats: no entries")
        else
            println(cache.name, " Stats:")
            println("  Entries:  ", misses, " (max ", cache.maxsize, "), built in ",
                    _fmt_time(sum(e -> e.buildtime, built)), ", of which ",
                    _fmt_time(sum(e -> e.compiletime, built)), " compiling")
            println("  Hits:     ", hits, " (", round(hits / (hits + misses) * 100, digits=1),
                    "% hit rate), keys took ", _fmt_time(sum(e -> e.keytime, built)))
        end

        withhits = sort!(filter(e -> e.hits > 0, built); by=e -> e.hits, rev=true)
        if !isempty(withhits)
            println("  Most hits (hits × build time):")
            for e in first(withhits, n)
                println("    ", lpad(e.hits, 5), " × ", rpad(_fmt_time(e.buildtime), 9), " ", e.label)
            end
            length(withhits) > n && println("    … and ", length(withhits) - n, " more")
        end
        nohits = filter(e -> e.hits == 0, built)
        if !isempty(nohits)
            labels = map(e -> e.label, first(nohits, n))
            more = length(nohits) > n ? ", …" : ""
            println("  Without hits (", length(nohits), "): ", join(labels, ", "), more)
        end
        if cache.evictions > 0
            println("  ", cache.evictions, " entries evicted at maxsize ", cache.maxsize,
                    ", consider a larger maxsize")
        end
    end
end
# Julia 1.10 doesn't report compile time
_compiletime(timed) = get(timed, :compile_time, 0.0)
_runtime(timed) = timed.time - _compiletime(timed)

function _fmt_time(t)
    str(x) = NetworkDynamics.str_significant(x; sigdigits=3)
    t < 1e-3 ? string(str(t * 1e6), " µs") : t < 1 ? string(str(t * 1e3), " ms") : string(str(t), " s")
end

# `generate_io_function`, reused for systems with the same structure. Callers read all
# per-instance data (defaults, guesses, metadata) from `sys`.
function generate_io_function_cached(sys, inputss, outputss; verbose=false, mtkcompile=false, kwargs...)
    build() = generate_io_function(sys, inputss, outputss; verbose, mtkcompile, kwargs...)
    # Verbose builds run fresh so that all their messages show up. The key only knows what ND's
    # own simplification reads, so `mtkcompile` builds are not cached.
    mode = MTK_CACHE_MODE[]
    if mode === :off || verbose || mtkcompile
        return build()
    end
    hit = true
    timedkey = @timed _structural_key(sys, inputss, outputss, kwargs)
    label() = _system_label(sys)
    gen = get!(COMPILE_CACHE, timedkey.value; label, keytime=timedkey) do
        hit = false
        build()
    end
    if hit && mode === :check
        rebuilt = build()
        diff = [k for k in keys(rebuilt) if !_samebuild(gen[k], rebuilt[k])]
        isempty(diff) || error("The MTK compile cache returned a different build for system $(getname(sys)), \
            differing in $(join(diff, ", ")). Please report this as a bug. \
            `NetworkDynamics.set_mtk_cache!(:off)` switches the cache off.")
    end
    # the containers which end up in the component are not shared with the cache
    fresh = map(copy, (; gen.obsrules, gen.aliasmap, gen.observed, gen.equations,
                         gen.full_equations, gen.outputeqs, gen.full_outputeqs, gen.mass_matrix))
    (; gen..., fresh...)
end

# The name of the first system with this structure, and its constructor if MTK recorded one.
function _system_label(sys)
    gui = ModelingToolkitBase.get_gui_metadata(sys)
    isnothing(gui) ? string(getname(sys)) : string(getname(sys), " (", gui.type.name, ")")
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
    NetworkDynamics.set_mtk_cache!(mode; maxsize=nothing)

Set the mode of the MTK cache. It has two parts: the compile cache lets a [`VertexModel`](@ref) or
[`EdgeModel`](@ref) reuse the generated functions of an earlier `System` with the same structure,
and the component cache lets constructors defined with [`@cached_component`](@ref) reuse the
systems they built.
- `:off` (default): everything is built from scratch. Switching off also empties both caches.
- `:on`: both caches are used.
- `:check`: like `:on`, but every cache hit is also built from scratch and compared with the cached
  result. Errors if they differ. This is slow and meant for debugging.

`true` and `false` stand for `:on` and `:off`. `maxsize` sets how many entries each of the two
caches keeps (500 by default).

See also [`NetworkDynamics.with_mtk_cache`](@ref).
"""
function NetworkDynamics.set_mtk_cache!(mode; maxsize=nothing)
    mode = cache_mode(mode)
    for cache in ALL_CACHES
        @lock cache.lock begin
            if !isnothing(maxsize)
                cache.maxsize = maxsize
            end
            mode === :off && empty!(cache.entries)
        end
    end
    MTK_CACHE_MODE[] = mode
end
cache_mode(use::Bool) = use ? :on : :off
function cache_mode(mode::Symbol)
    mode in (:off, :on, :check) || throw(ArgumentError("Cache mode must be :off, :on or :check, got :$mode."))
    mode
end

"""
    NetworkDynamics.wipe_mtk_cache!()

Empty both parts of the MTK cache (see [`NetworkDynamics.set_mtk_cache!`](@ref)). Needed when
the code generating models changes within a session, for example when editing NetworkDynamics or
a cached component with Revise.
"""
function NetworkDynamics.wipe_mtk_cache!()
    foreach(empty!, ALL_CACHES)
    nothing
end

"""
    NetworkDynamics.with_mtk_cache(f, mode=:on)

Run `f` with the MTK cache in the given mode (see [`NetworkDynamics.set_mtk_cache!`](@ref))
and restore the previous mode afterwards. Also resets the hit counts of
[`NetworkDynamics.mtk_cache_stats`](@ref). Since the stats describe the entries currently in the
cache, and switching the cache off empties it, print them inside the block.

```julia
NetworkDynamics.with_mtk_cache(false) do
    # code which builds mtk models without the cache
end
```

The cache is off by default. When it is on and a [`VertexModel`](@ref) or [`EdgeModel`](@ref) is built from
a `System` with the same equations, variables, bindings, init and guess formulas, subsystems
and inputs/outputs as an earlier one, the generated functions and equations are reused. Defaults,
guesses and other metadata are always read from the new system. Warnings of the build only show
up once, when the cache entry is filled. Builds with `verbose=true` or `mtkcompile=true` bypass
the cache. The cache is thread safe, so components can be built concurrently.

See also [`NetworkDynamics.set_mtk_cache!`](@ref), [`NetworkDynamics.wipe_mtk_cache!`](@ref).
"""
function NetworkDynamics.with_mtk_cache(f, mode=:on)
    oldmode = MTK_CACHE_MODE[]
    foreach(reset_counts!, ALL_CACHES)
    NetworkDynamics.set_mtk_cache!(mode)
    try
        f()
    finally
        NetworkDynamics.set_mtk_cache!(oldmode)
    end
end

"""
    NetworkDynamics.mtk_cache_stats()

Print statistics about both parts of the MTK cache: the entries with their build times and hits,
and the entries which were never hit. The stats describe the entries currently in the cache, so
print them inside a [`NetworkDynamics.with_mtk_cache`](@ref) block.
"""
function NetworkDynamics.mtk_cache_stats()
    foreach(print_stats, ALL_CACHES)
end

"""
    NetworkDynamics.mtk_cache_mode()

Return the current mode of the MTK cache, `:off`, `:on` or `:check`. See
[`NetworkDynamics.set_mtk_cache!`](@ref).
"""
NetworkDynamics.mtk_cache_mode() = MTK_CACHE_MODE[]
