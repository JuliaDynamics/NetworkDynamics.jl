####
#### Expansion of `@cached_component`
####

# Splits the constructor in two. The user's function, without its `@structural` line, becomes a
# plain build function. The constructor under the user's name also takes the defaults, computes
# the key and calls `load_cached_component`, which runs the build function on a miss.
#
#     @cached_component function Lag(; name, lag=true, T=1.0)
#         @structural iszero(T)
#         ...body...
#     end
#
# becomes roughly
#
#     function _lag_build(; name, lag=true, T=1.0)
#         ...body...
#     end
#     Base.@__doc__ function Lag(; name, lag=true, T=1.0, cache=NetworkDynamics.mtk_model_cache_mode(), _defaults...)
#         NetworkDynamics.load_cached_component(
#             () -> _lag_build(; name, lag, T);  # the builder function
#             key = (Lag, 0x5f3a..., (lag, iszero(T))),  # constructor, token, implicit and explicit structural
#             mode = cache,                              # :off, :on or :check
#             name,
#             kwargs = (; lag, T),                       # applied like defaults where they name a symbol
#             defaults = _defaults,
#             gui = Main.Lag
#         )
#     end
function NetworkDynamics.cached_component_expr(mod, fdef)
    def = MacroTools.splitdef(fdef)
    fname = def[:name]
    isempty(def[:args]) || error("@cached_component supports keyword arguments only, $fname has positional arguments.")
    isempty(get(def, :whereparams, ())) || error("@cached_component does not support `where` clauses in $fname.")
    haskey(def, :rtype) && error("@cached_component does not support a return type annotation in $fname.")

    kwargs = get(def, :kwargs, Any[])
    kwnames = Symbol[]
    for kw in kwargs
        kwname, _, slurp, _ = MacroTools.splitarg(kw)
        if slurp
            error("$fname slurps `$kwname...`. @cached_component takes the defaults itself and \
                applies them to the finished system, so remove them from the signature and drop the \
                call to `set_mtk_defaults`.")
        end
        push!(kwnames, kwname)
    end
    :name in kwnames || error("@cached_component constructor $fname needs a `name` keyword argument.")
    :cache in kwnames && error("`cache` is reserved by @cached_component, rename it in $fname.")

    # The key holds the kwargs which `@structural` doesn't mention, as they are, and the values of
    # the `@structural` expressions. A mentioned kwarg only enters the key through those.
    explicit_structural, body = _pop_structural(def[:body], fname)
    mentioned = filter(k -> _reads(Expr(:block, explicit_structural...), k), kwnames)
    :name in mentioned && error("`name` can't be structural in $fname.")
    implicit_structural = setdiff(kwnames, mentioned, (:name,))

    builddef = copy(def)
    buildname = Symbol("#", fname, "#build")
    builddef[:name] = buildname
    builddef[:body] = body

    # The defaults get a hidden name, so that the `@structural` expressions can't read them. A
    # random token changes whenever the definition is evaluated anew, e.g. by Revise, so that an
    # edited constructor does not hit entries of its old version.
    defaults = gensym(:defaults)
    key = :(($fname, $(rand(UInt)), ($(implicit_structural...), $(explicit_structural...))))
    cachekw = Expr(:kw, :cache, :($(GlobalRef(NetworkDynamics, :mtk_model_cache_mode))()))
    def[:kwargs] = Any[kwargs..., cachekw, :($defaults...)]
    def[:body] = quote
        $(GlobalRef(NetworkDynamics, :load_cached_component))(() -> $buildname(; $(kwnames...));
            key = $key,
            mode = cache,
            name,
            kwargs = (; $(setdiff(kwnames, (:name,))...)),
            defaults = $defaults,
            gui = $(QuoteNode(GlobalRef(mod, fname))))
    end
    quote
        $(MacroTools.combinedef(builddef))
        Base.@__doc__ $(MacroTools.combinedef(def))
    end
end

# Splits the `@structural` line off the body. Returns its expressions and the rest of the body.
# The line has to come first, so that its expressions can only see the kwargs.
function _pop_structural(body, fname)
    isstructural(ex) = Meta.isexpr(ex, :macrocall) && ex.args[1] === Symbol("@structural")
    body = Meta.isexpr(body, :block) ? copy(body) : Expr(:block, body)
    idx = findfirst(ex -> !(ex isa LineNumberNode), body.args)
    structural = Any[]
    if !isnothing(idx) && isstructural(body.args[idx])
        for ex in body.args[idx].args[2:end]
            # a `begin` block holds one expression per line
            exs = Meta.isexpr(ex, :block) ? ex.args : [ex]
            append!(structural, filter(x -> !(x isa LineNumberNode), exs))
        end
        deleteat!(body.args, idx)
    end
    MacroTools.postwalk(body) do ex
        isstructural(ex) && error("@structural must be the first statement in the body of $fname.")
        ex
    end
    structural, body
end

# Whether `ex` refers to the variable `sym`. Keyword names in calls like `f(; sym=1)` don't
# count, the shorthand `f(; sym)` does.
_reads(ex, sym) = ex === sym
function _reads(ex::Expr, sym)
    if Meta.isexpr(ex, :kw)
        _reads(ex.args[2], sym)
    else
        any(x -> _reads(x, sym), ex.args)
    end
end

####
#### Cache
####

# The constructor calls go into `COMPONENT_CACHE`, which is defined next to the compile cache.

# Called by the code `@cached_component` generates. `build()` runs the user's function, `key`
# holds the constructor, its token and the structural values, and `kwargs` the explicit kwargs.
function NetworkDynamics.load_cached_component(build; key, mode, name, kwargs, defaults, gui)
    mode = cache_mode(mode)
    modelname = nameof(first(key))
    if mode === :off
        return _finish(_built(build, gui), kwargs, defaults)
    end

    # values which are equal but of different type, like 1 and 1.0, are different keys
    typedkey = (key[1], key[2], map(v -> (typeof(v), v), key[3]))
    hit = true
    cached = get!(COMPONENT_CACHE, typedkey) do
        hit = false
        sys = _built(build, gui)
        _warn_overridden(sys, kwargs, modelname)
        sys
    end

    sys = cached
    hit && @reset sys.name = name
    sys = _finish(sys, kwargs, defaults)
    if hit && mode === :check
        _check_cached(sys, _finish(_built(build, gui), kwargs, defaults), modelname)
    end
    sys
end

# The user's function returns a plain system, which gets the GUI metadata `@component` would set.
function _built(build, gui)
    sys = build()
    if isnothing(ModelingToolkitBase.get_gui_metadata(sys))
        @reset sys.gui_metadata = ModelingToolkitBase.GUIMetadata(gui)
    end
    sys
end

# Applies the defaults to a built system. Explicit kwargs which name a symbol of the system count
# as defaults too, on every path, so that a hit is the same as a fresh build.
function _finish(sys, kwargs, defaults)
    applied = filter(k -> _names_symbol(sys, k), collect(keys(kwargs)))
    NetworkDynamics.set_mtk_defaults(sys, (; kwargs[Tuple(applied)]..., defaults...))
end

# An explicit kwarg whose value replaces a different value the body gave its symbol is probably
# a mistake, since the body's value never shows up.
function _warn_overridden(sys, kwargs, modelname)
    ics = initial_conditions(sys)
    for (k, v) in pairs(kwargs)
        _names_symbol(sys, k) || continue
        sym = getproperty_symbolic(sys, k; might_contain_toplevel_ns=false)
        built = haskey(ics, sym) ? unwrap_const(ics[sym]) : get(ModelingToolkitBase.get_bindings(sys), unwrap(sym), nothing)
        if !isequal(built, v)
            @warn "Keyword argument $k of cached component $modelname names a symbol of the system, \
                so it is applied as its default and replaces the value $built set in the body."
        end
    end
end

function _names_symbol(sys, name)
    sym = try
        getproperty_symbolic(sys, name; might_contain_toplevel_ns=false)
    catch
        return false
    end
    !(sym isa ModelingToolkitBase.AbstractSystem)
end

# Compares a cache hit with a fresh build of the same call: the structure, and the data which
# every instance has for itself. A difference means that the hit gave a wrong result.
function _check_cached(hit, fresh, modelname)
    diff = String[]
    if _structural_hash(hit, zero(UInt)) != _structural_hash(fresh, zero(UInt))
        push!(diff, "structure")
    elseif !_varmeta_isequal(hit, fresh)
        push!(diff, "variable metadata")
    end
    nameof(hit) === nameof(fresh) || push!(diff, "name")
    _dict_isequal(initial_conditions(hit), initial_conditions(fresh)) || push!(diff, "initial conditions")
    _dict_isequal(ModelingToolkitBase.guesses(hit), ModelingToolkitBase.guesses(fresh)) || push!(diff, "guesses")
    isempty(diff) || error("Cached component $modelname returned a cached system for $(nameof(fresh)) \
        which differs from a fresh build in its $(join(diff, ", ")). Probably a `@structural` \
        expression hides a structural choice, a keyword argument reaches the system other than by \
        name, or the code of a constructor changed since the cache was filled \
        (`wipe_mtk_model_cache!()`). `cache=:off` skips the cache for a single call.")
end
_dict_isequal(a, b) = isequal(Dict(a), Dict(b))

# Variable metadata compared in full, apart from the default value. A symbol keeps the default
# it was declared with, while the instance reads its defaults from the system.
function _varmeta_isequal(a, b)
    samemeta(x, y) = isequal(_meta_without_default(x), _meta_without_default(y))
    for get_vars in (ModelingToolkitBase.get_unknowns, ModelingToolkitBase.get_ps)
        _allpairs(samemeta, get_vars(a), get_vars(b)) || return false
    end
    _allpairs(_varmeta_isequal, ModelingToolkitBase.get_systems(a), ModelingToolkitBase.get_systems(b))
end
_allpairs(f, xs, ys) = length(xs) == length(ys) && all(i -> f(xs[i], ys[i]), eachindex(xs, ys))
function _meta_without_default(v)
    md = Symbolics.metadata(unwrap(v))
    isnothing(md) && return Dict()
    # look up each key, iterating would also return the values it shadows
    Dict(k => md[k] for (k, _) in md if k !== Symbolics.VariableDefaultValue)
end
