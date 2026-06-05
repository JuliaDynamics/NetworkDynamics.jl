####
#### Default-from (parameter default propagation)
####
# Parameters can carry `:default_from` metadata (see `ParameterDefaultFrom`) which
# causes them to take their *default* value from another parameter:
#
#   - `default_from = :busbar₊Vbase`          take it from a parameter of a sibling
#                                             subcomponent in the same component.
#   - `default_from = (:src, :busbar₊Vbase)`  take it from a parameter of the src/dst
#   - `default_from = (:dst, :busbar₊Vbase)`  vertex of an edge (resolved on Network
#                                             construction).
#
# The source spec is the *full* parameter name, exactly as it appears in `psym`
# of the source component (namespaced with `₊`). This is different from writing
# `Vbase = busbar.Vbase`: it only ever *copies* a default from a known, named
# parameter, it does not traverse/evaluate the model hierarchy.
#
# It only ever touches metadata defaults; it never modifies the flat parameter
# array and never throws on value conflicts. Structural misuse (default_from on a
# non-parameter, or src/dst on a vertex) is an error.

# find the source parameter `srcspec` (a full, namespaced parameter name) in `c`
function _find_default_from_source(c::ComponentModel, srcspec)
    srcspec in psym(c) ? srcspec : nothing
end

# emit a warning only the first time a given message is seen (per resolution run)
function _warn_once!(warned, msg)
    msg in warned && return nothing
    push!(warned, msg)
    @warn msg
    nothing
end

# A propagated default is tagged with `:default_from_value` metadata holding the
# value that was copied in. `default_from` "owns" a default (and may overwrite it)
# while there is no default yet, or the current default still equals that tag. The
# tag is self-validating: if the default was changed by other means (e.g.
# `set_default!` clears it, or a raw `symmetadata[p][:default] = …` leaves it
# mismatched) the parameter is demoted to an explicit default. This makes
# resolution idempotent and order-independent.
function _default_overwritable!(c::ComponentModel, p::Symbol)
    has_default(c, p) || return true             # nothing there yet -> we may set it
    has_metadata(c, p, :default_from_value) || return false  # explicit (user) default
    if isequal(get_metadata(c, p, :default_from_value), get_default(c, p))
        return true                              # still a propagated default
    else
        delete_metadata!(c, p, :default_from_value)  # stale tag -> demote to explicit
        return false
    end
end

# apply the merge rule. `default_from` (re)writes the default (and its
# `:default_from_value` tag) when it owns it; an explicit default is kept (warns
# on conflict if verbose).
function _apply_default_from!(c::ComponentModel, p::Symbol, v; verbose, warned)
    if _default_overwritable!(c, p)
        changed = !has_default(c, p) || !isequal(get_default(c, p), v)
        set_default!(c, p, v)                       # NB: clears the value tag …
        set_metadata!(c, p, :default_from_value, v) # … so re-tag with the copied value
        return changed
    elseif isequal(get_default(c, p), v)
        return false
    else
        verbose && _warn_once!(warned, "default_from: parameter :$p of $(c.name) has default \
            $(get_default(c, p)) which differs from source value $v. Keeping existing default.")
        return false
    end
end

# error on structural misuse of default_from metadata
function _validate_default_from(c::ComponentModel)
    pset = Set(psym(c))
    for (s, md) in symmetadata(c)
        haskey(md, :default_from) || continue
        if s ∉ pset
            throw(ArgumentError("default_from metadata on non-parameter symbol :$s in $(c.name): \
                it is only supported for parameters."))
        end
        spec = md[:default_from]
        if spec isa Tuple && c isa VertexModel
            throw(ArgumentError("default_from = $(repr(spec)) on vertex parameter :$s in $(c.name): \
                src/dst sources are only valid on edge parameters."))
        end
    end
    nothing
end

# single pass of local (same-component) resolution, returns number of updates
function _default_from_local_pass!(c::ComponentModel, warned; verbose)
    n = 0
    for p in psym(c)
        has_default_from(c, p) || continue
        spec = get_default_from(c, p)
        spec isa Tuple && continue  # cross-component (edge), resolved on Network level
        src = _find_default_from_source(c, spec)
        if isnothing(src)
            _warn_once!(warned, "default_from: could not resolve source parameter $(repr(spec)) \
                for parameter :$p in $(c.name). Skipping.")
            continue
        end
        has_default(c, src) || continue  # nothing to copy yet
        changed = _apply_default_from!(c, p, get_default(c, src); verbose, warned)
        if changed
            n += 1
        end
    end
    n
end

# resolve same-component default_from to convergence, returns number of updates
function _resolve_default_from_local!(c::ComponentModel, warned; verbose)
    total = 0
    while (n = _default_from_local_pass!(c, warned; verbose)) > 0
        total += n
    end
    total
end

"""
    resolve_default_from!(c::ComponentModel; verbose=false)
    resolve_default_from!(nw::Network; verbose=false)

Resolve `:default_from` metadata (see [`ParameterDefaultFrom`](@ref)) by copying
**default values** from the referenced source parameters. This is called
automatically when components and the `Network` are constructed.

On a [`VertexModel`](@ref)/[`EdgeModel`](@ref) only *same-component* sources
(`default_from = :compname₊param`) are resolved. On a [`Network`](@ref)
cross-component (`default_from = (:src, …)` / `(:dst, …)`) sources from an edge's
src/dst vertex are resolved as well. Chains propagate in a single call.

A default you set manually is never overwritten (a conflict warns only when
`verbose`). A missing source parameter always warns. Structural misuse throws:
`:default_from` on a non-parameter symbol, or `(:src, …)`/`(:dst, …)` on a vertex
parameter.

Returns the number of defaults that were set.

See also [`ParameterDefaultFrom`](@ref).
"""
function resolve_default_from!(c::ComponentModel; verbose=false)
    _validate_default_from(c)
    _resolve_default_from_local!(c, Set{String}(); verbose)
end

# resolve cross-component (src/dst) sources on all edges, returns # of updates.
# an edge model instance is shared if it backs more than one edge (alias group of
# size > 1, as tracked by the index manager)
function _edge_is_aliased(im, ef)
    haskey(im.aliased_edgems, ef) && length(im.aliased_edgems[ef].idxs) > 1
end

function _resolve_default_from_crosscomponent!(nw::Network, warned; verbose)
    im = nw.im
    total = 0
    for i in eachindex(im.edgem)
        ef = im.edgem[i]
        any(p -> has_default_from(ef, p) && get_default_from(ef, p) isa Tuple, psym(ef)) || continue
        # a shared instance cannot hold per-edge values, so reject it up front
        # (also prevents the fixed-point loop from oscillating when two edges
        # would write different values to the same object).
        if _edge_is_aliased(im, ef)
            throw(ArgumentError("Cross-component (`:src`/`:dst`) `default_from` on an \
                edge model that is shared across multiple edges is not supported, because each \
                edge may need different source values (e.g. a transformer template reused for \
                several lines). Construct the `Network` with `dealias=true` (or pass per-edge \
                `copy`s of the edge model) so each edge gets its own instance."))
        end
        e = im.edgevec[i]
        for p in psym(ef)
            has_default_from(ef, p) || continue
            spec = get_default_from(ef, p)
            spec isa Tuple || continue
            dir, srcspec = spec
            vidx = if dir === :src
                e.src
            elseif dir === :dst
                e.dst
            else
                throw(ArgumentError("default_from: unknown direction $(repr(dir)) for parameter \
                    :$p of edge $i. Must be :src or :dst."))
            end
            vm = im.vertexm[vidx]
            src = _find_default_from_source(vm, srcspec)
            if isnothing(src)
                _warn_once!(warned, "default_from: could not resolve source parameter $(repr(srcspec)) \
                    on $dir vertex $vidx for parameter :$p of edge $i. Skipping.")
                continue
            end
            has_default(vm, src) || continue
            _apply_default_from!(ef, p, get_default(vm, src); verbose, warned) && (total += 1)
        end
    end
    total
end

function resolve_default_from!(nw::Network; verbose=false)
    for c in nw.im.vertexm
        _validate_default_from(c)
    end
    for c in nw.im.edgem
        _validate_default_from(c)
    end
    warned = Set{String}()
    total = 0
    while true
        # alternate local and cross-component resolution until a fixed point so
        # that cyclic deps between the two are resolved (cross-component changes
        # may enable further same-component resolution and vice versa)
        n = 0
        for c in nw.im.vertexm
            n += _resolve_default_from_local!(c, warned; verbose)
        end
        for c in nw.im.edgem
            n += _resolve_default_from_local!(c, warned; verbose)
        end
        n += _resolve_default_from_crosscomponent!(nw, warned; verbose)
        total += n
        n == 0 && break
    end
    total
end
