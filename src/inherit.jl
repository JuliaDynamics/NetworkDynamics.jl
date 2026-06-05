####
#### Parameter inheritance
####
# Parameters can carry `:inherit` metadata (see `ParameterInherit`) which causes
# them to automatically take their *default* value from another component:
#
#   - `inherit = :busbar₊Vbase`          inherit from a parameter of a sibling
#                                        subcomponent in the same component.
#   - `inherit = (:src, :busbar₊Vbase)`  inherit from a parameter of the src/dst
#   - `inherit = (:dst, :busbar₊Vbase)`  vertex of an edge (resolved on Network
#                                        construction).
#
# The source spec is the *full* parameter name, exactly as it appears in `psym`
# of the source component (namespaced with `₊`). This is different from writing
# `Vbase = busbar.Vbase`: inheritance only ever *seeds* a default from a known,
# named parameter, it does not traverse/evaluate the model hierarchy.
#
# Inheritance only ever touches metadata defaults; it never modifies the flat
# parameter array and never throws on value conflicts. Structural misuse (inherit
# on a non-parameter, or src/dst inheritance on a vertex) is an error.

# find the source parameter `srcspec` (a full, namespaced parameter name) in `c`
function _find_inherit_source(c::ComponentModel, srcspec)
    srcspec in psym(c) ? srcspec : nothing
end

# emit a warning only the first time a given message is seen (per resolution run)
function _warn_once!(warned, msg)
    msg in warned && return nothing
    push!(warned, msg)
    @warn msg
    nothing
end

# An inherited default is tagged with an `:inherited` metadata holding the value
# that was inherited. Inheritance "owns" a default (and may overwrite it) while
# there is no default yet, or the current default still equals its `:inherited`
# tag. The tag is self-validating: if the default was changed by other means (e.g.
# `set_default!` clears it, or a raw `symmetadata[p][:default] = …` leaves it
# mismatched) the parameter is demoted to an explicit default. This makes
# resolution idempotent and order-independent.
function _default_overwritable!(c::ComponentModel, p::Symbol)
    has_default(c, p) || return true             # nothing there yet -> we may set it
    has_metadata(c, p, :inherited) || return false  # explicit (user) default
    if isequal(get_metadata(c, p, :inherited), get_default(c, p))
        return true                              # still an inherited default
    else
        delete_metadata!(c, p, :inherited)       # stale tag -> demote to explicit
        return false
    end
end

# apply the merge rule. Inheritance (re)writes the default (and its `:inherited`
# tag) when it owns it; an explicit default is kept (warns on conflict if verbose).
function _apply_inherit!(c::ComponentModel, p::Symbol, v; verbose, warned)
    if _default_overwritable!(c, p)
        changed = !has_default(c, p) || !isequal(get_default(c, p), v)
        set_default!(c, p, v)               # NB: clears the `:inherited` tag …
        set_metadata!(c, p, :inherited, v)  # … so re-tag with the inherited value
        return changed
    elseif isequal(get_default(c, p), v)
        return false
    else
        verbose && _warn_once!(warned, "inherit: parameter :$p of $(c.name) has default \
            $(get_default(c, p)) which differs from inherited value $v. Keeping existing default.")
        return false
    end
end

# error on structural misuse of inherit metadata
function _validate_inherit(c::ComponentModel)
    pset = Set(psym(c))
    for (s, md) in symmetadata(c)
        haskey(md, :inherit) || continue
        if s ∉ pset
            throw(ArgumentError("inherit metadata on non-parameter symbol :$s in $(c.name): \
                inheritance is only supported for parameters."))
        end
        spec = md[:inherit]
        if spec isa Tuple && c isa VertexModel
            throw(ArgumentError("inherit = $(repr(spec)) on vertex parameter :$s in $(c.name): \
                src/dst inheritance is only valid on edge parameters."))
        end
    end
    nothing
end

# single pass of local (same-component) inheritance, returns number of updates
function _inherit_local_pass!(c::ComponentModel, warned; verbose)
    n = 0
    for p in psym(c)
        has_inherit(c, p) || continue
        spec = get_inherit(c, p)
        spec isa Tuple && continue  # cross-component (edge), resolved on Network level
        src = _find_inherit_source(c, spec)
        if isnothing(src)
            _warn_once!(warned, "inherit: could not resolve source parameter $(repr(spec)) \
                for parameter :$p in $(c.name). Skipping.")
            continue
        end
        has_default(c, src) || continue  # nothing to inherit yet
        changed = _apply_inherit!(c, p, get_default(c, src); verbose, warned)
        if changed
            n += 1
        end
    end
    n
end

# resolve same-component inheritance to convergence, returns number of updates
function _inherit_local!(c::ComponentModel, warned; verbose)
    total = 0
    while (n = _inherit_local_pass!(c, warned; verbose)) > 0
        total += n
    end
    total
end

"""
    inherit_parameters!(c::ComponentModel; verbose=false)
    inherit_parameters!(nw::Network; verbose=false)

Resolve parameter inheritance (`:inherit` metadata, see [`ParameterInherit`](@ref))
by copying **default values** from the referenced source parameters. This is called
automatically when components and the `Network` are constructed.

On a [`VertexModel`](@ref)/[`EdgeModel`](@ref) only *same-component* inheritance
(`inherit = :compname₊param`) is resolved. On a [`Network`](@ref) cross-component
(`inherit = (:src, …)` / `(:dst, …)`) inheritance from an edge's src/dst vertex is
resolved as well. Chains propagate in a single call.

A default you set manually is never overwritten (a conflict warns only when
`verbose`). A missing source parameter always warns. Structural misuse throws:
`:inherit` on a non-parameter symbol, or `(:src, …)`/`(:dst, …)` on a vertex
parameter.

Returns the number of defaults that were set.

See also [`ParameterInherit`](@ref).
"""
function inherit_parameters!(c::ComponentModel; verbose=false)
    _validate_inherit(c)
    _inherit_local!(c, Set{String}(); verbose)
end

# resolve cross-component (src/dst) inheritance on all edges, returns # of updates
# an edge model instance is shared if it backs more than one edge (alias group of
# size > 1, as tracked by the index manager)
function _edge_is_aliased(im, ef)
    haskey(im.aliased_edgems, ef) && length(im.aliased_edgems[ef].idxs) > 1
end

function _inherit_crosscomponent_parameters!(nw::Network, warned; verbose)
    im = nw.im
    total = 0
    for i in eachindex(im.edgem)
        ef = im.edgem[i]
        any(p -> has_inherit(ef, p) && get_inherit(ef, p) isa Tuple, psym(ef)) || continue
        # a shared instance cannot hold per-edge inherited values, so reject it up
        # front (also prevents the fixed-point loop from oscillating when two edges
        # would write different values to the same object).
        if _edge_is_aliased(im, ef)
            throw(ArgumentError("Cross-component (`:src`/`:dst`) parameter inheritance on an \
                edge model that is shared across multiple edges is not supported, because each \
                edge may need different inherited values (e.g. a transformer template reused for \
                several lines). Construct the `Network` with `dealias=true` (or pass per-edge \
                `copy`s of the edge model) so each edge gets its own instance."))
        end
        e = im.edgevec[i]
        for p in psym(ef)
            has_inherit(ef, p) || continue
            spec = get_inherit(ef, p)
            spec isa Tuple || continue
            dir, srcspec = spec
            vidx = if dir === :src
                e.src
            elseif dir === :dst
                e.dst
            else
                throw(ArgumentError("inherit: unknown direction $(repr(dir)) for parameter \
                    :$p of edge $i. Must be :src or :dst."))
            end
            vm = im.vertexm[vidx]
            src = _find_inherit_source(vm, srcspec)
            if isnothing(src)
                _warn_once!(warned, "inherit: could not resolve source parameter $(repr(srcspec)) \
                    on $dir vertex $vidx for parameter :$p of edge $i. Skipping.")
                continue
            end
            has_default(vm, src) || continue
            _apply_inherit!(ef, p, get_default(vm, src); verbose, warned) && (total += 1)
        end
    end
    total
end

function inherit_parameters!(nw::Network; verbose=false)
    for c in nw.im.vertexm
        _validate_inherit(c)
    end
    for c in nw.im.edgem
        _validate_inherit(c)
    end
    warned = Set{String}()
    total = 0
    while true
        # alternate local and cross-component resolution until a fixed point so
        # that cyclic deps between the two are resolved (cross-component changes
        # may enable further same-component inheritance and vice versa)
        n = 0
        for c in nw.im.vertexm
            n += _inherit_local!(c, warned; verbose)
        end
        for c in nw.im.edgem
            n += _inherit_local!(c, warned; verbose)
        end
        n += _inherit_crosscomponent_parameters!(nw, warned; verbose)
        total += n
        n == 0 && break
    end
    total
end
