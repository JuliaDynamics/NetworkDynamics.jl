"""
    resolve_default_from(nw, additional_initformula, default_overrides; verbose=false)

Pre-pass for network initialization: scan every component for `:default_from` parameter metadata,
copy the source value (a composition-layer `default_overrides` entry wins over the source's
metadata default), bake it into a *weak* [`InitFormula`](@ref), and merge that formula into the
per-component `additional_initformula` map. Returns the (possibly augmented)
`additional_initformula`; if no component declares `default_from`, the argument is returned
unchanged.

The source value is read only from parameters/defaults — never from solved state — so the pass
runs before any per-component init task and is order-independent.

Failure modes: a `:default_from` on a non-parameter, a wrong scope for the component type
(`:src`/`:dst` on a vertex, `:hub` on an edge or on a non-injector vertex), or an unresolvable
source symbol all **error**. A resolvable source that simply carries no value yet is **skipped**
(no formula built), leaving the target parameter's own default in place.
"""
function resolve_default_from(nw::Network, additional_initformula, default_overrides; verbose=false)
    im = nw.im
    # cheap early-out: nothing to do unless some component declares default_from
    (any(_has_default_from, im.vertexm) || any(_has_default_from, im.edgem)) || return additional_initformula

    merged = _as_initformula_dict(additional_initformula)
    for vidx in eachindex(im.vertexm)
        _has_default_from(im.vertexm[vidx]) || continue
        _merge_formulas!(merged, VIndex(vidx),
            _default_from_formulas(nw, im.vertexm[vidx], vidx, false, default_overrides; verbose))
    end
    for eidx in eachindex(im.edgem)
        _has_default_from(im.edgem[eidx]) || continue
        _merge_formulas!(merged, EIndex(eidx),
            _default_from_formulas(nw, im.edgem[eidx], eidx, true, default_overrides; verbose))
    end
    merged
end

_has_default_from(c::ComponentModel) = any(md -> haskey(md, :default_from), values(symmetadata(c)))

# Build the baked weak InitFormulas for one component. Structural misuse errors here.
function _default_from_formulas(nw, c, cidx, is_edge, default_overrides; verbose)
    _validate_default_from_targets(c)
    formulas = InitFormula[]
    for p in psym(c)
        has_default_from(c, p) || continue
        spec = get_default_from(c, p)
        if !(spec isa Tuple && length(spec) == 2)
            throw(ArgumentError(
                "default_from on :$p of $(c.name) must be a `(scope, srcsym)` tuple with \
                 scope ∈ (:src, :dst, :hub), got $(repr(spec))."))
        end
        dir, srcsym = spec
        srcvidx = _resolve_default_from_source_vidx(nw, cidx, is_edge, dir, p, c.name)
        srcvm = nw.im.vertexm[srcvidx]
        if srcsym ∉ psym(srcvm)
            throw(ArgumentError(
                "default_from = $(repr(spec)) on :$p of $(c.name): source parameter :$srcsym \
                 not found on $dir vertex $srcvidx ($(srcvm.name)). Available parameters: \
                 $(sort(psym(srcvm)))."))
        end

        # Two-step read: a composition-layer `default_overrides` entry for the source index wins
        # over the source's metadata default (skip if neither carries a value yet). NB: overrides
        # are not alias-normalized (a known v1 gap) — the key must match `srcsym` exactly.
        override_idx = VIndex(srcvidx, srcsym)
        val = if default_overrides isa AbstractDict && haskey(default_overrides, override_idx)
            default_overrides[override_idx]
        elseif has_default(srcvm, srcsym)
            get_default(srcvm, srcsym)
        else
            verbose && printstyled(" - default_from: :$p of $(c.name) skipped — source :$srcsym \
                on $dir vertex $srcvidx has no value to copy yet\n")
            continue
        end

        # Bake the copied constant into a weak single-output InitFormula. The setter is built by
        label = "default_from($(repr(dir)), $(repr(srcsym)))"
        push!(formulas, InitFormula(_default_from_setter(p, val), [p], Symbol[]; weak=true, label))
    end
    formulas
end
_default_from_setter(p, val) = (out, u) -> (out[p] = val; nothing)

# `default_from` is only meaningful on a parameter (it weakly defaults it); anywhere else is a
# mistake that would otherwise be silently ignored.
function _validate_default_from_targets(c)
    pset = Set(psym(c))
    for (s, md) in symmetadata(c)
        haskey(md, :default_from) || continue
        s ∈ pset || throw(ArgumentError(
            "default_from metadata on non-parameter symbol :$s in $(c.name): it is only \
             supported for parameters."))
    end
    nothing
end

# Resolve the scope keyword to the source *vertex* index. src/dst read an edge's endpoints; hub
# reads the vertex an injector hangs off (the dst end of its LoopbackConnection edge). A wrong
# scope for the component type is an error.
function _resolve_default_from_source_vidx(nw, cidx, is_edge, dir, p, cname)
    im = nw.im
    if is_edge
        dir === :src && return im.edgevec[cidx].src
        dir === :dst && return im.edgevec[cidx].dst
        dir === :hub && throw(ArgumentError(
            "default_from = (:hub, …) on edge parameter :$p of $cname: :hub is only valid on \
             vertex parameters (injector ← hub); use :src/:dst on edge parameters."))
    else
        if dir === :hub
            is_injector(im, cidx) || throw(ArgumentError(
                "default_from = (:hub, …) on vertex $cidx ($cname): :hub sources are only valid \
                 on injector nodes (the src of a LoopbackConnection edge), but this vertex is \
                 not an injector."))
            for i in eachindex(im.edgem)
                is_loopback(im.edgem[i]) && im.edgevec[i].src == cidx && return im.edgevec[i].dst
            end
            error("internal: no loopback edge found for injector vertex $cidx")
        end
        (dir === :src || dir === :dst) && throw(ArgumentError(
            "default_from = ($(repr(dir)), …) on vertex parameter :$p of $cname: only :hub \
             sources are valid on vertex parameters; :src/:dst are only valid on edge parameters."))
    end
    throw(ArgumentError("default_from: unknown scope $(repr(dir)) for :$p of $cname. \
        Must be one of :src, :dst, :hub."))
end

_as_initformula_dict(::Nothing) = Dict{SymbolicIndex,Any}()
_as_initformula_dict(d::AbstractDict) = Dict{SymbolicIndex,Any}(d)

function _merge_formulas!(merged, idx, formulas)
    isempty(formulas) && return merged
    merged[idx] = collect_formulas(get(merged, idx, nothing), formulas)
    merged
end
