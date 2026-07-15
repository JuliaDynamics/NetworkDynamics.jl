"""
    const AliasMap = Dict{Symbol, Tuple{Float64, Symbol}}

Maps an *alias* symbol to a `(factor, canonical)` pair with the semantics

    value(alias) == factor * value(canonical)

Components compiled from ModelingToolkit undergo alias elimination: many user-visible
symbols survive only as *observed* symbols which are pure (possibly sign-flipped or
scaled) aliases of a settable symbol. The aliasmap records that relationship so that the
initialization pipeline can canonicalize user input regardless of which member of an
alias class it was written against.

Canonical symbols are always settable symbols of the component, i.e. states, parameters,
inputs or outputs. Identity entries (`s => (1.0, s)`) are never stored; absence of a key
means "already canonical, or not an alias".

Stored as component metadata under the key `:aliasmap`, see [`set_aliasmap!`](@ref).
"""
const AliasMap = Dict{Symbol, Tuple{Float64, Symbol}}

"""
    settable_symbols(c::ComponentModel)

Returns the `Set` of symbols which can be *set* on the component, i.e. states,
parameters, inputs and outputs. Observables are excluded: they are computed, not set.

This is the universe of valid canonical symbols for an `AliasMap`.
"""
function settable_symbols(c::ComponentModel)
    s = Set{Symbol}()
    union!(s, sym(c))
    union!(s, psym(c))
    hasinsym(c) && union!(s, insym_flat(c))
    union!(s, outsym_flat(c))
    s
end

"""
    assert_aliasmap_compat(c::ComponentModel, am::AliasMap)

Validates an [`AliasMap`](@ref) against a component: factors must be finite and nonzero,
every canonical target must be settable, and no alias key may itself be settable (a
settable symbol must never be recorded as an alias of another settable symbol).

Returns `am` on success, throws an `ArgumentError` otherwise.
"""
function assert_aliasmap_compat(c::ComponentModel, am::AliasMap)
    settable = settable_symbols(c)
    for (alias, (factor, canonical)) in am
        if !isfinite(factor) || iszero(factor)
            throw(ArgumentError("AliasMap factor for :$alias must be finite and nonzero, got $factor."))
        end
        if canonical ∉ settable
            throw(ArgumentError("AliasMap maps :$alias to :$canonical, which is not a settable \
                                 symbol of the component model."))
        end
        if alias ∈ settable
            throw(ArgumentError("AliasMap key :$alias is itself a settable symbol of the \
                                 component model and must not be aliased to :$canonical."))
        end
    end
    am
end

# Tolerances for comparing two values which land on the same canonical symbol.
const ALIAS_RTOL = 1e-10
const ALIAS_ATOL = 1e-12

"""
    canonicalize(am::AliasMap, s::Symbol)::Tuple{Float64, Symbol}

Resolves `s` to its canonical settable symbol, returning `(factor, canonical)` with

    value(s) == factor * value(canonical)

Symbols which are not alias keys pass through as `(1.0, s)`. That includes symbols which
are already canonical, observables which are not pure aliases, and symbols unknown to the
component — canonicalization never validates the symbol universe.

A single lookup suffices: alias chains are resolved transitively when the map is
extracted, so a canonical symbol is never itself an alias key.
"""
canonicalize(am::AliasMap, s::Symbol) = get(am, s, (1.0, s))

"""
    normalize_valuedict(am::AliasMap, d; what="value", verbose=false, io=stdout)

Moves every value written on an alias symbol onto its canonical symbol, i.e. `d[alias] = v`
becomes `d[canonical] = v / factor`. Used for the `default`, `guess` and `init` dicts of
the initialization pipeline.

Values on multiple members of one alias class must agree *after* transformation, i.e.
`:a => 1.0` and `:b => -1.0` merge silently if `a ~ -b`. Disagreement throws an
`ArgumentError`. `what` names the kind of value in messages.

Returns a new dict; `d` is never mutated. `nothing` values (removal markers) travel with
their key untransformed.

See also: [`normalize_bounds`](@ref), [`canonicalize`](@ref).
"""
function normalize_valuedict(am::AliasMap, d; what="value", verbose=false, io=stdout)
    _normalize_symdict(am, d, _transform_value, what, verbose, io)
end

"""
    normalize_bounds(am::AliasMap, d; verbose=false, io=stdout)

Like [`normalize_valuedict`](@ref), but for a dict of `(lb, ub)` bounds tuples: both
endpoints are divided by the factor, and a *negative* factor swaps them, since scaling by
a negative number flips the interval.

Bounds on multiple members of one alias class must agree in both endpoints after
transformation.
"""
function normalize_bounds(am::AliasMap, d; verbose=false, io=stdout)
    _normalize_symdict(am, d, _transform_bounds, "bound", verbose, io)
end

# Shared skeleton of the two `normalize_*` functions above: transform each aliased key onto
# its canonical symbol via `transform`, checking values which collide there for agreement.
#
# Iteration order matters for reproducibility: colliding values need only agree
# approximately, so *which* one survives (and which symbols an error names) would otherwise
# depend on hash order. Hence two passes — canonical entries first, aliases after in sorted
# order — which makes a value written on the canonical symbol itself always win, and sorted
# order decide between competing aliases.
function _normalize_symdict(am::AliasMap, d, transform, what, verbose, io)
    isempty(am) && return d

    res = empty(d)
    aliases = Symbol[]
    for (s, v) in d
        haskey(am, s) ? push!(aliases, s) : (res[s] = v)
    end
    isempty(aliases) && return res

    # canonical => symbol the value in `res` came from; absent means it came from the
    # canonical symbol itself in the pass above
    origin = Dict{Symbol,Symbol}()
    for s in sort!(aliases)
        factor, canonical = am[s]
        val = transform(d[s], factor)
        if haskey(res, canonical)
            other = get(origin, canonical, canonical)
            _assert_agreement(canonical, what, factor,
                              other => d[other] => res[canonical],
                              s => d[s] => val)
            continue # first writer wins, both agree anyways
        end
        res[canonical] = val
        origin[canonical] = s
        verbose && printstyled(io, " - Move $what :$s (=$(_valstring(d[s]))) → :$canonical \
                                    (=$(_valstring(val))) via factor $factor.\n")
    end
    res
end

_transform_value(v, factor) = isnothing(v) ? v : v / factor

function _transform_bounds(v, factor)
    isnothing(v) && return v
    lb, ub = v ./ factor
    factor > 0 ? (lb, ub) : (ub, lb) # a negative factor flips the interval
end

# each entry is `symbol => raw_value => transformed_value`
function _assert_agreement(canonical, what, factor, entry1, entry2)
    s1, (raw1, v1) = entry1
    s2, (raw2, v2) = entry2
    _agree(v1, v2) && return nothing
    throw(ArgumentError("Conflicting $what values in the alias class of :$canonical: \
        :$s1 = $(_valstring(raw1)) (→ $(_valstring(v1))) and :$s2 = $(_valstring(raw2)) \
        (→ $(_valstring(v2)) via factor $factor). Values written on members of one alias \
        class must agree after transformation to the canonical symbol."))
end

_agree(v1::Real, v2::Real) = isapprox(v1, v2; rtol=ALIAS_RTOL, atol=ALIAS_ATOL)
_agree(v1::Tuple, v2::Tuple) = all(splat(_agree), zip(v1, v2))
_agree(v1, v2) = isequal(v1, v2) # `nothing` markers and anything else exotic

_valstring(v::Real) = str_significant(v; sigdigits=5)
_valstring(v::Tuple) = "(" * join(_valstring.(v), ", ") * ")"
_valstring(v) = repr(v)
