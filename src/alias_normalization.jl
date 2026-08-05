"""
    const AliasMap = Dict{Symbol, Symbol}

Maps an *alias* symbol to its canonical symbol, with the semantics

    value(alias) == value(canonical)

An alias class is several *names for one variable*, so canonicalizing is a rename and never
touches a value.

Components compiled from ModelingToolkit undergo alias elimination: many user-visible
symbols survive only as *observed* symbols which are pure aliases of a settable symbol. The
aliasmap records that relationship so that the initialization pipeline can canonicalize user
input regardless of which member of an alias class it was written against.

A canonical symbol is a settable symbol of the component — a state, parameter, input or
output — whenever the alias class has one; it *must* be that member, else
[`normalize_valuedict`](@ref) would move a default onto a symbol with no slot behind it. A
class with no settable member (two names for one observable, e.g. an output port wired to an
input port) canonicalizes onto its *terminal observable* instead: nothing can be stored
there, but the names unify, so a formula writing one member and one reading another meet on
the same symbol.

Stored as component metadata under the key `:aliasmap`, see [`set_aliasmap!`](@ref).
"""
const AliasMap = Dict{Symbol, Symbol}

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

Validates an [`AliasMap`](@ref) against a component: every canonical target must be settable
or an observable of the component, no canonical may itself be an alias key (a class has one
representative, so one lookup must suffice), and no alias key may itself be settable (a
settable symbol must never be recorded as an alias of another settable symbol).

Returns `am` on success, throws an `ArgumentError` otherwise.
"""
function assert_aliasmap_compat(c::ComponentModel, am::AliasMap)
    settable = settable_symbols(c)
    obs = obssym(c)
    for (alias, canonical) in am
        if canonical ∉ settable && canonical ∉ obs
            throw(ArgumentError("AliasMap maps :$alias to :$canonical, which is neither a settable \
                                 symbol nor an observable of the component model."))
        end
        if canonical ∉ settable && haskey(am, canonical)
            throw(ArgumentError("AliasMap maps :$alias to the observable :$canonical, which is \
                                 itself an alias of :$(am[canonical]). An observable canonical \
                                 must be terminal, i.e. the end of its alias chain."))
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
    canonicalize(am::AliasMap, s::Symbol)::Symbol

Resolves `s` to the canonical name of its alias class, i.e. the one symbol of the class that
has a storage slot behind it (or, failing that, its terminal observable).

Symbols which are not alias keys pass through unchanged. That includes symbols which are
already canonical, observables which are not aliases, and symbols unknown to the component —
canonicalization never validates the symbol universe.

A single lookup suffices: a class is consolidated onto one representative, so a canonical
symbol is never itself an alias key.
"""
canonicalize(am::AliasMap, s::Symbol) = get(am, s, s)

"""
    normalize_valuedict(am::AliasMap, d; what=:value, on_conflict=:error, verbose=false, io=stdout)

Moves every value written on an alias symbol onto its canonical symbol, i.e. `d[alias] = v`
becomes `d[canonical] = v`. Used for the `default`, `guess` and `init` dicts of the
initialization pipeline.

Values written on several members of one class must agree, since the members are one
variable. `what` names the kind of value in messages. `on_conflict` decides what happens when
two members disagree:

- `:error` (default) — throw an `ArgumentError`. Correct for asserted values (defaults):
  two inconsistent defaults on one variable are a model contradiction.
- `:keepfirst` — keep the established value and drop the conflicting one, noting it under
  `verbose`. Correct for solver seeds (guesses): the class only needs *a* starting point,
  and the deterministic winner is the value on the canonical symbol, else the first alias
  in sorted order.

Returns a new dict; `d` is never mutated. A `nothing` value (a removal marker) moves onto
the canonical symbol like any other and takes part in the collision check as an ordinary
value, so a marker and a real value on one class disagree: setting and removing a class in a
single dict is a contradiction.

See also: [`normalize_bounds`](@ref), [`canonicalize`](@ref).
"""
function normalize_valuedict(am::AliasMap, d; what::Symbol=:value, on_conflict::Symbol=:error, verbose=false, io=stdout)
    _normalize_symdict(am, d, what, on_conflict, verbose, io)
end

"""
    normalize_bounds(am::AliasMap, d; verbose=false, io=stdout)

Like [`normalize_valuedict`](@ref), but for a dict of `(lb, ub)` bounds tuples. Bounds on
multiple members of one alias class must agree in both endpoints.
"""
function normalize_bounds(am::AliasMap, d; verbose=false, io=stdout)
    _normalize_symdict(am, d, :bound, :error, verbose, io)
end

# Shared skeleton of the two `normalize_*` functions above: move each aliased key onto its
# canonical symbol, checking values which collide there for agreement.
#
# Iteration order matters for reproducibility: colliding values need only agree
# approximately, so *which* one survives (and which symbols an error names) would otherwise
# depend on hash order. Hence two passes — canonical entries first, aliases after in sorted
# order — which makes a value written on the canonical symbol itself always win, and sorted
# order decide between competing aliases.
function _normalize_symdict(am::AliasMap, d, what, on_conflict, verbose, io)
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
    moves = String[]   # rows for the verbose block, one per de-aliased symbol
    ndrop = 0
    for s in sort!(aliases)
        canonical = am[s]
        val = d[s]
        if haskey(res, canonical)
            other = get(origin, canonical, canonical)
            if !_agree(res[canonical], val)
                # asserted values (defaults/bounds) must not disagree; solver seeds
                # (guesses) may — keep the established one and drop this member.
                on_conflict === :error && _conflict_error(canonical, what,
                                                          other => res[canonical], s => val)
                ndrop += 1
                verbose && push!(moves, ":$s &⇒ :$canonical &($(_valstring(val))) \
                                         &dropped, class holds $(_valstring(res[canonical]))")
            end
            continue # first writer wins
        end
        res[canonical] = val
        origin[canonical] = s
        verbose && push!(moves, ":$s &⇒ :$canonical &($(_valstring(val))) &")
    end
    verbose && print_aligned_group(io, _dealias_title(what, ndrop), moves)
    res
end

function _dealias_title(what, ndrop)
    kind = what === :default ? "defaults" :
           what === :guess   ? "guesses"  :
           what === :bound   ? "bounds"   : "$(what)s"
    ndrop > 0 ? "De-aliased $kind ($ndrop conflicting dropped):" : "De-aliased $kind:"
end

# each entry is `symbol => value`; called only once the two values are known to disagree
function _conflict_error(canonical, what, entry1, entry2)
    s1, v1 = entry1
    s2, v2 = entry2
    throw(ArgumentError("Conflicting $what values in the alias class of :$canonical: \
        :$s1 = $(_valstring(v1)) and :$s2 = $(_valstring(v2)). Members of one alias class \
        are one variable, so values written on them must agree."))
end

# NaN counts as equal to NaN here: `resolve_rules` asks `!_agree` whether a value moved, and
# without this clause a NaN would look like it keeps changing forever.
_agree(v1::Real, v2::Real) = isapprox(v1, v2; rtol=ALIAS_RTOL, atol=ALIAS_ATOL) || (isnan(v1) && isnan(v2))
_agree(v1::Tuple, v2::Tuple) = all(splat(_agree), zip(v1, v2))
_agree(v1, v2) = isequal(v1, v2) # `nothing` markers and anything else exotic

_valstring(v::Real) = str_significant(v; sigdigits=5, phantom_minus=true)
_valstring(v::Tuple) = "(" * join(_valstring.(v), ", ") * ")"
_valstring(v) = repr(v)

