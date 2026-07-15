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
