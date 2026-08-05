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
there, but the names unify, which is what lets a pin travel between them (see
[`pinned_obssyms`](@ref)).

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

"""
    generate_obs_expansion(cf::ComponentModel, syms::Vector{Symbol}; stop_at=Set{Symbol}()) -> (roots, f)

Expresses `syms` through the component's observed equations in terms of *settable* symbols.
Returns the `roots` the expansion bottoms out on and a closure

    f(rootvals::AbstractVector{<:Real}, t::Real)::Vector{Float64}

giving the values of `syms`, in order, from the values of `roots`, in order.

A symbol without an observed equation is its own root and passes through untouched. That
makes the whole input list of a formula expandable in one sweep, no matter which entries
happen to be observables — see [`normalize`](@ref).

Symbols in `stop_at` are treated as if they had no observed equation: they are their own
root, and expansion of other symbols does not substitute through their defining equation.
This is how pinned observables (observables written by an `InitFormula`) become readable
dataflow nodes instead of being expanded away — see [`pinned_obssyms`](@ref).

The method itself lives in the ModelingToolkit extension.
"""
function generate_obs_expansion end

# Only the MTK extension implements `generate_obs_expansion`. Without it there are no
# symbolic observed equations to expand, so name the fix rather than fail on dispatch.
function _obs_expansion(cf::ComponentModel, syms::Vector{Symbol}; stop_at=Set{Symbol}())
    if !hasmethod(generate_obs_expansion, Tuple{typeof(cf),Vector{Symbol}})
        throw(ArgumentError("Expanding the observed symbol(s) $(intersect(syms, obssym(cf))) \
            to their settable roots requires the ModelingToolkit extension. Load \
            ModelingToolkit (or ModelingToolkitBase) to enable it."))
    end
    generate_obs_expansion(cf, syms; stop_at)
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

# The NaN clause is not cosmetic: `resolve_rules` decides "has this value moved" with `!_agree`,
# and `isapprox(NaN, NaN)` is `false`. Without it a NaN inside a rule block would re-arm its
# readers on every pass and the block would never settle.
_agree(v1::Real, v2::Real) = isapprox(v1, v2; rtol=ALIAS_RTOL, atol=ALIAS_ATOL) || (isnan(v1) && isnan(v2))
_agree(v1::Tuple, v2::Tuple) = all(splat(_agree), zip(v1, v2))
_agree(v1, v2) = isequal(v1, v2) # `nothing` markers and anything else exotic

_valstring(v::Real) = str_significant(v; sigdigits=5, phantom_minus=true)
_valstring(v::Tuple) = "(" * join(_valstring.(v), ", ") * ")"
_valstring(v) = repr(v)

"""
    pinned_obssyms(formulas, cf::ComponentModel) -> Set{Symbol}
    pinned_obssyms(cf::ComponentModel; guess=false) -> Set{Symbol}

The set of *pinned observables* of a formula set: observable symbols some formula writes, in
canonical names. A formula output which is neither settable nor a pure alias of a settable
symbol has no storage slot behind it — but when a formula explicitly writes it, it becomes a
legitimate dataflow node of that initialization: the written value is held in the working
dicts while formulas run, downstream formulas read it directly instead of expanding through
its defining equation (see [`generate_obs_expansion`](@ref) and [`normalize`](@ref)). It is
never stored on the component.

The pin is collected under the *canonical* name of the written symbol, which matters when an
observable-only alias class holds several names for one node (an output port wired to an
input port). Expansion bottoms out on the canonical, so a formula writing one member and one
reading another meet on it — which of the two interchangeable names each side happened to be
written against must not decide whether the pin connects.

The two formula kinds pin with different strength, staged the way they run:

- an `InitFormula` pin is a *commitment*: the value lands in the working defaults and is
  checked post-solve against the value the final state actually produces.
- a `GuessFormula` pin is a *hint*: the value lands in the working guesses, seeds downstream
  guess formulas and the nonlinear solve, and is never checked — the solver may move away
  from it freely.

Since init formulas run before guess formulas, they must never depend on a guess writer:
init formulas are normalized against the init pin-set only, guess formulas against the
union of both. A guess reader of an init-pinned symbol finds the value through the
documented defaults-before-guesses lookup, and an init pin shadows a guess pin on the same
symbol the same way.

This is a static property of the formula *set*: it must be computed over all formulas of an
initialization before any of them is normalized, so that readers and writers of a pin agree
on the frontier no matter in which order they run.

The component method computes the frontier from the formulas attached to the component's
metadata: the init pin-set by default, the guess-formula frontier with `guess=true`.
"""
function pinned_obssyms(formulas, cf::ComponentModel)
    pins = Set{Symbol}()
    isnothing(formulas) && return pins
    am = get_aliasmap(cf)
    settable = settable_symbols(cf)
    obs = obssym(cf)
    for f in formulas
        for s in f.outsym
            s ∈ obs || continue
            c = canonicalize(am, s)
            c ∉ settable && push!(pins, c)
        end
    end
    pins
end
function pinned_obssyms(cf::ComponentModel; guess=false)
    pins = pinned_obssyms(has_initformula(cf) ? get_initformulas(cf) : nothing, cf)
    guess || return pins
    pins ∪ pinned_obssyms(has_guessformula(cf) ? get_guessformulas(cf) : nothing, cf)
end

"""
    normalize(f::InitFormula, am::AliasMap, cf::ComponentModel; t=NaN, pinned=Set{Symbol}())
    normalize(f::GuessFormula, am::AliasMap, cf::ComponentModel; t=NaN, pinned=Set{Symbol}())

Rewrites a formula so that it speaks in *settable* symbols only, without touching what the
user wrote: the returned formula carries the original in its `derived_from` field, and calls
the original's closure unchanged, under the symbol names it was written against.

Both ends are rewritten:

- **outputs** are canonicalized through `am`, so the value written travels untouched.
  Writing to a non-translatable symbol (an observable which is not a pure alias) is
  meaningless — the value of an observable *is* its expression — so that throws, unless the
  symbol is in `pinned`: a pinned observable is an init-time dataflow node a formula may
  legitimately write (see [`pinned_obssyms`](@ref), also for how the frontier differs
  between the two formula kinds).
- **inputs** are expanded to their settable roots via [`generate_obs_expansion`](@ref),
  unconditionally: observables are never storage, so a default pinned on one is a
  consistency claim about the model, not a value a formula may read. A pure alias input
  needs no special case, it is the degenerate expansion onto its root. Expansion stops at
  symbols in `pinned` — those are written by a formula of the same init, so they are read
  directly instead of being expanded through their defining equation.

The resulting input list (the roots) and canonical output list are what the topological
sort, the duplicate-writer detection and the skip logic downstream operate on, which is the
point of the exercise: a formula writing `:machine₊Efd` and one reading `:avr₊Efd` are one
dependency edge once both are canonical, and so are a formula writing `:x` and one reading
an observable `a = x + y` — or one pinning observable `:y` and one reading `:y`.

A formula with no aliased outputs and no expandable observable inputs is returned `===`.

`t` is passed to the expansion for observables which depend explicitly on time.

`InitConstraint`s need none of this: they are evaluated against a full candidate state where
the observable mapping makes every symbol readable already.
"""
function normalize(f::Union{InitFormula,GuessFormula}, am::AliasMap, cf::ComponentModel;
                   t=NaN, pinned=Set{Symbol}())
    canonout = _canonicalize_outputs(f, am, cf, pinned)
    # pinned obs inputs are read directly from the working dict, only the rest expands
    expandable = setdiff(intersect(f.sym, obssym(cf)), pinned)

    # nothing to rewrite: hand back the very same formula, and stay independent of the ext
    canonout == f.outsym && isempty(expandable) && return f

    # only expand when there is something to expand: an aliased output alone must not drag
    # in the MTK extension
    roots, expand = isempty(expandable) ? (copy(f.sym), nothing) : _obs_expansion(cf, f.sym; stop_at=pinned)
    _assert_no_self_dependency(f, roots, canonout, cf, expandable, pinned)

    wrapped = function (out, u)
        rootvals = Float64[u[r] for r in roots]
        invals = if isnothing(expand)
            rootvals
        else
            vals = expand(rootvals, t)
            _assert_expansion_resolved(vals, f.sym, roots, t)
            vals
        end
        tmp = SymbolicView(zeros(length(f.outsym)), f.outsym)
        f(tmp, invals) # the user's closure, under the names the user wrote
        for (i, s) in enumerate(f.outsym)
            out[canonout[i]] = tmp[s]
        end
        nothing
    end
    # Rebuild carrying the original as `derived_from`: swap in the rewritten closure/outputs/
    # inputs; every other field (prettyprint, and InitFormula's `weak`) rides along untouched.
    # `setproperties` covers both formula types and still runs the self-dependency check.
    setproperties(f, (; f=wrapped, outsym=canonout, sym=roots, derived_from=f))
end

_formulatype(::InitFormula) = InitFormula
_formulatype(::GuessFormula) = GuessFormula

function _canonicalize_outputs(f, am, cf, pinned=Set{Symbol}())
    canonout = [canonicalize(am, s) for s in f.outsym]

    # a pinned obs canonical is a legitimate target: it IS the dataflow node, even though
    # there is no storage slot behind it
    settable = settable_symbols(cf)
    bad = [s => c for (s, c) in zip(f.outsym, canonout) if c ∉ settable && c ∉ pinned]
    if !isempty(bad)
        throw(ArgumentError("$(_formulatype(f)) cannot write to $(first.(bad)): \
            $(length(bad) == 1 ? "it is an observable" : "they are observables") which \
            $(length(bad) == 1 ? "is" : "are") not a pure alias of a settable symbol. The \
            value of such an observable is its expression over other symbols, so there is \
            no slot to write it to. Target the underlying settable symbol(s) instead."))
    end

    if !allunique(canonout)
        dupes = unique(c for c in canonout if count(isequal(c), canonout) > 1)
        throw(ArgumentError("$(_formulatype(f)) writes $(f.outsym), which collapses onto \
            the canonical symbol(s) $dupes — the same variable would be written twice. \
            Members of one alias class are one variable."))
    end
    canonout
end

# The raw lists can hide a self-dependency which only normalization makes visible: `:θ` and
# `:u_r` are one variable, and reading an observable `a = x + y` is reading `:x`. The
# formula constructor would catch it too, but only knows the canonical names.
function _assert_no_self_dependency(f, roots, canonout, cf, expandable, pinned)
    self_deps = intersect(roots, canonout)
    isempty(self_deps) && return nothing
    # one fact per line: the symbol lists get long in deeply nested models, and a reader who
    # has to scan them wants an obvious place to jump to next
    lines = ["$(_formulatype(f)) $(_formula_ref(f)) depends on its own output:",
             " - writes $(_symrepr(f.outsym))" *
                 (canonout == f.outsym ? "" : ", i.e. $(_symrepr(canonout))"),
             " - reads $(_symrepr(f.sym))" *
                 (roots == f.sym ? "" : ", expanded to $(_symrepr(roots))")]
    append!(lines, _self_dependency_lines(f, cf, expandable, pinned, roots, self_deps))
    throw(ArgumentError(join(lines, "\n")))
end

function _self_dependency_lines(f, cf, expandable, pinned, roots, self_deps)
    stop_at = pinned ∪ setdiff(roots, settable_symbols(cf))
    lines = String[]
    fixes = String[]
    for s in expandable
        sroots = first(_obs_expansion(cf, [s]; stop_at))
        overlap = intersect(sroots, self_deps)
        isempty(overlap) && continue
        push!(lines, " - :$s is an observable depending on $(_symrepr(sroots)), which \
                       introduces $(_symrepr(overlap))")
        # suggesting the formula's own kind fits either: an init formula must be pinned by an
        # init formula, a guess formula reads both pin-sets (see `pinned_obssyms`)
        kind = f isa InitFormula ? "an InitFormula" : "a GuessFormula"
        push!(fixes, isdisjoint(overlap, pinned) ?
            " - stop expanding :$s by declaring $kind for it, then it is read as-is" :
            " - $(_symrepr(overlap)) is declared by this formula itself, so no further formula \
               breaks the loop — read what $(_symrepr(overlap)) is computed from instead")
    end
    # no observable input explains the overlap, so it came in through canonicalization
    isempty(lines) && return [" - $(_symrepr(self_deps)) is read and written: members of one \
                               alias class are one variable"]
    append!(lines, fixes)
end

# `:a` for a single symbol, `[:a, :b]` for several — brackets around a lone symbol are noise
_symrepr(syms) = length(syms) == 1 ? ":$(only(syms))" : repr(syms)

"""
    UnresolvableExpansionError <: Exception

Thrown by a normalized formula when expanding its observable inputs yields NaN although
every root resolved. This is an *input* which could not be resolved, just detected one layer
deeper than the plain dict lookup — so the callers translate it according to their
`error_unresolvable` setting rather than letting it escape.
"""
struct UnresolvableExpansionError <: Exception
    msg::String
end
Base.showerror(io::IO, e::UnresolvableExpansionError) = print(io, e.msg)

# The roots were all resolvable (`apply_init_formulas!` checks that before calling us), so a
# NaN can only have come out of the expansion itself: an observable depending explicitly on
# time evaluated at `t=NaN`, or one leaving its domain (`log` of a negative root, …) at the
# values at hand.
function _assert_expansion_resolved(invals, origsym, roots, t)
    any(isnan, invals) || return nothing
    bad = origsym[findall(isnan, invals)]
    throw(UnresolvableExpansionError("expanding $bad to the settable root(s) $roots \
        produced NaN even though every root resolved. Either those observables depend \
        explicitly on time and no initialization time was given (t=$t), or they leave \
        their domain at the values at hand."))
end
