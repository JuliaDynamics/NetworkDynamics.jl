####
#### per sym metadata
####
function _assert_symbol_exists(c::ComponentModel, s::Symbol)
    contains = s ∈ sym(c) ||
               s ∈ psym(c) ||
               (hasinsym(c) && s ∈ insym_all(c)) ||
               s ∈ outsym_flat(c) ||
               s ∈ obssym(c)
    contains || throw(ArgumentError("Symbol $s does not exist in component model."))
    return nothing
end

# to account for both _metadata(::ComponentModel, ::Symbol) and _metadata(::Network, ::SymbolicIndex)
const Comp_or_NW = Union{ComponentModel, Network}

"""
    has_metadata(c::ComponentModel, sym::Symbol, key::Symbol)
    has_metadata(nw::Network, sni::SymbolicIndex, key::Symbol)

Checks if symbol metadata `key` is present for symbol `sym` in a component model,
or for a symbol referenced by `sni` in a network.
Throws an error if the symbol does not exist in the component model.
"""
function has_metadata(c::ComponentModel, sym::Symbol, key::Symbol)
    _assert_symbol_exists(c, sym)
    md = symmetadata(c)
    haskey(md, sym) && haskey(md[sym], key)
end
function has_metadata(nw::Network, sym::SymbolicIndex, key::Symbol)
    has_metadata(getcomp(nw, sym), sym.subidx, key)
end

"""
    get_metadata(c::ComponentModel, sym::Symbol, key::Symbol)
    get_metadata(nw::Network, sni::SymbolicIndex, key::Symbol)

Retrieves the metadata `key` for symbol `sym` in a component model,
or for a symbol referenced by `sni` in a network.
Throws an error if the symbol does not exist in the component model.
"""
function get_metadata(c::ComponentModel, sym::Symbol, key::Symbol)
    _assert_symbol_exists(c, sym)
    symmetadata(c)[sym][key]
end
function get_metadata(nw::Network, sym::SymbolicIndex, key::Symbol)
    get_metadata(getcomp(nw, sym), sym.subidx, key)
end

"""
    set_metadata!(c::ComponentModel, sym::Symbol, key::Symbol, value)
    set_metadata!(nw::Network, sni::SymbolicIndex, key::Symbol, value)
    set_metadata!(c::ComponentModel, sym::Symbol, pair::Pair)
    set_metadata!(nw::Network, sni::SymbolicIndex, pair::Pair)

Sets the metadata `key` for symbol `sym` to `value` in a component model,
or for a symbol referenced by `sni` in a network.
Throws an error if the symbol does not exist in the component model.
"""
function set_metadata!(c::ComponentModel, sym::Symbol, key::Symbol, value)
    _assert_symbol_exists(c, sym)
    d = get!(symmetadata(c), sym, Dict{Symbol,Any}())
    d[key] = value
end
function set_metadata!(nw::Network, sym::SymbolicIndex, key::Symbol, value)
    set_metadata!(getcomp(nw, sym), sym.subidx, key, value)
end

function set_metadata!(c::ComponentModel, sym::Symbol, pair::Pair)
    _assert_symbol_exists(c, sym)
    set_metadata!(c, sym, pair.first, pair.second)
end
function set_metadata!(nw::Network, sym::SymbolicIndex, pair::Pair)
    set_metadata!(getcomp(nw, sym), sym.subidx, pair)
end

"""
    delete_metadata!(c::ComponentModel, sym::Symbol, key::Symbol)
    delete_metadata!(nw::Network, sni::SymbolicIndex, key::Symbol)

Removes the metadata `key` for symbol `sym` in a component model,
or for a symbol referenced by `sni` in a network.
Returns `true` if the metadata existed and was removed, `false` otherwise.
Throws an error if the symbol does not exist in the component model.
"""
function delete_metadata!(c::ComponentModel, sym::Symbol, key::Symbol)
    _assert_symbol_exists(c, sym)
    md = symmetadata(c)
    if haskey(md, sym) && haskey(md[sym], key)
        delete!(md[sym], key)
        if isempty(md[sym])
            delete!(md, sym)
        end
        return true
    end
    return false
end
delete_metadata!(nw::Network, sym::SymbolicIndex, key::Symbol) = delete_metadata!(getcomp(nw, sym), sym.subidx, key)

# generate default methods for some per-symbol metadata fields
for md in [:default, :guess, :init, :bounds]
    fname_has = Symbol(:has_, md)
    fname_get = Symbol(:get_, md)
    fname_set = Symbol(:set_, md, :!)
    fname_del = Symbol(:delete_, md, :!)
    @eval begin
        """
            has_$($(QuoteNode(md)))(c::ComponentModel, sym::Symbol)
            has_$($(QuoteNode(md)))(nw::Network, sni::SymbolicIndex)

        Checks if a `$($(QuoteNode(md)))` value is present for symbol `sym` in a component model,
        or for a symbol referenced by `sni` in a network.

        See also [`get_$($(QuoteNode(md)))`](@ref), [`set_$($(QuoteNode(md)))!`](@ref).
        """
        $fname_has(c::Comp_or_NW, sym) = has_metadata(c, sym, $(QuoteNode(md)))

        """
            get_$($(QuoteNode(md)))(c::ComponentModel, sym::Symbol)
            get_$($(QuoteNode(md)))(nw::Network, sni::SymbolicIndex)

        Returns the `$($(QuoteNode(md)))` value for symbol `sym` in a component model,
        or for a symbol referenced by `sni` in a network.

        See also [`has_$($(QuoteNode(md)))`](@ref), [`set_$($(QuoteNode(md)))!`](@ref).
        """
        $fname_get(c::Comp_or_NW, sym) = get_metadata(c, sym, $(QuoteNode(md)))

        """
            set_$($(QuoteNode(md)))!(c::ComponentModel, sym::Symbol, value)
            set_$($(QuoteNode(md)))!(nw::Network, sni::SymbolicIndex, value)

        Sets the `$($(QuoteNode(md)))` value for symbol `sym` to `value` in a component model,
        or for a symbol referenced by `sni` in a network.

        See also [`has_$($(QuoteNode(md)))`](@ref), [`get_$($(QuoteNode(md)))`](@ref).
        """
        $fname_set(c::Comp_or_NW, sym, val) = set_metadata!(c, sym, $(QuoteNode(md)), val)

        """
            delete_$($(QuoteNode(md)))!(c::ComponentModel, sym::Symbol)
            delete_$($(QuoteNode(md)))!(nw::Network, sni::SymbolicIndex)

        Removes the `$($(QuoteNode(md)))` value for symbol `sym` in a component model,
        or for a symbol referenced by `sni` in a network.

        See also [`has_$($(QuoteNode(md)))`](@ref), [`set_$($(QuoteNode(md)))!`](@ref).
        """
        $fname_del(c::Comp_or_NW, sym) = delete_metadata!(c, sym, $(QuoteNode(md)))
    end
end

"""
    is_unused(c::ComponentModel, sym::Symbol)
    is_unused(nw::Network, sni::SymbolicIndex)

Checks if symbol `sym` is marked as unused (i.e. it does not appear in the equations of f and g explicitly)
in a component model, or if a symbol referenced by `sni` is marked as unused in a network.
"""
function is_unused(c::Comp_or_NW, sym)
    if has_metadata(c, sym, :unused)
        get_metadata(c, sym, :unused)
    else
        false
    end
end

#### default or init
"""
    has_default_or_init(c::ComponentModel, sym::Symbol)
    has_default_or_init(nw::Network, sni::SymbolicIndex)

Checks if a `default` or `init` value is present for symbol `sym` in a component model,
or for a symbol referenced by `sni` in a network.
"""
has_default_or_init(c::Comp_or_NW, sym) = has_default(c, sym) || has_init(c, sym)

"""
    get_default_or_init(c::ComponentModel, sym::Symbol)
    get_default_or_init(nw::Network, sni::SymbolicIndex)

Returns if a `default` value if available, otherwise returns `init` value for symbol `sym` in a component model,
or for a symbol referenced by `sni` in a network.
"""
get_default_or_init(c::Comp_or_NW, sym) = has_default(c, sym) ? get_default(c, sym) : get_init(c, sym)

#### default or guess
"""
    has_default_or_guess(c::ComponentModel, sym::Symbol)
    has_default_or_guess(nw::Network, sni::SymbolicIndex)

Checks if a `default` or `guess` value is present for symbol `sym` in a component model,
or for a symbol referenced by `sni` in a network.
"""
has_default_or_guess(c::Comp_or_NW, sym) = has_default(c, sym) || has_guess(c, sym)

"""
    get_default_or_guess(c::ComponentModel, sym::Symbol)
    get_default_or_guess(nw::Network, sni::SymbolicIndex)

Returns if a `default` value if available, otherwise returns `guess` value for symbol `sym` in a component model,
or for a symbol referenced by `sni` in a network.
"""
get_default_or_guess(c::Comp_or_NW, sym) = has_default(c, sym) ? get_default(c, sym) : get_guess(c, sym)


"""
    get_defaults(c::ComponentModel, syms; missing_val=nothing)
    get_defaults(nw::Network, idx::Union{VIndex,EIndex}, syms; missing_val=nothing)

Gets all default values for the specified symbols in the component.
Returns `missing_val` for symbols without default values.
"""
function get_defaults(c::ComponentModel, syms; missing_val=nothing)
    [has_default(c, sym) ? get_default(c, sym) : missing_val for sym in syms]
end
get_defaults(nw::Network, idx::VCIndex, syms; kw...) = get_defaults(getcomp(nw, idx), syms; kw...)
get_defaults(nw::Network, idx::ECIndex, syms; kw...) = get_defaults(getcomp(nw, idx), syms; kw...)

"""
    get_guesses(c::ComponentModel, syms; missing_val=nothing)
    get_guesses(nw::Network, idx::Union{VIndex,EIndex}, syms; missing_val=nothing)

Gets all guess values for the specified symbols in the component.
Returns `missing_val` for symbols without guess values.
"""
function get_guesses(c::ComponentModel, syms; missing_val=nothing)
    [has_guess(c, sym) ? get_guess(c, sym) : missing_val for sym in syms]
end
get_guesses(nw::Network, idx::VCIndex, syms; kw...) = get_guesses(getcomp(nw, idx), syms; kw...)
get_guesses(nw::Network, idx::ECIndex, syms; kw...) = get_guesses(getcomp(nw, idx), syms; kw...)

"""
    get_defaults_or_inits(c::ComponentModel, syms; missing_val=nothing)
    get_defaults_or_inits(nw::Network, idx::Union{VIndex,EIndex}, syms; missing_val=nothing)

Gets all default or init values for the specified symbols in the component.
Returns `missing_val` for symbols without default or init values.
"""
function get_defaults_or_inits(c::ComponentModel, syms; missing_val=nothing)
    [has_default_or_init(c, sym) ? get_default_or_init(c, sym) : missing_val for sym in syms]
end
get_defaults_or_inits(nw::Network, idx::VCIndex, syms; kw...) = get_defaults_or_inits(getcomp(nw, idx), syms; kw...)
get_defaults_or_inits(nw::Network, idx::ECIndex, syms; kw...) = get_defaults_or_inits(getcomp(nw, idx), syms; kw...)

# extract varmaps
function _get_metadata_dict(c::ComponentModel, key; T=Any)
    dict = Dict{Symbol, T}()
    for (sym, smd) in c.symmetadata
        if haskey(smd, key)
            dict[sym] = smd[key]
        end
    end
    dict
end
"""
    get_defaults_dict(c::ComponentModel)

Returns a dictionary mapping symbols to their default values.
Only includes symbols that have default values set.

See also: [`get_guesses_dict`](@ref), [`get_inits_dict`](@ref)
"""
get_defaults_dict(c::ComponentModel) = _get_metadata_dict(c, :default, Float64)
"""
    get_guesses_dict(c::ComponentModel)

Returns a dictionary mapping symbols to their guess values.
Only includes symbols that have guess values set.

See also: [`get_defaults_dict`](@ref), [`get_inits_dict`](@ref)
"""
get_guesses_dict(c::ComponentModel) = _get_metadata_dict(c, :guess, Float64)
"""
    get_inits_dict(c::ComponentModel)

Returns a dictionary mapping symbols to their initialization values.
Only includes symbols that have initialization values set.

See also: [`get_defaults_dict`](@ref), [`get_guesses_dict`](@ref)
"""
get_inits_dict(c::ComponentModel) = _get_metadata_dict(c, :init, Float64)
"""
    get_bounds_dict(c::ComponentModel)

Returns a dictionary mapping symbols to their bounds values.
Only includes symbols that have bounds values set.

See also: [`get_defaults_dict`](@ref), [`get_guesses_dict`](@ref), [`get_inits_dict`](@ref)
"""
get_bounds_dict(c::ComponentModel) = _get_metadata_dict(c, :bounds, Tuple{Float64,Float64})
"""
    get_defaults_or_inits_dict(c::ComponentModel)

Returns a dictionary mapping symbols to their default values or initialization values (if no default exists).
Only includes symbols that have either default or init values set.

See also: [`get_defaults_dict`](@ref), [`get_guesses_dict`](@ref), [`get_inits_dict`](@ref)
"""
function get_defaults_or_inits_dict(c::ComponentModel)
    inits = get_inits_dict(c)
    defaults = get_defaults_dict(c)
    merge!(inits, defaults) # add defaults to inits, overwriting existing inits
end


"""
    set_defaults!(nw::Network, s::NWState)

Set the default values of the network to the values of the given state.
Can be used to "store" the found fixpoint in the network metadata.

Values of `missing`, `nothing` or `NaN` are ignored.
"""
function set_defaults!(nw::Network, s::NWState)
    for (sni, val) in zip(SII.variable_symbols(nw), uflat(s))
        val isa Number || continue
        isnan(val) && continue
        set_default!(nw, sni, val)
    end
    set_defaults!(nw, s.p)
    nw
end

"""
    set_defaults!(nw::Network, p::NWParameter)

Set the parameter default values of the network to the values of the given parameter object.

Values of `missing`, `nothing` or `NaN` are ignored.
"""
function set_defaults!(nw::Network, p::NWParameter)
    for (sni, val) in zip(SII.parameter_symbols(nw), pflat(p))
        val isa Number || continue
        isnan(val) && continue
        set_default!(nw, sni, val)
    end
    nw
end

####
#### Component metadata
####
# types for for referencing a single component
# not union to resolve dispatch ambiguity
const VCIndex = VIndex{<:Any, Nothing}
const ECIndex = EIndex{<:Any, Nothing}

"""
    has_metadata(c::ComponentModel, key::Symbol)
    has_metadata(nw::Network, idx::Union{VIndex,EIndex}, key::Symbol)

Checks if metadata `key` is present for the component.
"""
function has_metadata(c::ComponentModel, key)
    haskey(metadata(c), key)
end
has_metadata(nw::Network, idx::VCIndex, key::Symbol) = has_metadata(getcomp(nw, idx), key)
has_metadata(nw::Network, idx::ECIndex, key::Symbol) = has_metadata(getcomp(nw, idx), key)

"""
    get_metadata(c::ComponentModel, key::Symbol)
    get_metadata(nw::Network, idx::Union{VIndex,EIndex}, key::Symbol)

Retrieves the metadata `key` for the component.
"""
get_metadata(c::ComponentModel, key::Symbol) = metadata(c)[key]
get_metadata(nw::Network, idx::VCIndex, key::Symbol) = get_metadata(getcomp(nw, idx), key)
get_metadata(nw::Network, idx::ECIndex, key::Symbol) = get_metadata(getcomp(nw, idx), key)

"""
    set_metadata!(c::ComponentModel, key::Symbol, value)
    set_metadata!(nw::Network, idx::Union{VIndex,EIndex}, key::Symbol, value)

Sets the metadata `key` for the component to `value`.
"""
set_metadata!(c::ComponentModel, key::Symbol, val) = setindex!(metadata(c), val, key)
set_metadata!(nw::Network, idx::VCIndex, key::Symbol, val) = set_metadata!(getcomp(nw, idx), key, val)
set_metadata!(nw::Network, idx::ECIndex, key::Symbol, val) = set_metadata!(getcomp(nw, idx), key, val)

"""
    delete_metadata!(c::ComponentModel, key::Symbol)
    delete_metadata!(nw::Network, idx::Union{VIndex,EIndex}, key::Symbol)

Removes the component-wide metadata `key` from the component model,
or from a component referenced by `idx` in a network.
Returns `true` if the metadata existed and was removed, `false` otherwise.
"""
function delete_metadata!(c::ComponentModel, key::Symbol)
    md = metadata(c)
    if haskey(md, key)
        delete!(md, key)
        return true
    end
    return false
end
delete_metadata!(nw::Network, idx::VCIndex, key::Symbol) = delete_metadata!(getcomp(nw, idx), key)
delete_metadata!(nw::Network, idx::ECIndex, key::Symbol) = delete_metadata!(getcomp(nw, idx), key)


#### graphelement field for edges and vertices
"""
    has_graphelement(c::ComponentModel)
    has_graphelement(nw::Network, idx::Union{VIndex,EIndex})

Checks if the edge or vertex function has the `graphelement` metadata.
"""
has_graphelement(c::ComponentModel) = has_metadata(c, :graphelement)
has_graphelement(nw::Network, idx::VCIndex) = has_graphelement(getcomp(nw, idx))
has_graphelement(nw::Network, idx::ECIndex) = has_graphelement(getcomp(nw, idx))

"""
    get_graphelement(c::EdgeModel)::@NamedTuple{src::T, dst::T}
    get_graphelement(c::VertexModel)::Int
    get_graphelement(nw::Network, idx::Union{VIndex,EIndex})

Retrieves the `graphelement` metadata for the component model. For edges this
returns a named tuple `(;src, dst)` where both are either integers (vertex index)
or symbols (vertex name).
"""
get_graphelement(c::EdgeModel) = get_metadata(c, :graphelement)::@NamedTuple{src::T, dst::T} where {T<:Union{Int,Symbol}}
get_graphelement(c::VertexModel) = get_metadata(c, :graphelement)::Int
get_graphelement(nw::Network, idx::VCIndex) = get_graphelement(getcomp(nw, idx))
get_graphelement(nw::Network, idx::ECIndex) = get_graphelement(getcomp(nw, idx))

"""
    set_graphelement!(c::EdgeModel, nt::@NamedTuple{src::T, dst::T})
    set_graphelement!(c::EdgeModel, p::Pair)
    set_graphelement!(c::VertexModel, vidx::Int)
    set_graphelement!(nw::Network, idx::Union{VIndex,EIndex}, value)

Sets the `graphelement` metadata for the component. For edges this takes a
named tuple `(;src, dst)` where both are either integer (vertex index) or symbol
(vertex name). For vertices it takes a single integer `vidx`.
"""
set_graphelement!(c::EdgeModel, nt::@NamedTuple{src::T, dst::T}) where {T<:Union{Int,Symbol}} = set_metadata!(c, :graphelement, nt)
set_graphelement!(c::EdgeModel, p::Pair) = set_graphelement!(c, (;src=p.first, dst=p.second))
set_graphelement!(c::VertexModel, vidx::Int) = set_metadata!(c, :graphelement, vidx)
set_graphelement!(nw::Network, idx::VCIndex, value) = set_graphelement!(getcomp(nw, idx), value)
set_graphelement!(nw::Network, idx::ECIndex, value) = set_graphelement!(getcomp(nw, idx), value)

"""
    has_callback(c::ComponentModel)
    has_callback(nw::Network, idx::Union{VIndex,EIndex})

Checks if the component has a callback function in metadata.
"""
has_callback(c::ComponentModel) = has_metadata(c, :callback)
has_callback(nw::Network, idx::VCIndex) = has_callback(getcomp(nw, idx))
has_callback(nw::Network, idx::ECIndex) = has_callback(getcomp(nw, idx))

"""
    get_callbacks(c::ComponentModel)
    get_callbacks(nw::Network, idx::Union{VIndex,EIndex})

Gets all callback functions for the component. Wraps in tuple, even if there is only a single one.
"""
function get_callbacks(c::ComponentModel)
    cb = get_metadata(c, :callback)
    cb isa ComponentCallback ? (cb,) : cb
end
get_callbacks(nw::Network, idx::VCIndex) = get_callbacks(getcomp(nw, idx))
get_callbacks(nw::Network, idx::ECIndex) = get_callbacks(getcomp(nw, idx))

"""
    set_callback!(c::ComponentModel, cb; check=true)
    set_callback!(nw::Network, idx::Union{VIndex,EIndex}, cb; check=true)

Sets the callback function for the component. Overwrites any existing callback.
See also [`add_callback!`](@ref).
"""
function set_callback!(c::ComponentModel, cb; check=true)
    if !(cb isa ComponentCallback) && !(cb isa Tuple && all(c -> c isa ComponentCallback, cb))
        throw(ArgumentError("Callback must be a ComponentCallback or a tuple of ComponentCallbacks, got $(typeof(cb))."))
    end
    check && assert_cb_compat(c, cb)
    set_metadata!(c, :callback, cb)
end
set_callback!(nw::Network, idx::VCIndex, cb; kw...) = set_callback!(getcomp(nw, idx), cb; kw...)
set_callback!(nw::Network, idx::ECIndex, cb; kw...) = set_callback!(getcomp(nw, idx), cb; kw...)

"""
    add_callback!(c::ComponentModel, cb; check=true)
    add_callback!(nw::Network, idx::Union{VIndex,EIndex}, cb; check=true)

Adds a callback function to the component. Does not overwrite existing callbacks.
See also [`set_callback!`](@ref).
"""
function add_callback!(c::ComponentModel, cb; check=true)
    check && assert_cb_compat(c, cb)
    newcb = has_callback(c) ? (get_callbacks(c)..., cb) : (cb, )
    set_metadata!(c, :callback, newcb)
end
add_callback!(nw::Network, idx::VCIndex, cb; kw...) = add_callback!(getcomp(nw, idx), cb; kw...)
add_callback!(nw::Network, idx::ECIndex, cb; kw...) = add_callback!(getcomp(nw, idx), cb; kw...)

# generate methods and docstrings for position and marker
for md in [:position, :marker]
    fname_has = Symbol(:has_, md)
    fname_get = Symbol(:get_, md)
    fname_set = Symbol(:set_, md, :!)
    @eval begin
        """
            has_$($(QuoteNode(md)))(v::VertexModel)
            has_$($(QuoteNode(md)))(nw::Network, vidx::VIndex)

        Checks if vertex `v` has `$($(QuoteNode(md)))` metadata.

        See also: [`get_$($(QuoteNode(md)))`](@ref), [`set_$($(QuoteNode(md)))!`](@ref).
        """
        $fname_has(c::VertexModel) = has_metadata(c, $(QuoteNode(md)))
        $fname_has(nw::Network, idx::VCIndex) = $fname_has(getcomp(nw, idx))

        """
            get_$($(QuoteNode(md)))(v::VertexModel)
            get_$($(QuoteNode(md)))(nw::Network, vidx::VIndex)

        Returns the `$($(QuoteNode(md)))` metadata of vertex `v`. Might error if not present.

        See also: [`has_$($(QuoteNode(md)))`](@ref), [`set_$($(QuoteNode(md)))!`](@ref).
        """
        $fname_get(c::VertexModel) = get_metadata(c, $(QuoteNode(md)))
        $fname_get(nw::Network, idx::VCIndex) = $fname_get(getcomp(nw, idx))

        """
            set_$($(QuoteNode(md)))!(v::VertexModel, val)
            set_$($(QuoteNode(md)))!(nw::Network, vidx::VIndex, val)

        Sets the `$($(QuoteNode(md)))` metadata of vertex `v` to `val`.

        See also: [`has_$($(QuoteNode(md)))`](@ref), [`get_$($(QuoteNode(md)))`](@ref).
        """
        $fname_set(c::VertexModel, val) = set_metadata!(c, $(QuoteNode(md)), val)
        $fname_set(nw::Network, idx::VCIndex, val) = $fname_set(getcomp(nw, idx), val)
    end
end

####
#### Extract initial state from component
####
"""
    get_initial_state(c::ComponentModel, syms; missing_val=nothing)
    get_initial_state(nw::Network, sni::SymbolicIndex; missing_val=nothing)

Returns the initial state for symbol `sym` (single symbol or vector) of the component model `c`.
Returns `missing_val` if the symbol is not initialized. Also works for observed symbols.

See also: [`dump_initial_state`](@ref).
"""
get_initial_state(c::ComponentModel, s::Symbol; kw...) = only(get_initial_state(c, (s,); kw...))
function get_initial_state(cf::ComponentModel, syms; missing_val=nothing)
    obsbuf = any(in(syms), obssym(cf)) ? _get_initial_observed(cf) : nothing
    map(syms) do sym
        if (idx = findfirst(isequal(sym), obssym(cf))) !== nothing
            obs = obsbuf[idx]
            isnan(obs) ? missing_val : obs
        elseif has_default_or_init(cf, sym)
            Float64(get_default_or_init(cf, sym))
        else
            missing_val
        end
    end
end
get_initial_state(nw::Network, sni; kw...) = get_initial_state(getcomp(nw, sni), sni.subidx; kw...)

function _get_initial_observed(cf)
    missing_val = NaN
    obs = Vector{Float64}(undef, length(obssym(cf)))
    u = get_defaults_or_inits(cf, sym(cf); missing_val)
    ins = if cf isa EdgeModel
        (get_defaults_or_inits(cf, insym(cf).src; missing_val),
         get_defaults_or_inits(cf, insym(cf).dst; missing_val))
    else
        (get_defaults_or_inits(cf, insym(cf); missing_val), )
    end
    p = get_defaults_or_inits(cf, psym(cf); missing_val)
    obsf(cf)(obs, u, ins..., p, NaN)
    obs
end

"""
    dump_initial_state([IO=stdout], cf::ComponentModel; sigdigits=5, p=true, obs=true)

Prints the initial state of the component model `cf` to `IO` (defaults to stdout). Optionally
contains parameters and observed.

See also: [`get_initial_state`](@ref) and [`dump_state`](@ref).
"""
dump_initial_state(cf::ComponentModel; kwargs...) = dump_initial_state(stdout, cf; kwargs...)
function dump_initial_state(io, cf::ComponentModel; sigdigits=5, p=true, obs=true)
    lns = AnnotatedString[]
    symidx  = _append_states!(lns, cf, sort(sym(cf)); sigdigits)
    psymidx = _append_states!(lns, cf, sort(psym(cf)); sigdigits)
    insymidx = _append_states!(lns, cf, sort(insym_flat(cf)); sigdigits)
    outsymidx = _append_states!(lns, cf, sort(outsym_flat(cf)); sigdigits)

    obsidx = _append_observed!(lns, cf; sigdigits)
    aligned = align_strings(lns)

    printstyled(io, "Inputs:\n", bold=true)
    _printlines(io, aligned, insymidx)
    printstyled(io, "States:\n", bold=true)
    _printlines(io, aligned, symidx)
    printstyled(io, "Outputs:\n", bold=true)
    _printlines(io, aligned, outsymidx)
    if p
        printstyled(io, "Parameters:\n", bold=true)
        _printlines(io, aligned, psymidx)
    end
    if obs
        if length(obsidx) == 0
            printstyled(io, "$(length(obssym(cf))) Observed symbols uninitialized.", bold=true)
        elseif length(obsidx) == length(obssym(cf))
            printstyled(io, "Observed:\n", bold=true)
        else
            diff = length(obssym(cf)) - length(obsidx)
            printstyled(io, "Observed ($diff additional uninitialized):\n", bold=true)
        end
        _printlines(io, aligned, obsidx; newline=false)
    end
end
function _append_states!(lns, cf, syms; sigdigits)
    fidx = length(lns)+1
    for sym in syms
        str = "  &"
        if is_unused(cf, sym)
            str *= styled"{gray:$(string(sym)) (unused)}"
        else
            str *= string(sym)
        end
        str *= " &&= "
        if has_default_or_init(cf, sym)
            val = get_default_or_init(cf, sym)
            val_str = str_significant(val; sigdigits, phantom_minus=true)
            if has_default(cf, sym)
                str*= styled"{blue:$(val_str)}"
            else
                str *= styled"{yellow:$(val_str)}"
            end
        else
            val = nothing
            str *= styled"{red: uninitialized}"
        end
        str *= "&&"
        if has_guess(cf, sym)
            guess = str_significant(get_guess(cf, sym); sigdigits, phantom_minus=true)
            str *= " (guess $guess)"
        end
        str *= "&&"
        if has_bounds(cf, sym)
            lb, ub = get_bounds(cf, sym)
            if isnothing(val) || bounds_satisfied(val, (lb, ub))
                str *= " (bounds $lb..$ub)"
            else
                str *= styled" {red:(bounds $lb..$ub not satisfied)}"
            end
        end
        push!(lns, str)
    end
    fidx:length(lns)
end
function _append_observed!(lns, cf; sigdigits)
    fidx = length(lns)+1
    syms = obssym(cf)
    obs = _get_initial_observed(cf)
    perm = sortperm(syms)
    for (sym, val) in zip(syms[perm], obs[perm])
        isnan(val) && continue
        str = "  &" * string(sym) * " &&= " * str_significant(val; sigdigits, phantom_minus=true)
        str *= "&& &&"
        if has_bounds(cf, sym)
            lb, ub = get_bounds(cf, sym)
            if bounds_satisfied(val, (lb, ub))
                str *= " (bounds $lb..$ub)"
            else
                str *= styled" {red:(bounds $lb..$ub not satisfied)}"
            end
        end
        push!(lns, str)
    end
    fidx:length(lns)
end
function _printlines(io, aligned, range; newline=true)
    lines = @views aligned[range]
    for i in eachindex(lines)
        if newline == false && i == lastindex(lines)
            print(io, lines[i])
        else
            println(io, lines[i])
        end
    end
end

"""
    dump_state([IO=stdout], sol, t, idx; sigdigits=5)

Takes a Network solution `sol` and prints the state at `t` as well as the initial
state of the specified component model to `IO` (defaults to `stdout`).

`idx` musst a valid component index, i.e. `VIndex` or `EIndex` without symbol specification.

    dump_state(sol, 1.0, VIndex(4))
    dump_state(sol, 1.0, EIndex(2))

See also: [`dump_initial_state`](@ref).
"""
dump_state(sol, t, idx; kwargs...) = dump_state(stdout, sol, t, idx; kwargs...)
function dump_state(io, sol, t, idx; sigdigits=3)
    cf = extract_nw(sol)[idx]
    _sym = sort(sym(cf))
    _outsym = sort(outsym_flat(cf))
    _insym = sort(insym_flat(cf))
    _psym = sort(psym(cf))
    _obssym = sort(obssym(cf))

    groups = [("Outputs", _outsym), ("States", _sym), ("Inputs", _insym), ("Parameters", _psym), ("Observables", _obssym)]
    filter!(g -> !isempty(g[2]), groups)
    allsym = reduce(vcat, g[2] for g in groups)
    nwidxs = idxtype(idx).(idx.compidx, allsym)
    u0s, uts = sol([sol.t[begin], t]; idxs=nwidxs)
    guesses = [has_guess(cf, sym) ? get_guess(cf, sym) : nothing for sym in allsym]
    inits = [has_init(cf, sym) ? get_init(cf, sym) : nothing for sym in allsym]
    defaults = [has_default(cf, sym) ? get_default(cf, sym) : nothing for sym in allsym]

    t_str = str_significant(t; sigdigits)
    t0_str = str_significant(sol.t[begin]; sigdigits)
    strings = AnnotatedString[styled"{bold:& &&  t=$t_str &&  t=$t0_str &&  Default/Init}"]
    for (sym, g, d, i, u0, ut) in zip(allsym, guesses, defaults, inits, u0s, uts)
        buf = AnnotatedIOBuffer()
        print(buf, " &"*string(sym)*" && ")

        if !isnothing(ut)
            print(buf, str_significant(ut; sigdigits, phantom_minus=true))
        end
        print(buf, " && ")

        if !isnothing(u0)
            print(buf, str_significant(u0; sigdigits, phantom_minus=true))
        end
        print(buf, " && ")

        if !isnothing(d)
            print(buf, str_significant(d; sigdigits, phantom_minus=true))
        elseif !isnothing(i)
            ival = str_significant(i; sigdigits, phantom_minus=true)
            print(buf, styled"{blue: $ival}")
            if !isnothing(g)
                gval = str_significant(i; sigdigits, phantom_minus=true)
                print(buf, " (from $gval)")
            end
        end

        s = read(seekstart(buf), AnnotatedString)
        push!(strings, s)
    end
    aligned = align_strings(strings)

    let i = 2
        println(io, aligned[1])
        for (label, syms) in groups
            println(io, styled"{bold:$label}")
            for sym in syms
                println(io, aligned[i])
                i += 1
            end
        end
    end
end


####
#### Metadata Accessors through Network
####
function aliased_changed(nw::Network; warn=true)
    vchanged = _has_changed_hash(nw.im.aliased_vertexms)
    echanged = _has_changed_hash(nw.im.aliased_edgems)
    changed = vchanged || echanged
    if changed && warn
        s = if vchanged && echanged
            "vertices and edges"
        elseif vchanged
            "vertices"
        else
            "edges"
        end
        @warn """
        The metadata of at least one of your aliased $s changed! Proceed with caution!

        Some edgem/vertexm provided to to the `Network` constructor alias eachother.
        Which means, the Network object references the same component model in
        multiple places. Thus, metadata changes (such as changing of default values or
        component initialization) will be reflected in multiple components. To prevent
        this use the `dealias=true` keyword or manualy `copy` edge/vertex models
        before creating the network.
        """
    end
    changed
end
function _has_changed_hash(aliased_cfs)
    isempty(aliased_cfs) && return false
    changed = false
    for (k, v) in aliased_cfs
        if hash(k) != v.hash
            changed = true
            break
        end
    end
    changed
end


function describe_vertices end
function describe_edges end
