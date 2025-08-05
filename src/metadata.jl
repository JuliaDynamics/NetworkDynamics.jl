# to account for both _metadata(::ComponentModel, ::Symbol) and _metadata(::Network, ::SymbolicIndex)
const Comp_or_NW = Union{ComponentModel, Network}
# types for for referencing a single component
# not union to resolve dispatch ambiguity
const VCIndex = VIndex{<:Any, Nothing}
const ECIndex = EIndex{<:Any, Nothing}

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

function _match_symbol_name(c::ComponentModel, pattern)
    allsym = Symbol[]
    append!(allsym, sym(c))
    append!(allsym, psym(c))
    hasinsym(c) && append!(allsym, insym_all(c))
    append!(allsym, outsym_flat(c))
    append!(allsym, obssym(c))
    allssym = unique!(allsym)

    fidx = findall(sym -> contains(string(sym), pattern), allsym)
    if isempty(fidx)
        throw(ArgumentError("No symbol matching pattern $pattern found in component model."))
    elseif length(fidx) > 1
        throw(ArgumentError("Multiple symbols matching pattern $pattern found in component model: $(allsym[fidx]). Please use a more specific pattern."))
    end
    return Symbol(allsym[only(fidx)])
end

"""
    has_metadata(c::ComponentModel, sym::Symbol, key::Symbol)
    has_metadata(nw::Network, sni::SymbolicIndex, key::Symbol)

Checks if symbol metadata `key` is present for symbol `sym` in a component model,
or for a symbol referenced by `sni` in a network.

`sym` can also be a String or Regex, to address the only symbol containing the pattern, see [`set_metadata!`](@ref) for details.

Throws an error if the symbol does not exist in the component model.
"""
function has_metadata(c::ComponentModel, sym::Symbol, key::Symbol)
    _assert_symbol_exists(c, sym)
    md = symmetadata(c)
    haskey(md, sym) && haskey(md[sym], key)
end
function has_metadata(c::ComponentModel, pattern::Union{String,Regex}, key::Symbol)
    sym = _match_symbol_name(c, pattern)
    has_metadata(c, sym, key)
end
function has_metadata(nw::Network, sym::SymbolicIndex, key::Symbol)
    has_metadata(getcomp(nw, sym), sym.subidx, key)
end

"""
    get_metadata(c::ComponentModel, sym::Symbol, key::Symbol)
    get_metadata(nw::Network, sni::SymbolicIndex, key::Symbol)

Retrieves the metadata `key` for symbol `sym` in a component model,
or for a symbol referenced by `sni` in a network.

`sym` can also be a String or Regex, to address the only symbol containing the pattern, see [`set_metadata!`](@ref) for details.

Throws an error if the symbol does not exist in the component model.
"""
function get_metadata(c::ComponentModel, sym::Symbol, key::Symbol)
    _assert_symbol_exists(c, sym)
    symmetadata(c)[sym][key]
end
function get_metadata(c::ComponentModel, pattern::Union{String,Regex}, key::Symbol)
    sym = _match_symbol_name(c, pattern)
    get_metadata(c, sym, key)
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

For component models, you can also use a `String` or `Regex` pattern to match symbol names:
- String patterns use substring matching (e.g., `"δ"` matches `machine₊δ`)
- Regex patterns use full regex matching (e.g., `r"P\$"` matches symbols ending with "P")
This will error if there are none or multiple matches.

If the pattern matches multiple symbols, an error is thrown. Use a more specific pattern.
Throws an error if the symbol does not exist in the component model.
"""
function set_metadata!(c::ComponentModel, sym::Symbol, key::Symbol, value)
    _assert_symbol_exists(c, sym)
    d = get!(symmetadata(c), sym, Dict{Symbol,Any}())
    d[key] = value
end
function set_metadata!(c::ComponentModel, pattern::Union{String,Regex}, key::Symbol, value)
    sym = _match_symbol_name(c, pattern)
    set_metadata!(c, sym, key, value)
end
function set_metadata!(nw::Network, sym::SymbolicIndex, key::Symbol, value)
    set_metadata!(getcomp(nw, sym), sym.subidx, key, value)
end

function set_metadata!(c::ComponentModel, sym::Symbol, pair::Pair)
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

`sym` can also be a String or Regex, to address the only symbol containing the pattern, see [`set_metadata!`](@ref) for details.

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
function delete_metadata!(c::ComponentModel, pattern::Union{String,Regex}, key::Symbol)
    sym = _match_symbol_name(c, pattern)
    delete_metadata!(c, sym, key)
end
delete_metadata!(nw::Network, sym::SymbolicIndex, key::Symbol) = delete_metadata!(getcomp(nw, sym), sym.subidx, key)

"""
    strip_metadata!(c::ComponentModel, key::Symbol)

Remove all metadata of type `key` from the component.
"""
function strip_metadata!(c::ComponentModel, key::Symbol)
    md_dict = symmetadata(c)
    for (sym, sym_md) in md_dict
        delete!(sym_md, key)
    end
    return c
end
strip_metadata!(nw::Network, sym::SymbolicIndex, key::Symbol) = strip_metadata!(getcomp(nw, sym), key)

# generate default methods for some per-symbol metadata fields
for md in [:default, :guess, :init, :bounds]
    fname_has = Symbol(:has_, md)
    fname_get = Symbol(:get_, md)
    fname_set = Symbol(:set_, md, :!)
    fname_del = Symbol(:delete_, md, :!)
    fname_strip = if md == :guess
        :strip_guesses!
    elseif md == :bounds
        :strip_bounds!
    else
        Symbol(:strip_, string(md)*"s", :!)
    end
    @eval begin
        """
            has_$($(QuoteNode(md)))(c::ComponentModel, sym::Symbol)
            has_$($(QuoteNode(md)))(nw::Network, sni::SymbolicIndex)

        Checks if a `$($(QuoteNode(md)))` value is present for symbol `sym` in a component model,
        or for a symbol referenced by `sni` in a network.

        `sym` can be a String or Regex, to address the only symbol containing the pattern, see [`set_metadata!`](@ref) for details.

        See also [`get_$($(QuoteNode(md)))`](@ref), [`set_$($(QuoteNode(md)))!`](@ref).
        """
        $fname_has(c::Comp_or_NW, sym) = has_metadata(c, sym, $(QuoteNode(md)))

        """
            get_$($(QuoteNode(md)))(c::ComponentModel, sym::Symbol)
            get_$($(QuoteNode(md)))(nw::Network, sni::SymbolicIndex)

        Returns the `$($(QuoteNode(md)))` value for symbol `sym` in a component model,
        or for a symbol referenced by `sni` in a network.

        `sym` can be a String or Regex, to address the only symbol containing the pattern, see [`set_metadata!`](@ref) for details.

        See also [`has_$($(QuoteNode(md)))`](@ref), [`set_$($(QuoteNode(md)))!`](@ref).
        """
        $fname_get(c::Comp_or_NW, sym) = get_metadata(c, sym, $(QuoteNode(md)))

        """
            set_$($(QuoteNode(md)))!(c::ComponentModel, sym::Symbol, value)
            set_$($(QuoteNode(md)))!(nw::Network, sni::SymbolicIndex, value)

        Sets the `$($(QuoteNode(md)))` value for symbol `sym` to `value` in a component model,
        or for a symbol referenced by `sni` in a network.

        If `value` is `nothing` or `missing`, the metadata is removed.

        `sym` can be a String or Regex, to address the only symbol containing the pattern, see [`set_metadata!`](@ref) for details.

        See also [`has_$($(QuoteNode(md)))`](@ref), [`get_$($(QuoteNode(md)))`](@ref).
        """
        function $fname_set(c::Comp_or_NW, sym, val)
            if !isnothing(val) && !ismissing(val)
                set_metadata!(c, sym, $(QuoteNode(md)), val)
            else
                delete_metadata!(c, sym, $(QuoteNode(md)))
            end
        end

        """
            delete_$($(QuoteNode(md)))!(c::ComponentModel, sym::Symbol)
            delete_$($(QuoteNode(md)))!(nw::Network, sni::SymbolicIndex)

        Removes the `$($(QuoteNode(md)))` value for symbol `sym` in a component model,
        or for a symbol referenced by `sni` in a network.

        `sym` can be a String or Regex, to address the only symbol containing the pattern, see [`set_metadata!`](@ref) for details.

        See also [`has_$($(QuoteNode(md)))`](@ref), [`set_$($(QuoteNode(md)))!`](@ref).
        """
        $fname_del(c::Comp_or_NW, sym) = delete_metadata!(c, sym, $(QuoteNode(md)))

        """
            strip_$($(QuoteNode(md)))!(c::ComponentModel)
            strip_$($(QuoteNode(md)))!(nw::Network, idx::Union{VIndex,EIndex})

        Removes all `$($(QuoteNode(md)))` values from a component model,
        or from a component referenced by `idx` in a network.

        See also [`delete_$($(QuoteNode(md)))!`](@ref), [`set_$($(QuoteNode(md)))!`](@ref).
        """
        $fname_strip(c::Comp_or_NW) = strip_metadata!(c, $(QuoteNode(md)))
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
function _get_metadata_dict(c::ComponentModel, key, T; conv=identity)
    dict = Dict{Symbol, T}()
    for (sym, smd) in c.symmetadata
        if haskey(smd, key)
            dict[sym] = conv(smd[key])
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
get_bounds_dict(c::ComponentModel) = _get_metadata_dict(c, :bounds, Tuple{Float64,Float64}, conv=Tuple)
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

####
#### Callbacks
####

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

####
#### Init constraints
####

"""
    has_initconstraint(c::ComponentModel)
    has_initconstraint(nw::Network, idx::Union{VIndex,EIndex})

Checks if the component has an initialization constraint in metadata.

See also: [`get_initconstraints`](@ref), [`set_initconstraint!`](@ref).
"""
has_initconstraint(c::ComponentModel) = has_metadata(c, :initconstraint)
has_initconstraint(nw::Network, idx::VCIndex) = has_initconstraint(getcomp(nw, idx))
has_initconstraint(nw::Network, idx::ECIndex) = has_initconstraint(getcomp(nw, idx))

"""
    get_initconstraints(c::ComponentModel)
    get_initconstraints(nw::Network, idx::Union{VIndex,EIndex})

Gets all initialization constraints for the component. Returns a tuple, even if there is only a single constraint.

See also: [`has_initconstraint`](@ref), [`set_initconstraint!`](@ref), [`add_initconstraint!`](@ref).
"""
function get_initconstraints(c::ComponentModel)
    constraint = get_metadata(c, :initconstraint)
    constraint isa InitConstraint ? (constraint,) : constraint
end
get_initconstraints(nw::Network, idx::VCIndex) = get_initconstraints(getcomp(nw, idx))
get_initconstraints(nw::Network, idx::ECIndex) = get_initconstraints(getcomp(nw, idx))

"""
    set_initconstraint!(c::ComponentModel, constraint; check=true)
    set_initconstraint!(nw::Network, idx::Union{VIndex,EIndex}, constraint; check=true)

Sets the initialization constraint(s) for the component. Overwrites any existing constraints.
`constraint` can be a single `InitConstraint` or a tuple of `InitConstraint` objects.

See also: [`add_initconstraint!`](@ref), [`get_initconstraints`](@ref), [`delete_initconstraints!`](@ref).
"""
function set_initconstraint!(c::ComponentModel, constraint; check=true)
    if !(constraint isa InitConstraint) && !(constraint isa Tuple && all(c -> c isa InitConstraint, constraint))
        throw(ArgumentError("Constraint must be an InitConstraint or a tuple of InitConstraint objects, got $(typeof(constraint))."))
    end
    if check
        if constraint isa InitConstraint
            assert_initconstraint_compat(c, constraint)
        else
            for constr in constraint
                assert_initconstraint_compat(c, constr)
            end
        end
    end
    set_metadata!(c, :initconstraint, constraint)
end
set_initconstraint!(nw::Network, idx::VCIndex, constraint; kw...) = set_initconstraint!(getcomp(nw, idx), constraint; kw...)
set_initconstraint!(nw::Network, idx::ECIndex, constraint; kw...) = set_initconstraint!(getcomp(nw, idx), constraint; kw...)

"""
    add_initconstraint!(c::ComponentModel, constraint; check=true)
    add_initconstraint!(nw::Network, idx::Union{VIndex,EIndex}, constraint; check=true)

Adds an initialization constraint to the component. Does not overwrite existing constraints.
`constraint` should be a single `InitConstraint` object.

See also: [`set_initconstraint!`](@ref), [`get_initconstraints`](@ref).
"""
function add_initconstraint!(c::ComponentModel, constraint; check=true)
    if !(constraint isa InitConstraint)
        throw(ArgumentError("Constraint must be an InitConstraint, got $(typeof(constraint))."))
    end
    check && assert_initconstraint_compat(c, constraint)

    if has_initconstraint(c)
        existing_constraints = get_initconstraints(c)

        constraint ∈ existing_constraints && return false

        new_constraints = (existing_constraints..., constraint)
        set_metadata!(c, :initconstraint, new_constraints)
    else
        set_metadata!(c, :initconstraint, (constraint,))
    end
    return true
end
add_initconstraint!(nw::Network, idx::VCIndex, constraint; kw...) = add_initconstraint!(getcomp(nw, idx), constraint; kw...)
add_initconstraint!(nw::Network, idx::ECIndex, constraint; kw...) = add_initconstraint!(getcomp(nw, idx), constraint; kw...)

"""
    delete_initconstraints!(c::ComponentModel)
    delete_initconstraints!(nw::Network, idx::Union{VIndex,EIndex})

Removes the initialization constraint from the component model,
or from a component referenced by `idx` in a network.
Returns `true` if the constraint existed and was removed, `false` otherwise.

See also: [`set_initconstraint!`](@ref).
"""
delete_initconstraints!(c::ComponentModel) = delete_metadata!(c, :initconstraint)
delete_initconstraints!(nw::Network, idx::VCIndex) = delete_initconstraints!(getcomp(nw, idx))
delete_initconstraints!(nw::Network, idx::ECIndex) = delete_initconstraints!(getcomp(nw, idx))

####
#### Init formulas
####

"""
    has_initformula(c::ComponentModel)
    has_initformula(nw::Network, idx::Union{VIndex,EIndex})

Checks if the component has initialization formulas in metadata.

See also: [`get_initformulas`](@ref), [`set_initformula!`](@ref), [`add_initformula!`](@ref).
"""
has_initformula(c::ComponentModel) = has_metadata(c, :initformula)
has_initformula(nw::Network, idx::VCIndex) = has_initformula(getcomp(nw, idx))
has_initformula(nw::Network, idx::ECIndex) = has_initformula(getcomp(nw, idx))

"""
    get_initformulas(c::ComponentModel)
    get_initformulas(nw::Network, idx::Union{VIndex,EIndex})

Gets all initialization formulas for the component. Returns a tuple, even if there is only a single formula.

See also: [`has_initformula`](@ref), [`set_initformula!`](@ref), [`add_initformula!`](@ref).
"""
function get_initformulas(c::ComponentModel)
    formula = get_metadata(c, :initformula)
    formula isa InitFormula ? (formula,) : formula
end
get_initformulas(nw::Network, idx::VCIndex) = get_initformulas(getcomp(nw, idx))
get_initformulas(nw::Network, idx::ECIndex) = get_initformulas(getcomp(nw, idx))

"""
    set_initformula!(c::ComponentModel, formula; check=true)
    set_initformula!(nw::Network, idx::Union{VIndex,EIndex}, formula; check=true)

Sets the initialization formula(s) for the component. Overwrites any existing formulas.
`formula` can be a single `InitFormula` or a tuple of `InitFormula` objects.

See also: [`add_initformula!`](@ref), [`get_initformulas`](@ref), [`delete_initformulas!`](@ref).
"""
function set_initformula!(c::ComponentModel, formula; check=true)
    if !(formula isa InitFormula) && !(formula isa Tuple && all(f -> f isa InitFormula, formula))
        throw(ArgumentError("Formula must be an InitFormula or a tuple of InitFormula objects, got $(typeof(formula))."))
    end
    if check
        if formula isa InitFormula
            assert_initformula_compat(c, formula)
        else
            for f in formula
                assert_initformula_compat(c, f)
            end
        end
    end
    set_metadata!(c, :initformula, formula)
end
set_initformula!(nw::Network, idx::VCIndex, formula; kw...) = set_initformula!(getcomp(nw, idx), formula; kw...)
set_initformula!(nw::Network, idx::ECIndex, formula; kw...) = set_initformula!(getcomp(nw, idx), formula; kw...)

"""
    add_initformula!(c::ComponentModel, formula; check=true)
    add_initformula!(nw::Network, idx::Union{VIndex,EIndex}, formula; check=true)

Adds an initialization formula to the component. Does not overwrite existing formulas.
`formula` should be a single `InitFormula` object.

See also: [`set_initformula!`](@ref), [`get_initformulas`](@ref).
"""
function add_initformula!(c::ComponentModel, formula; check=true)
    if !(formula isa InitFormula)
        throw(ArgumentError("Formula must be an InitFormula, got $(typeof(formula))."))
    end
    check && assert_initformula_compat(c, formula)

    if has_initformula(c)
        existing_formulas = get_initformulas(c)

        formula ∈ existing_formulas && return false

        new_formulas = (existing_formulas..., formula)
        set_metadata!(c, :initformula, new_formulas)
    else
        set_metadata!(c, :initformula, (formula,))
    end
    return true
end
add_initformula!(nw::Network, idx::VCIndex, formula; kw...) = add_initformula!(getcomp(nw, idx), formula; kw...)
add_initformula!(nw::Network, idx::ECIndex, formula; kw...) = add_initformula!(getcomp(nw, idx), formula; kw...)

"""
    delete_initformulas!(c::ComponentModel)
    delete_initformulas!(nw::Network, idx::Union{VIndex,EIndex})

Removes all initialization formulas from the component model,
or from a component referenced by `idx` in a network.
Returns `true` if formulas existed and were removed, `false` otherwise.

See also: [`set_initformula!`](@ref), [`add_initformula!`](@ref).
"""
delete_initformulas!(c::ComponentModel) = delete_metadata!(c, :initformula)
delete_initformulas!(nw::Network, idx::VCIndex) = delete_initformulas!(getcomp(nw, idx))
delete_initformulas!(nw::Network, idx::ECIndex) = delete_initformulas!(getcomp(nw, idx))

# Deprecated backward compatibility
@deprecate delete_initconstraint!(args...; kwargs...) delete_initconstraints!(args...; kwargs...)
@deprecate delete_initformula!(args...; kwargs...) delete_initformulas!(args...; kwargs...)


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
    get_initial_state(c::ComponentModel, [state=get_defaults_or_inits_dict(c)], syms; missing_val=nothing)
    get_initial_state(nw::Network, sni::SymbolicIndex; missing_val=nothing)

Returns the initial state for symbol `sym` (single symbol or vector) of the component model `c`.
Returns `missing_val` if the symbol is not initialized. Also works for observed symbols.

See also: [`dump_initial_state`](@ref).
"""
get_initial_state(c::ComponentModel, state, s::Symbol; kw...) = only(get_initial_state(c, state, (s,); kw...))
get_initial_state(c::ComponentModel, s::Symbol; kw...) = only(get_initial_state(c, (s,); kw...))
function get_initial_state(cf::ComponentModel, syms; kw...)
    get_initial_state(cf, get_defaults_or_inits_dict(cf), syms; kw...)
end
function get_initial_state(cf::ComponentModel, state, syms; missing_val=nothing)
    obsbuf = if any(in(syms), obssym(cf))
        _get_component_observed(cf, state)
    else
        nothing
    end
    ret = Vector{Union{Float64, typeof(missing_val)}}(undef, length(syms))
    for (i, sym) in enumerate(syms)
        if (idx = findfirst(isequal(sym), obssym(cf))) !== nothing
            obs = obsbuf[idx]
            ret[i] = isnan(obs) ? missing_val : obs
        else haskey(state, sym)
            ret[i] = get(state, sym, missing_val)
        end
    end
    ret
end
get_initial_state(nw::Network, sni; kw...) = get_initial_state(getcomp(nw, sni), sni.subidx; kw...)

function _get_component_observed(cf, state=get_defaults_or_inits_dict(cf))
    obs = Vector{Float64}(undef, length(obssym(cf)))
    u = get.(Ref(state), sym(cf), NaN)
    ins = if cf isa EdgeModel
        (get.(Ref(state), insym(cf).src, NaN),
         get.(Ref(state), insym(cf).dst, NaN))
    else
        (get.(Ref(state), insym(cf), NaN), )
    end
    p = get.(Ref(state), psym(cf), NaN)
    obsf(cf)(obs, u, ins..., p, NaN)
    obs
end

"""
    dump_initial_state([IO=stdout], cf::ComponentModel,
                       [defaults=get_defaults_dict(cf)],
                       [inits=get_inits_dict(cf)],
                       [guesses=get_guesses_dict(cf)],
                       [bounds=get_bounds_dict(cf)];
                       sigdigits=5, p=true, obs=true)

Prints the initial state of the component model `cf` to `IO` (defaults to stdout). Optionally
contains parameters and observed.

See also: [`get_initial_state`](@ref) and [`dump_state`](@ref).
"""
dump_initial_state(cf::ComponentModel, args...; kwargs...) = dump_initial_state(stdout, cf, args...; kwargs...)
function dump_initial_state(io, cf::ComponentModel,
                            defaults=get_defaults_dict(cf),
                            inits=get_inits_dict(cf),
                            guesses=get_guesses_dict(cf),
                            bounds=get_bounds_dict(cf);
                            sigdigits=5, p=true, obs=true)
    # Create combined state dict for observed calculation
    dicts = (; defaults, inits, guesses, bounds)
    lns = AnnotatedString[]
    symidx  = _append_states!(lns, cf, sort(sym(cf)), dicts; sigdigits)
    psymidx = _append_states!(lns, cf, sort(psym(cf)), dicts; sigdigits)
    insymidx = _append_states!(lns, cf, sort(insym_flat(cf)), dicts; sigdigits)
    outsymidx = _append_states!(lns, cf, sort(outsym_flat(cf)), dicts; sigdigits)

    obsidx = _append_observed!(lns, cf, dicts; sigdigits)
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
function _append_states!(lns, cf, syms, dicts; sigdigits)
    fidx = length(lns)+1
    for sym in syms
        str = "  &"
        if is_unused(cf, sym)
            str *= styled"{gray:$(string(sym)) (unused)}"
        else
            str *= string(sym)
        end
        str *= " &&= "
        if haskey(dicts.defaults, sym)
            val = dicts.defaults[sym]
            val_str = str_significant(val; sigdigits, phantom_minus=true)
            str *= styled"{blue:$(val_str)}"
        elseif haskey(dicts.inits, sym)
            val = dicts.inits[sym]
            val_str = str_significant(val; sigdigits, phantom_minus=true)
            str *= styled"{yellow:$(val_str)}"
        else
            val = nothing
            str *= styled"{red: uninitialized}"
        end
        str *= "&&"
        if haskey(dicts.guesses, sym)
            guess = str_significant(dicts.guesses[sym]; sigdigits, phantom_minus=true)
            str *= " (guess $guess)"
        end
        str *= "&&"
        if haskey(dicts.bounds, sym)
            lb, ub = dicts.bounds[sym]
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
function _append_observed!(lns, cf, dicts; sigdigits)
    fidx = length(lns)+1
    syms = obssym(cf)
    obs = _get_component_observed(cf, merge(dicts.inits, dicts.defaults))
    perm = sortperm(syms)
    for (sym, val) in zip(syms[perm], obs[perm])
        isnan(val) && continue
        str = "  &" * string(sym) * " &&= " * str_significant(val; sigdigits, phantom_minus=true)
        str *= "&& &&"
        if haskey(dicts.bounds, sym)
            lb, ub = dicts.bounds[sym]
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

"""
    free_p(cf::ComponentModel)
    free_p(nw::Network)

Returns the "free" parameters (parameters without default values) for the given system.

# Returns
- Vector of parameter symbols that do not have default values set

See also: [`free_u`](@ref), [`has_default`](@ref), [`set_default!`](@ref)
"""
function free_p(cf::NetworkDynamics.ComponentModel)
    return filter(p -> !NetworkDynamics.has_default(cf, p), NetworkDynamics.psym(cf))
end
function free_p(nw::Network)
    setdiff(NetworkDynamics.SII.parameter_symbols(nw), keys(NetworkDynamics.SII.default_values(nw)))
end

"""
    free_u(nw::Network)
    free_u(cf::ComponentModel)

Returns the "free" variables/states (variables without default values) for the given system.

# Returns
- Vector of variable/state symbols that do not have default values set

See also: [`free_p`](@ref), [`has_default`](@ref), [`set_default!`](@ref)
"""
function free_u(nw::Network)
    setdiff(NetworkDynamics.SII.variable_symbols(nw), keys(NetworkDynamics.SII.default_values(nw)))
end
function free_u(cf::NetworkDynamics.ComponentModel)
    return filter(u -> !NetworkDynamics.has_default(cf, u), NetworkDynamics.sym(cf))
end
