####
#### per sym metadata
####
"""
    has_metadata(c::ComponentModel, sym::Symbol, key::Symbol)

Checks if symbol metadata `key` is present for symbol `sym`.
"""
function has_metadata(c::ComponentModel, sym::Symbol, key::Symbol)
    md = symmetadata(c)
    haskey(md, sym) && haskey(md[sym], key)
end
"""
    get_metadata(c::ComponentModel, sym::Symbol, key::Symbol)

Retrievs the metadata `key` for symbol `sym`.
"""
get_metadata(c::ComponentModel, sym::Symbol, key::Symbol) = symmetadata(c)[sym][key]

"""
    set_metadata!(c::ComponentModel, sym::Symbol, key::Symbol, value)
    set_metadata!(c::ComponentModel, sym::Symbol, pair)

Sets the metadata `key` for symbol `sym` to `value`.
"""
function set_metadata!(c::ComponentModel, sym::Symbol, key::Symbol, value)
    d = get!(symmetadata(c), sym, Dict{Symbol,Any}())
    d[key] = value
end
set_metadata!(c::ComponentModel, sym::Symbol, pair::Pair) = set_metadata!(c, sym, pair.first, pair.second)

# generate default methods for some per-symbol metadata fields
for md in [:default, :guess, :init, :bounds]
    fname_has = Symbol(:has_, md)
    fname_get = Symbol(:get_, md)
    fname_set = Symbol(:set_, md, :!)
    @eval begin
        """
            has_$($(QuoteNode(md)))(c::ComponentModel, sym::Symbol)

        Checks if a `$($(QuoteNode(md)))` value is present for symbol `sym`.

        See also [`get_$($(QuoteNode(md)))`](@ref), [`set_$($(QuoteNode(md)))!`](@ref).
        """
        $fname_has(c::ComponentModel, sym::Symbol) = has_metadata(c, sym, $(QuoteNode(md)))

        """
            get_$($(QuoteNode(md)))(c::ComponentModel, sym::Symbol)

        Returns the `$($(QuoteNode(md)))` value for symbol `sym`.

        See also [`has_$($(QuoteNode(md)))`](@ref), [`set_$($(QuoteNode(md)))!`](@ref).
        """
        $fname_get(c::ComponentModel, sym::Symbol) = get_metadata(c, sym, $(QuoteNode(md)))


        """
            set_$($(QuoteNode(md)))(c::ComponentModel, sym::Symbol, value)

        Sets the `$($(QuoteNode(md)))` value for symbol `sym` to `value`.

        See also [`has_$($(QuoteNode(md)))`](@ref), [`get_$($(QuoteNode(md)))`](@ref).
        """
        $fname_set(c::ComponentModel, sym::Symbol, val) = set_metadata!(c, sym, $(QuoteNode(md)), val)
    end
end
set_graphelement!(c::EdgeModel, p::Pair) = set_graphelement!(c, (;src=p.first, dst=p.second))


#### default or init
"""
    has_default_or_init(c::ComponentModel, sym::Symbol)

Checks if a `default` or `init` value is present for symbol `sym`.
"""
has_default_or_init(c::ComponentModel, sym::Symbol) = has_default(c, sym) || has_init(c, sym)
"""
    get_default_or_init(c::ComponentModel, sym::Symbol)

Returns if a `default` value if available, otherwise returns `init` value for symbol `sym`.
"""
get_default_or_init(c::ComponentModel, sym::Symbol) = has_default(c, sym) ? get_default(c, sym) : get_init(c, sym)

#### default or guess
"""
    has_default_or_guess(c::ComponentModel, sym::Symbol)

Checks if a `default` or `guess` value is present for symbol `sym`.
"""
has_default_or_guess(c::ComponentModel, sym::Symbol) = has_default(c, sym) || has_guess(c, sym)
"""
    get_default_or_guess(c::ComponentModel, sym::Symbol)

Returns if a `default` value if available, otherwise returns `guess` value for symbol `sym`.
"""
get_default_or_guess(c::ComponentModel, sym::Symbol) = has_default(c, sym) ? get_default(c, sym) : get_guess(c, sym)


# TODO: legacy, only used within show methods
function def(c::ComponentModel)::Vector{Union{Nothing,Float64}}
    map(c.sym) do s
        has_default_or_init(c, s) ? get_default_or_init(c, s) : nothing
    end
end
function guess(c::ComponentModel)::Vector{Union{Nothing,Float64}}
    map(c.sym) do s
        has_guess(c, s) ? get_guess(c, s) : nothing
    end
end
function pdef(c::ComponentModel)::Vector{Union{Nothing,Float64}}
    map(c.psym) do s
        has_default_or_init(c, s) ? get_default_or_init(c, s) : nothing
    end
end
function pguess(c::ComponentModel)::Vector{Union{Nothing,Float64}}
    map(c.psym) do s
        has_guess(c, s) ? get_guess(c, s) : nothing
    end
end

####
#### Component metadata
####
"""
    has_metadata(c::ComponentModel, key::Symbol)

Checks if metadata `key` is present for the component.
"""
function has_metadata(c::ComponentModel, key)
    haskey(metadata(c), key)
end
"""
    get_metadata(c::ComponentModel, key::Symbol)

Retrieves the metadata `key` for the component.
"""
get_metadata(c::ComponentModel, key::Symbol) = metadata(c)[key]
"""
    set_metadata!(c::ComponentModel, key::Symbol, value)

Sets the metadata `key` for the component to `value`.
"""
set_metadata!(c::ComponentModel, key::Symbol, val) = setindex!(metadata(c), val, key)

#### graphelement field for edges and vertices
"""
    has_graphelement(c)

Checks if the edge or vetex function function has the `graphelement` metadata.
"""
has_graphelement(c::ComponentModel) = has_metadata(c, :graphelement)
"""
    get_graphelement(c::EdgeModel)::@NamedTuple{src::T, dst::T}
    get_graphelement(c::VertexModel)::Int

Retrieves the `graphelement` metadata for the component model. For edges this
returns a named tupe `(;src, dst)` where both are either integers (vertex index)
or symbols (vertex name).
"""
get_graphelement(c::EdgeModel) = get_metadata(c, :graphelement)::@NamedTuple{src::T, dst::T} where {T<:Union{Int,Symbol}}
get_graphelement(c::VertexModel) = get_metadata(c, :graphelement)::Int
"""
    set_graphelement!(c::EdgeModel, src, dst)
    set_graphelement!(c::VertexModel, vidx)

Sets the `graphelement` metadata for the edge model. For edges this takes two
arguments `src` and `dst` which are either integer (vertex index) or symbol
(vertex name). For vertices it takes a single integer `vidx`.
"""
set_graphelement!(c::EdgeModel, nt::@NamedTuple{src::T, dst::T}) where {T<:Union{Int,Symbol}} = set_metadata!(c, :graphelement, nt)
set_graphelement!(c::VertexModel, vidx::Int) = set_metadata!(c, :graphelement, vidx)


function get_defaults(c::ComponentModel, syms; missing_val=nothing)
    [has_default(c, sym) ? get_default(c, sym) : missing_val for sym in syms]
end
function get_guesses(c::ComponentModel, syms; missing_val=nothing)
    [has_guess(c, sym) ? get_guess(c, sym) : missing_val for sym in syms]
end
function get_defaults_or_inits(c::ComponentModel, syms; missing_val=nothing)
    [has_default_or_init(c, sym) ? get_default_or_init(c, sym) : missing_val for sym in syms]
end

####
#### Extract initial state from component
####
"""
    get_initial_state(c::ComponentModel, syms; missing_val=nothing)

Returns the initial state for symbol `sym` (single symbol of vector) of the component model `c`.
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
    cf.obsf(obs, u, ins..., p, NaN)
    obs
end

"""
    dump_initial_state(cf::ComponentModel; sigdigits=5, p=true, obs=true)

Prints the initial state of the component model `cf` to the console. Optionally
contains parameters and observed.

See also: [`get_initial_state`](@ref).
"""
function dump_initial_state(cf::ComponentModel; sigdigits=5, p=true, obs=true)
    lns = AnnotatedString[]
    symidx  = _append_states!(lns, cf, sort(sym(cf)); sigdigits)
    psymidx = _append_states!(lns, cf, sort(psym(cf)); sigdigits)
    insymidx = _append_states!(lns, cf, sort(insym_all(cf)); sigdigits)
    outsymidx = _append_states!(lns, cf, sort(outsym_flat(cf)); sigdigits)

    obsidx = _append_observed!(lns, cf; sigdigits)
    aligned = align_strings(lns)

    printstyled("Inputs:\n", bold=true)
    _printlines(aligned, insymidx)
    printstyled("States:\n", bold=true)
    _printlines(aligned, symidx)
    printstyled("Outputs:\n", bold=true)
    _printlines(aligned, outsymidx)
    if p
        printstyled("Parameters:\n", bold=true)
        _printlines(aligned, psymidx)
    end
    if obs
        if length(obsidx) == 0
            printstyled("$(length(obssym(cf))) Observed symbols uninitialized.", bold=true)
        elseif length(obsidx) == length(obssym(cf))
            printstyled("Observed:\n", bold=true)
        else
            diff = length(obssym(cf)) - length(obsidx)
            printstyled("Observed ($diff additional uninitialized):\n", bold=true)
        end
        _printlines(aligned, obsidx; newline=false)
    end
end
function _append_states!(lns, cf, syms; sigdigits)
    fidx = length(lns)+1
    for sym in syms
        str = "  &" * string(sym) * " &&= "
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
function _printlines(aligned, range; newline=true)
    lines = @views aligned[range]
    for i in eachindex(lines)
        if newline == false && i == lastindex(lines)
            print(lines[i])
        else
            println(lines[i])
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
