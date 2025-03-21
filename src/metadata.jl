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
            has_$($(QuoteNode(md)))(nw::Network, sni::SymbolicIndex)

        Checks if a `$($(QuoteNode(md)))` value is present for symbol `sym`.

        See also [`get_$($(QuoteNode(md)))`](@ref), [`set_$($(QuoteNode(md)))!`](@ref).
        """
        $fname_has(c::ComponentModel, sym::Symbol) = has_metadata(c, sym, $(QuoteNode(md)))
        $fname_has(nw::Network, sni::SymbolicIndex) = $fname_has(getcomp(nw, sni), sni.subidx)

        """
            get_$($(QuoteNode(md)))(c::ComponentModel, sym::Symbol)
            get_$($(QuoteNode(md)))(nw::Network, sni::SymbolicIndex)

        Returns the `$($(QuoteNode(md)))` value for symbol `sym`.

        See also [`has_$($(QuoteNode(md)))`](@ref), [`set_$($(QuoteNode(md)))!`](@ref).
        """
        $fname_get(c::ComponentModel, sym::Symbol) = get_metadata(c, sym, $(QuoteNode(md)))
        $fname_get(nw::Network, sni::SymbolicIndex) = $fname_get(getcomp(nw, sni), sni.subidx)


        """
            set_$($(QuoteNode(md)))!(c::ComponentModel, sym::Symbol, value)
            set_$($(QuoteNode(md)))!(nw::Network, sni::SymbolicIndex, value)

        Sets the `$($(QuoteNode(md)))` value for symbol `sym` to `value`.

        See also [`has_$($(QuoteNode(md)))`](@ref), [`get_$($(QuoteNode(md)))`](@ref).
        """
        $fname_set(c::ComponentModel, sym::Symbol, val) = set_metadata!(c, sym, $(QuoteNode(md)), val)
        $fname_set(nw::Network, sni::SymbolicIndex, val) = $fname_set(getcomp(nw, sni), sni.subidx, val)
    end
end
set_graphelement!(c::EdgeModel, p::Pair) = set_graphelement!(c, (;src=p.first, dst=p.second))

"""
    is_unused(c::ComponentModel, sym::Symbol)

Checks if symbol `sym` is marked as unused (i.e. it does not appear in the equations of f and g explicitly).
"""
function is_unused(c::ComponentModel, sym::Symbol)
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

Checks if a `default` or `init` value is present for symbol `sym`.
"""
has_default_or_init(c::ComponentModel, sym::Symbol) = has_default(c, sym) || has_init(c, sym)
has_default_or_init(nw::Network, sni::SymbolicIndex) = has_default_or_init(getcomp(nw, sni), sni.subidx)
"""
    get_default_or_init(c::ComponentModel, sym::Symbol)
    get_default_or_init(nw::Network, sni::SymbolicIndex)

Returns if a `default` value if available, otherwise returns `init` value for symbol `sym`.
"""
get_default_or_init(c::ComponentModel, sym::Symbol) = has_default(c, sym) ? get_default(c, sym) : get_init(c, sym)
get_default_or_init(nw::Network, sni::SymbolicIndex) = get_default_or_init(getcomp(nw, sni), sni.subidx)

#### default or guess
"""
    has_default_or_guess(c::ComponentModel, sym::Symbol)
    has_default_or_guess(nw::Network, sni::SymbolicIndex)

Checks if a `default` or `guess` value is present for symbol `sym`.
"""
has_default_or_guess(c::ComponentModel, sym::Symbol) = has_default(c, sym) || has_guess(c, sym)
has_default_or_guess(nw::Network, sni::SymbolicIndex) = has_default_or_guess(getcomp(nw, sni), sni.subidx)
"""
    get_default_or_guess(c::ComponentModel, sym::Symbol)
    get_default_or_guess(nw::Network, sni::SymbolicIndex)

Returns if a `default` value if available, otherwise returns `guess` value for symbol `sym`.
"""
get_default_or_guess(c::ComponentModel, sym::Symbol) = has_default(c, sym) ? get_default(c, sym) : get_guess(c, sym)
get_default_or_guess(nw::Network, sni::SymbolicIndex) = get_default_or_guess(getcomp(nw, sni), sni.subidx)


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

"""
    has_callback(c::ComponentModel)

Checks if the component has a callback function in metadata.
"""
has_callback(c::ComponentModel) = has_metadata(c, :callback)
"""
    get_callback(c::ComponentModel)

Gets all callback functions for the component. Wraps in tuple, even if there is only a single one.
"""
function get_callbacks(c::ComponentModel)
    cb = get_metadata(c, :callback)
    cb isa ComponentCallback ? (cb,) : cb
end
"""
    set_callback!(c::ComponentModel, cb)

Sets the callback function for the component. Overwrites any existing callback.
See also [`add_callback!`](@ref).
"""
function set_callback!(c::ComponentModel, cb; check=true)
    if !(cb isa ComponentCallback) && !(cb isa NTuple{N, <:ComponentCallback} where N)
        throw(ArgumentError("Callback must be a ComponentCallback or a tuple of ComponentCallbacks, got $(typeof(cb))."))
    end
    check && assert_cb_compat(c, cb)
    set_metadata!(c, :callback, cb)
end
"""
    add_callback!(c::ComponentModel, cb)

Adds a callback function to the component. Does not overwrite existing callbacks.
See also [`set_callback!`](@ref).
"""
function add_callback!(c::ComponentModel, cb; check=true)
    check && assert_cb_compat(c, cb)
    newcb = has_callback(c) ? (get_callbacks(c)..., cb) : (cb, )
    set_metadata!(c, :callback, newcb)
end

function get_defaults(c::ComponentModel, syms; missing_val=nothing)
    [has_default(c, sym) ? get_default(c, sym) : missing_val for sym in syms]
end
function get_guesses(c::ComponentModel, syms; missing_val=nothing)
    [has_guess(c, sym) ? get_guess(c, sym) : missing_val for sym in syms]
end
function get_defaults_or_inits(c::ComponentModel, syms; missing_val=nothing)
    [has_default_or_init(c, sym) ? get_default_or_init(c, sym) : missing_val for sym in syms]
end

# generate methods and docstrings for position and marker
for md in [:position, :marker]
    fname_has = Symbol(:has_, md)
    fname_get = Symbol(:get_, md)
    fname_set = Symbol(:set_, md, :!)
    @eval begin
        """
            has_$($(QuoteNode(md)))(v::VertexModel)

        Checks if vertex `v` has `$($(QuoteNode(md)))` metadata.

        See also: [`get_$($(QuoteNode(md)))`](@ref), [`set_$($(QuoteNode(md)))!`](@ref).
        """
        $fname_has(c::VertexModel) = has_metadata(c, $(QuoteNode(md)))

        """
            get_$($(QuoteNode(md)))(v::VertexModel)

        Returns the `$($(QuoteNode(md)))` metadata of vertex `v`. Might error if not present.

        See also: [`has_$($(QuoteNode(md)))`](@ref), [`set_$($(QuoteNode(md)))!`](@ref).
        """
        $fname_get(c::VertexModel) = get_metadata(c, $(QuoteNode(md)))

        """
            set_$($(QuoteNode(md)))!(v::VertexModel, val)

        Sets the `$($(QuoteNode(md)))` metadata of vertex `v` to `val`.

        See also: [`has_$($(QuoteNode(md)))`](@ref), [`get_$($(QuoteNode(md)))`](@ref).
        """
        $fname_set(c::VertexModel, val) = set_metadata!(c, $(QuoteNode(md)), val)
    end
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
    insymidx = _append_states!(lns, cf, sort(insym_all(cf)); sigdigits)
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
