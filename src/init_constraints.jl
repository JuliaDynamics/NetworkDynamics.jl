"""
    struct InitConstraint{F}
    InitConstraint(f, sym, dim)

A representation of an additional constraint that is applied during the initialization phase of a component.
It contains a function `f` that defines the constraint, a vector of symbols `sym` that are involved in the constraint,
and the dimension `dim` of the constraint.

    InitConstraint([:x, :y], 2) do res, u
        res[1] = u[:x]^2 + u[:y]^2 - 1
    end

See also [`@initconstraint`](@ref) for a macro to create such constraints.
"""
struct InitConstraint{F}
    f::F
    sym::Vector{Symbol}
    dim::Int
    prettyprint::Union{Nothing,String}
end
InitConstraint(f, sym, dim) = InitConstraint(f, sym, dim, nothing)

"""
    InitConstraint(subconstraints::InitConstraint...)

Combine multiple `InitConstraint` objects into a single constraint.

The resulting constraint will have:
- Combined symbols from all subconstraints (deduplicated)
- Total dimension equal to sum of individual dimensions
- Function that evaluates all subconstraints sequentially

Limitation: All subconstraints need to use `Symbol` indices internally, i.e. u[:x] rather than u[2]!
"""
function InitConstraint(subconstraints::InitConstraint...)
    isempty(subconstraints) && throw(ArgumentError("At least one subconstraint must be provided"))
    if length(subconstraints) == 1
        return only(subconstraints)
    end

    all_syms = unique!(reduce(vcat, c.sym for c in subconstraints))
    total_dim = sum(dim(c) for c in subconstraints)
    f_range = []
    offset = 0
    for c in subconstraints
        tup = (f=c.f, range=(offset+1):(offset+dim(c)))
        offset += dim(c)
        push!(f_range, tup)
    end

    f_range_tup = Tuple(f_range)
    function combined_f(res, u)
        offset = 0
        unrolled_foreach(f_range_tup) do (cf, cr)
            res_view = view(res, cr)
            cf(res_view, u)
        end
        nothing
    end

    try
        su = SymbolicView(zeros(length(all_syms)), all_syms, false)
        res = zeros(total_dim)
        combined_f(res, su)
    catch e
        if e isa IllegalIntIndexingError
            throw(ArgumentError(
                "Cannot combine multiple init constraints, because at least \
                 on of them uses `u[::Int]` indexing internally. Use `u[::Symbol]` \
                 within the `InitConstraint` to fix this error."
            ))
        else
            rethrow(e)
        end
    end

    if any(isnothing, c.prettyprint for c in subconstraints)
        prettyprint = nothing
    else
        header = "InitConstraint($all_syms) do out, u"
        footer = "end"
        bodylines = mapreduce(vcat, subconstraints) do c
            full = c.prettyprint
            split(full, '\n')[2:end-1] # remove header and footer
        end
        outidx = 1
        for i in eachindex(bodylines)
            bodylines[i] = replace(bodylines[i], r"out\[.*?\]" => "out[$(outidx)]")
            outidx += 1
        end
        prettyprint = join([header, join(bodylines, "\n"), footer], "\n")
    end

    InitConstraint(combined_f, all_syms, total_dim, prettyprint)
end

dim(c::InitConstraint) = c.dim

(c::InitConstraint)(res, u) = c.f(res, SymbolicView(u, c.sym))

function Base.show(io::IO, ::MIME"text/plain", @nospecialize(c::InitConstraint))
    if c.prettyprint === nothing
        r = repr(c)
        s = replace(r, ", nothing)"=>")")
        print(io, s)
    else
        print(io, c.prettyprint)
    end
end

"""
   @initconstraint

Generate an [`InitConstraint`](@ref) from an expression using symbols.

    @initconstraint begin
        :x + :y
        :z^2
    end

is equal to

    InitConstraint([:x, :y, :z], 2) do out, u
        out[1] = u[:x] + u[:y]
        out[2] = u[:z]^2
    end
"""
macro initconstraint(ex)
    if ex isa QuoteNode || ex.head != :block
        ex = Base.remove_linenums!(Expr(:block, ex))
    end

    sym = Symbol[]
    out = gensym(:out)
    u = gensym(:u)
    body = Expr[]

    dim = 0
    for constraint in ex.args
        constraint isa Union{Expr, QuoteNode} || continue # skip line number nodes
        dim += 1
        wrapped = _wrap_symbols!(constraint, sym, u)
        push!(body, :($(esc(out))[$dim] = $wrapped))
    end
    unique!(sym)

    s = join(string.(body), "\n")
    s = replace(s, "(\$(Expr(:escape, Symbol(\"$(string(out))\"))))" => "    out")
    s = replace(s, "(\$(Expr(:escape, Symbol(\"$(string(u))\"))))" => "u")
    s = "InitConstraint($sym, $dim) do out, u\n" * s * "\nend"

    quote
        InitConstraint($sym, $dim, $s) do $(esc(out)), $(esc(u))
            $(body...)
            nothing
        end
    end
end
function _wrap_symbols!(ex, sym, u)
    postwalk(ex) do x
        if x isa QuoteNode && x.value isa Symbol
            push!(sym, x.value)
            :($(esc(u))[$x])
        else
            x
        end
    end
end

# for metadata check, just passes down the
function assert_initconstraint_compat(cf::ComponentModel, c::InitConstraint)
    allowed_symbols = Set(vcat(
        sym(cf),
        psym(cf),
        insym_flat(cf),
        outsym_flat(cf),
        obssym(cf)
    ))
    missmatch = setdiff(c.sym, allowed_symbols)
    if !isempty(missmatch)
        throw(ArgumentError("InitConstraint requires symbols that are not part of the component model: $missmatch"))
    end
    c
end

"""
prep_initconstrasint allocates all the buffers and so on in order to be called
within the NonlienarProblem during initialization.
"""
function prep_initiconstraint(cm::ComponentModel, c::InitConstraint, chunksize)
    obscache = if !isempty(c.sym ∩ obssym(cm))
        DiffCache(zeros(length(obssym(cm))), chunksize)
    else
        nothing
    end
    symcache = DiffCache(zeros(length(c.sym)), chunksize)
    symtup = Tuple(c.sym)
    obsf! = obsf(cm)

    mapping! = generate_init_input_mapping(cm, c)
    initf! = c.f

    (res, outbufs, ubuf, inbufs, pbuf, t) -> begin
        if !isnothing(obscache)
            obsbuf = PreallocationTools.get_tmp(obscache, res)
            obsf!(obsbuf, ubuf, inbufs..., pbuf, t)
            obsf!
        else
            obsbuf = nothing
        end

        symbuf = PreallocationTools.get_tmp(symcache, res)
        mapping!(symbuf, outbufs, ubuf, inbufs, pbuf, obsbuf)

        symv = SymbolicView(symbuf, symtup)
        initf!(res, symv)

        nothing
    end
end
prep_initiconstraint(cm::ComponentModel, c::Nothing, _) = (args...) -> nothing

"""
Returns a function, which collects all the symbols required for InitConstraint.
For each symbol, it checks if it should be copied from flat outputs, inputs, parameters, or state
buffers.
"""
function generate_init_input_mapping(cm::ComponentModel, c::InitConstraint)
    outmapping = NTuple{3,Int}[]
    umapping   = NTuple{2,Int}[]
    inmapping  = NTuple{3,Int}[]
    pmapping   = NTuple{2,Int}[]
    obsmapping = NTuple{2,Int}[]
    for (iidx, s) in enumerate(c.sym)
        if s in outsym_flat(cm)
            candidates = findfirst.(Ref(isequal(s)), outsym_normalized(cm))
            candidx = findfirst(!isnothing, candidates)
            push!(outmapping, (iidx, candidx, candidates[candidx]))
        elseif s in sym(cm)
            unidx = findfirst(isequal(s), sym(cm))
            push!(umapping, (iidx, unidx))
        elseif s in insym_flat(cm)
            candidates = findfirst.(Ref(isequal(s)), insym_normalized(cm))
            candidx = findfirst(!isnothing, candidates)
            push!(inmapping, (iidx, candidx, candidates[candidx]))
        elseif s in psym(cm)
            pnidx = findfirst(isequal(s), psym(cm))
            push!(pmapping, (iidx, pnidx))
        elseif s in obssym(cm)
            obsidx = findfirst(isequal(s), obssym(cm))
            push!(obsmapping, (iidx, obsidx))
        else
            error("Could not locate symbol $s needed for constraint.")
        end
    end

    (syms, outs, u, ins, p, obs) -> begin
        for (iidx, outnr, outidx) in outmapping
            syms[iidx] = outs[outnr][outidx]
        end
        for (iidx, unidx) in umapping
            syms[iidx] = u[unidx]
        end
        for (iidx, inr, inidx) in inmapping
            syms[iidx] = ins[inr][inidx]
        end
        for (iidx, pnidx) in pmapping
            syms[iidx] = p[pnidx]
        end
        for (iidx, obsidx) in obsmapping
            syms[iidx] = obs[obsidx]
        end
    end
end

####
#### InitFormula
####
"""
    InitFormula(f, outsym, sym)

A representation of initialization formulas that are applied during the initialization phase of a component.
InitFormulas act earlier in the initialization pipeline than InitConstraints - they essentially set additional
defaults rather than adding equations to the nonlinear system.

It contains a function `f` that defines the formulas, a vector of output symbols `outsym` that will be set by the formulas,
a vector of input symbols `sym` that are used in the formulas, and an optional pretty-print string.

    InitFormula([:Vset], [:u_r, :u_i]) do out, u
        out[:Vset] = sqrt(u[:u_r]^2 + u[:u_i]^2)
    end

See also [`@initformula`](@ref) for a macro to create such formulas.
"""
struct InitFormula{F}
    f::F
    outsym::Vector{Symbol}   # output symbols (from LHS of assignments)
    sym::Vector{Symbol}      # input symbols (from RHS of assignments)
    prettyprint::Union{Nothing,String}

    function InitFormula(f::F, outsym::Vector{Symbol}, sym::Vector{Symbol}, prettyprint::Union{Nothing,String}) where F
        # Check for self-dependencies (formula depending on its own output)
        self_deps = intersect(sym, outsym)
        if !isempty(self_deps)
            throw(ArgumentError("InitFormula cannot depend on its own output symbols: $self_deps"))
        end
        new{F}(f, outsym, sym, prettyprint)
    end
end
InitFormula(f, outsym, sym) = InitFormula(f, outsym, sym, nothing)

dim(c::InitFormula) = length(c.outsym)

(c::InitFormula)(out, u) = c.f(out, SymbolicView(u, c.sym))

function Base.show(io::IO, ::MIME"text/plain", @nospecialize(c::InitFormula))
    if c.prettyprint === nothing
        r = repr(c)
        s = replace(r, ", nothing)"=>")")
        print(io, s)
    else
        print(io, c.prettyprint)
    end
end

"""
    GuessFormula(f, outsym, sym)

A representation of guess formulas that improve initial guesses during component initialization.
GuessFormulas are applied after InitFormulas in the pipeline, reading from both defaults and
guesses to compute improved guess values for free variables.

Unlike InitFormulas which set defaults (reducing free variables), GuessFormulas only update
the `guesses` dictionary to improve solver convergence without changing the problem dimension.

It contains a function `f` that defines the formulas, a vector of output symbols `outsym` that
will be set by the formulas, a vector of input symbols `sym` that are used in the formulas,
and an optional pretty-print string.

# Input Lookup Behavior
When evaluating a GuessFormula, input symbols are looked up with this priority:
1. First check `defaults` dict (fixed/known values take precedence)
2. Then check `guesses` dict (free variable current guesses)
3. Error if symbol not found in either

# Examples

    GuessFormula([:V, :theta], [:u_r, :u_i]) do out, u
        out[:V] = sqrt(u[:u_r]^2 + u[:u_i]^2)
        out[:theta] = atan(u[:u_i], u[:u_r])
    end

See also [`@guessformula`](@ref) for a macro to create such formulas.
"""
struct GuessFormula{F}
    f::F
    outsym::Vector{Symbol}   # output symbols (from LHS of assignments)
    sym::Vector{Symbol}      # input symbols (from RHS of assignments)
    prettyprint::Union{Nothing,String}

    function GuessFormula(f::F, outsym::Vector{Symbol}, sym::Vector{Symbol}, prettyprint::Union{Nothing,String}) where F
        # Check for self-dependencies (formula depending on its own output)
        self_deps = intersect(sym, outsym)
        if !isempty(self_deps)
            throw(ArgumentError("GuessFormula cannot depend on its own output symbols: $self_deps"))
        end
        new{F}(f, outsym, sym, prettyprint)
    end
end
GuessFormula(f, outsym, sym) = GuessFormula(f, outsym, sym, nothing)

dim(c::GuessFormula) = length(c.outsym)

(c::GuessFormula)(out, u) = c.f(out, SymbolicView(u, c.sym))

function Base.show(io::IO, ::MIME"text/plain", @nospecialize(c::GuessFormula))
    if c.prettyprint === nothing
        r = repr(c)
        s = replace(r, ", nothing)"=>")")
        print(io, s)
    else
        print(io, c.prettyprint)
    end
end

"""
   @initformula

Generate an [`InitFormula`](@ref) from an expression using symbols.

    @initformula begin
        :Vset = sqrt(:u_r^2 + :u_i^2)
        :Pset = :u_r * :i_r + :u_i * :i_i
    end

is equal to

    InitFormula([:Vset, :Pset], [:u_r, :u_i, :i_r, :i_i]) do out, u
        out[:Vset] = sqrt(u[:u_r]^2 + u[:u_i]^2)
        out[:Pset] = u[:u_r] * u[:i_r] + u[:u_i] * u[:i_i]
    end
"""
macro initformula(ex)
    _formula_macro(InitFormula, ex)
end


"""
   @guessformula

Generate a [`GuessFormula`](@ref) from an expression using symbols.

This macro provides convenient syntax for creating guess formulas using assignment
expressions with quoted symbols. Each assignment computes a guess value for a free variable.

# Syntax

```julia
@guessformula begin
    :output1 = expression_with(:input_symbols)
    :output2 = other_expression(:more_inputs)
end
```

# Input Symbol Lookup
Input symbols (on RHS) are looked up with this priority:
1. First from `defaults` dict (fixed values)
2. Then from `guesses` dict (current guesses)
3. Error if not found in either

# Output Symbol Target
Output symbols (on LHS) must be:
- Valid component symbols (states, parameters, inputs, outputs)
- NOT observables (observables are computed, not guessed)

See also: [`GuessFormula`](@ref), [`@initformula`](@ref), [`initialize_component`](@ref)
"""
macro guessformula(ex)
    _formula_macro(GuessFormula, ex)
end

# combined macro code for InitFormula and GuessFormula
#=
In a nutshell: wrap every QuoteNote symbol either in u[:sym] or out[:sym]
- if it appears alone in the LHS of an assignment, it is an output symbol
- otherwise it is an input symbol
some thinks will break!
for example: set!(:out, :in) -> set!(u[:out], u[:in])
=#
function _formula_macro(type, ex)
    if ex isa QuoteNode || ex.head != :block
        ex = Base.remove_linenums!(Expr(:block, ex))
    end

    input_syms = Symbol[]    # RHS symbols
    output_syms = Symbol[]   # LHS symbols
    out_var = gensym(:out)
    u = gensym(:u)
    body = Expr[]

    for formula in ex.args
        formula isa Union{Expr, QuoteNode} || continue # skip line number nodes
        if formula isa Expr && formula.head == :(=)
            rhs = formula.args[2]
            # Wrap symbols in the RHS
            wrapped_rhs = _wrap_symbols!(rhs, input_syms, u)

            lhs = formula.args[1]
            # Extract the target symbol from the LHS (should be a QuoteNode)
            if lhs isa QuoteNode && lhs.value isa Symbol
                # :output = ... assigmend
                target_sym = lhs.value
                push!(output_syms, target_sym)
                push!(body, :($(esc(out_var))[$(QuoteNode(target_sym))] = $wrapped_rhs))
            else
                # "normal" assigmend
                push!(body, :($lhs = $wrapped_rhs))
            end
        else
            wrapped_expr = _wrap_symbols!(formula, input_syms, u)
            push!(body, :($wrapped_expr))  # standalone expression
        end
    end
    unique!(input_syms)

    # Generate pretty print string
    s = "    " * join(string.(body), "\n    ")
    s = replace(s, "(\$(Expr(:escape, Symbol(\"$(string(out_var))\"))))" => "out")
    s = replace(s, "(\$(Expr(:escape, Symbol(\"$(string(u))\"))))" => "u")
    s = "$(type)($output_syms, $input_syms) do out, u\n" * s * "\nend"

    quote
        $(type)($output_syms, $input_syms, $s) do $(esc(out_var)), $(esc(u))
            $(body...)
            nothing
        end
    end
end

# for metadata check, validates both input and output symbols
function assert_initformula_compat(cf::ComponentModel, c::InitFormula)
    settable_symbols = Set(vcat(
        sym(cf),
        psym(cf),
        insym_flat(cf),
        outsym_flat(cf)
    ))
    input_symbols = settable_symbols ∪ obssym(cf)

    input_mismatch = setdiff(c.sym, input_symbols)
    if !isempty(input_mismatch)
        throw(ArgumentError("InitFormula uses input symbols not part of component model: $input_mismatch (observables not allowed)"))
    end

    missing_symbols = setdiff(c.outsym, settable_symbols)
    if !isempty(missing_symbols)
        throw(ArgumentError("InitFormula output symbols must be existing component symbols (not observables), but these are not: $missing_symbols"))
    end

    c
end
function assert_guessformula_compat(cf::ComponentModel, c::GuessFormula)
    readable_symbols = Set(vcat(
        sym(cf),
        psym(cf),
        insym_flat(cf),
        outsym_flat(cf),
    ))

    guessable_symbols = Set(vcat(
        sym(cf),
        psym(cf),
        insym_flat(cf),
        outsym_flat(cf)
    ))

    input_mismatch = setdiff(c.sym, readable_symbols)
    if !isempty(input_mismatch)
        throw(ArgumentError("GuessFormula uses input symbols not part of component model: $input_mismatch (observables not allowed)"))
    end

    output_mismatch = setdiff(c.outsym, guessable_symbols)
    if !isempty(output_mismatch)
        throw(ArgumentError("GuessFormula output symbols must be guessable (not observables): $output_mismatch"))
    end

    c
end

"""
    topological_sort_formulas(formulas::Vector{Init/GuessFormula}) -> Vector{Init/GuessFormula}

Sort formulas in topological order based on their dependencies. This ensures that
formulas are applied in the correct order, where each formula's input symbols are
available before it executes.

A formula B depends on formula A if any of B's input symbols are produced by A's output symbols.
The function returns the formulas reordered such that dependencies are satisfied.
"""
function topological_sort_formulas(formulas)
    n = length(formulas)
    n <= 1 && return copy(formulas)

    type = if first(formulas) isa InitFormula
        "InitFormula"
    elseif first(formulas) isa GuessFormula
        "GuessFormula"
    else
        throw(ArgumentError("Expected a vector of InitFormula or GuessFormula, got $(typeof(first(formulas)))"))
    end

    # Check for output symbol conflicts
    all_outputs = Symbol[]
    for formula in formulas
        append!(all_outputs, formula.outsym)
    end

    if !allunique(all_outputs)
        conflicts = [s for s in unique(all_outputs) if count(==(s), all_outputs) > 1]
        throw(ArgumentError("Multiple $(type)s set the same symbol(s): $conflicts"))
    end

    # Build dependency graph using Graphs.jl
    g = SimpleDiGraph(n)

    # Add edges: i → j means formula j depends on formula i
    # (formula j must run after formula i)
    for i in 1:n, j in 1:n
        if i != j
            # Check if formula j depends on formula i
            # i.e., j's input symbols intersect with i's output symbols
            if !isdisjoint(formulas[j].sym, formulas[i].outsym)
                add_edge!(g, i, j)
            end
        end
    end

    # Perform topological sort
    try
        sorted_indices = Graphs.topological_sort(g)
        if length(sorted_indices) != n
            # This shouldn't happen with a proper topological sort implementation,
            # but let's be safe
            throw(ArgumentError("Circular dependency detected in $(type)s"))
        end
        return formulas[sorted_indices]
    catch e
        if e isa ErrorException && occursin("graph contains at least one loop", string(e))
            throw(ArgumentError("Circular dependency detected in $(type)s"))
        else
            rethrow(e)
        end
    end
end

function apply_init_formulas!(defaults, formulas_unsorted; verbose=false)
    # Convert tuple to vector if necessary
    formulas_vec = formulas_unsorted isa Tuple ? collect(formulas_unsorted) : formulas_unsorted
    formulas = topological_sort_formulas(formulas_vec)

    for f in formulas
        out = SymbolicView(zeros(length(f.outsym)), f.outsym)
        # ensure all input symbols are in defaults
        invals = map(f.sym) do s
            if !haskey(defaults, s)
                throw(ArgumentError("InitFormula requires input symbol $s to be defined in defaults"))
            end
            defaults[s]
        end
        if any(v -> ismissing(v) || isnothing(v) || isnan(v), invals)
            throw(ArgumentError("InitFormula requires all input symbols to be initialized, but found NaN/missing/nothing in inputs: $(f.sym .=> invals)"))
        end
        in = SymbolicView(invals, f.sym)
        f(out, in)
        for s in f.outsym
            val = out[s]
            if verbose
                if haskey(defaults, s)
                    if defaults[s] ≈ val
                        println("InitFomula: keeping default for :$s at $(val)")
                    else
                        println("InitFomula: updating default for :$s from $(defaults[s]) to $(val)")
                    end
                else
                    println("InitFomula: setting default for :$s to $(val)")
                end
            end
            defaults[s] = val
        end
    end
    return defaults
end
function apply_guess_formulas!(guesses, defaults, formulas_unsorted; verbose=false)
    # Convert tuple to vector if necessary
    formulas_vec = formulas_unsorted isa Tuple ? collect(formulas_unsorted) : formulas_unsorted
    formulas = topological_sort_formulas(formulas_vec)

    for f in formulas
        out = SymbolicView(zeros(length(f.outsym)), f.outsym)
        # Layered lookup: defaults (fixed) take precedence over guesses
        invals = map(f.sym) do s
            if haskey(defaults, s)
                defaults[s]  # Use fixed default value if available
            elseif haskey(guesses, s)
                guesses[s]   # Otherwise use guess value
            else
                throw(ArgumentError("GuessFormula requires input symbol $s to be defined in defaults or guesses"))
            end
        end
        # Validate inputs are not NaN/missing/nothing
        if any(v -> ismissing(v) || isnothing(v) || isnan(v), invals)
            throw(ArgumentError("GuessFormula requires all input symbols to be initialized, but found NaN/missing/nothing in inputs: $(f.sym .=> invals)"))
        end
        in = SymbolicView(invals, f.sym)
        f(out, in)
        # Update guesses dictionary (NOT defaults!)
        for s in f.outsym
            val = out[s]
            if verbose
                if haskey(defaults, s)
                    # This symbol is fixed, so updating guess won't affect the solve
                    println("GuessFormula: symbol :$s has default $(defaults[s]), guess update to $(val) will have no effect")
                elseif haskey(guesses, s)
                    if guesses[s] ≈ val
                        println("GuessFormula: keeping guess for :$s at $(val)")
                    else
                        println("GuessFormula: updating guess for :$s from $(guesses[s]) to $(val)")
                    end
                else
                    println("GuessFormula: setting guess for :$s to $(val)")
                end
            end
            guesses[s] = val
        end
    end
    return guesses
end

_vcattable(t::Tuple) = collect(t)
_vcattable(x) = x # Formula types and AbstractVector

# Generic formula collection function (works for InitFormula, GuessFormula, etc.)
collect_formulas(::Nothing, ::Nothing) = nothing
collect_formulas(::Nothing, f) = vcat(_vcattable(f))
collect_formulas(f, ::Nothing) = vcat(_vcattable(f))
collect_formulas(a, b) = vcat(_vcattable(a), _vcattable(b))

merge_initconstraints(::Nothing, ::Nothing) = nothing
merge_initconstraints(::Nothing, c::InitConstraint) = c
merge_initconstraints(::Nothing, cs) = InitConstraint(cs...)

merge_initconstraints(c ::InitConstraint, ::Nothing) = c
merge_initconstraints(cA::InitConstraint, cB::InitConstraint) = InitConstraint(cA, cB)
merge_initconstraints(c ::InitConstraint, cs) = InitConstraint(c, cs...)

merge_initconstraints(cs, ::Nothing) = InitConstraint(cs...)
merge_initconstraints(cs, c::InitConstraint) = InitConstraint(cs..., c)
merge_initconstraints(csA, csB) = InitConstraint(csA..., csB...)
