"""
    struct InitConstraint{F}
    InitConstraint(f, sym, dim)

A representation of an additionalconstraint that is applied during the initialization phase of a component.
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
