# AbstractTrees.TreeCharSet("├", "└","┌" "│", "─", "⋮", " ⇒ ")

# overload to reduce excessive printing
function Base.show(io::IO, @nospecialize(nw::Network))
    print(io, "Network($(nv(nw.im.g)) vertices, $(ne(nw.im.g)) edges)")
end

function Base.show(io::IO, ::MIME"text/plain", @nospecialize(nw::Network))
    compact = get(io, :compact, false)::Bool
    if compact
        print(io, "Network ($(nv(nw.im.g)) vertices, $(ne(nw.im.g)) edges)")
    else
        print(io, "Network")
        num, word = maybe_plural(dim(nw), "state")
        print(io, " with $num $word")
        num, word = maybe_plural(pdim(nw), "parameter")
        print(io, " and $num $word")
        num, word = maybe_plural(length(nw.vertexbatches), "type")
        print(io, "\n ├─ $(nv(nw.im.g)) vertices ($num unique $word)")
        num, word = maybe_plural(length(nw.layer.edgebatches), "type")
        print(io, "\n └─ $(ne(nw.im.g)) edges ($num unique $word)")
        print(io, "\nEdge-Aggregation using ")
        print(io, nw.layer.aggregator)
        # print(io, "\n ├─ vertex output dimension: $(nw.im.vdepth)")
        # print(io, "\n └─   edge output dimension: $(nw.im.edepth)")

        if !isnothing(nw.jac_prototype)
            jac = nw.jac_prototype
            rate = SparseArrays.nnz(jac) / dim(nw)^2
            print(io, "\nJacobian prototype defined: ")
            printstyled(io, (round(rate*100, digits=2)), "% sparsity", color=:blue)
        end

        Ncb = length(wrap_component_callbacks(nw))
        if Ncb ≥ 1
            nvert = mapreduce(has_callback, +, nw.im.vertexm)
            nedge = mapreduce(has_callback, +, nw.im.edgem)
            _, setword = maybe_plural(Ncb, "callback set")
            _, vertword = maybe_plural(nvert, "vertex", "vertices")
            _, edgeword = maybe_plural(nedge, "edge")
            print(io, "\n$Ncb $setword across $nvert $vertword and $nedge $edgeword")
        end
    end
end

Base.show(io::IO, s::NaiveAggregator) = print(io, "NaiveAggregator($(repr(s.f)))")
Base.show(io::IO, s::KAAggregator) = print(io, "KAAggregator($(repr(s.f)))")
Base.show(io::IO, s::SequentialAggregator) = print(io, "SequentialAggregator($(repr(s.f)))")
Base.show(io::IO, s::PolyesterAggregator) = print(io, "PolyesterAggregator($(repr(s.f)))")

function Base.show(io::IO, mime::MIME"text/plain", c::ComponentModel)
    type = string(typeof(c))
    print(io, type, styled" {NetworkDynamics_name::$(c.name)}")
    print(io, styled" {NetworkDynamics_fftype:$(fftype(c))}")
    if has_graphelement(c)
        ge = get_graphelement(c)
        if c isa VertexModel
            print(io, " @ Vertex $ge")
        else
            print(io, " @ Edge $(ge.src)=>$(ge.dst)")
        end
    end

    styling = Dict{Int,Symbol}()

    print_states_params(io, c, styling)

    # special known metadata of PowerDynamics
    if has_metadata(c, :pfmodel)
        pfm = get_metadata(c, :pfmodel)
        printstyled(io, "\nPowerflow model ", color=:blue, bold=true)
        print(io, ":", pfm.name)
        pdim(pfm) > 0 && print(io, styled" with $(stylesymbolarray(pfm.psym, _pdef(pfm), _pguess(pfm)))")
    end
end
# enable compact printing in lists of component models
function Base.show(io::IO, c::ComponentModel)
    type = string(typeof(c))
    print(io, type, styled" {NetworkDynamics_name::$(c.name)}")
    if has_graphelement(c)
        ge = get_graphelement(c)
        if c isa VertexModel
            print(io, " @ Vertex $ge")
        else
            print(io, " @ Edge $(ge.src)=>$(ge.dst)")
        end
    end
end

function print_states_params(io, @nospecialize(c::ComponentModel), styling)
    info = AnnotatedString{String}[]

    if hasinsym(c)
        push!(info, _inout_string(c, insym, "input"))
    end

    num, word = maybe_plural(dim(c), "state")
    push!(info, styled"$num &$word: &&$(stylesymbolarray(c.sym, _def(c), _guess(c), styling))")

    if hasproperty(c, :mass_matrix) && c.mass_matrix != LinearAlgebra.I
        if LinearAlgebra.isdiag(c.mass_matrix) && !(c.mass_matrix isa UniformScaling)
            info[end] *= "\n&with diagonal mass matrix $(LinearAlgebra.diag(c.mass_matrix))"
        else
            info[end] *= "\n&with mass matrix $(c.mass_matrix)"
        end
    end

    push!(info, _inout_string(c, outsym, "output"))

    num, word = maybe_plural(pdim(c), "param")
    pdim(c) > 0 && push!(info, styled"$num &$word: &&$(stylesymbolarray(c.psym, _pdef(c), _pguess(c)))")

    if has_external_input(c)
        num = extdim(c)
        arr = match(r"(\[.*\])", repr(extin(c)))[1]
        push!(info, styled"$num &ext in: &&$arr")
    end

    if has_initformula(c)
        formulas = get_initformulas(c)
        total_eqs = sum(length(formula.outsym) for formula in formulas)
        num, word = maybe_plural(total_eqs, "eq.", "eqs.")
        all_outsyms = reduce(vcat, [f.outsym for f in formulas])
        formula_word = length(formulas) == 1 ? "formula" : "formulas"
        str = "$num &add. init $word from $(length(formulas)) $formula_word setting $(all_outsyms)"
        push!(info, str)
    end

    if has_initconstraint(c)
        constraints = get_initconstraints(c)
        total_eqs = sum(constraint.dim for constraint in constraints)
        num, word = maybe_plural(total_eqs, "eq.", "eqs.")
        all_syms = unique(reduce(vcat, [constraint.sym for constraint in constraints]))
        constraint_word = length(constraints) == 1 ? "constraint" : "constraints"
        str = "$num &add. init $word from $(length(constraints)) $constraint_word for $(all_syms)"
        push!(info, str)
    end

    # PowerDynamics metadata display
    if has_metadata(c, :pfinitformula)
        formulas = get_metadata(c, :pfinitformula)
        formulas_tuple = formulas isa Tuple ? formulas : (formulas,)
        total_eqs = sum(length(formula.outsym) for formula in formulas_tuple)
        num, word = maybe_plural(total_eqs, "pf init eq.", "pf init eqs.")
        formula_word = length(formulas_tuple) == 1 ? "formula" : "formulas"
        all_outsyms = unique(reduce(vcat, [formula.outsym for formula in formulas_tuple]))
        all_pfsyms = unique(reduce(vcat, [formula.pfsym for formula in formulas_tuple]))
        str = "$num &add. $word from $(length(formulas_tuple)) $formula_word setting $(all_outsyms) using @pf($(all_pfsyms))"
        push!(info, str)
    end

    if has_metadata(c, :pfinitconstraint)
        constraints = get_metadata(c, :pfinitconstraint)
        constraints_tuple = constraints isa Tuple ? constraints : (constraints,)
        total_eqs = sum(constraint.dim for constraint in constraints_tuple)
        num, word = maybe_plural(total_eqs, "add. pf init eq.", "add. pf init eqs.")
        constraint_word = length(constraints_tuple) == 1 ? "constraint" : "constraints"
        all_syms = unique(reduce(vcat, [constraint.sym for constraint in constraints_tuple]))
        all_pfsyms = unique(reduce(vcat, [constraint.pfsym for constraint in constraints_tuple]))
        str = "$num &add. $word from $(length(constraints_tuple)) $constraint_word for $(all_syms) using @pf($(all_pfsyms))"
        push!(info, str)
    end

    if has_callback(c)
        cbs = get_callbacks(c)
        num, word = maybe_plural(length(cbs), "callback")
        str = ""
        str *= styled"$num &$word: &&$(shortrepr(first(cbs)))"
        for cb in Base.tail(cbs)
            str *= styled"\n&&&$(shortrepr(cb))"
        end
        push!(info, str)
    end

    print_treelike(io, align_strings(info))
end
function _inout_string(c::VertexModel, f, name)
    sym = f(c)
    num, word = maybe_plural(length(sym), name)
    defs = get_defaults_or_inits(c, sym)
    guesses = get_guesses(c, sym)
    styled"$num &$word: &&$(stylesymbolarray(sym, defs, guesses))"
end
function _inout_string(c::EdgeModel, f, name)
    sym = f(c)
    word = name*"s"
    srcnum = length(sym.src)
    dstnum = length(sym.dst)
    srcdefs = get_defaults_or_inits(c, sym.src)
    dstdefs = get_defaults_or_inits(c, sym.dst)
    srcguesses = get_guesses(c, sym.src)
    dstguesses = get_guesses(c, sym.dst)
    styled"$srcnum/$dstnum &$word: &&src=&$(stylesymbolarray(sym.src, srcdefs, srcguesses)) \
            dst=$(stylesymbolarray(sym.dst, dstdefs, dstguesses))"
end
function _def(c::ComponentModel)::Vector{Union{Nothing,Float64}}
    map(c.sym) do s
        has_default_or_init(c, s) ? get_default_or_init(c, s) : nothing
    end
end
function _guess(c::ComponentModel)::Vector{Union{Nothing,Float64}}
    map(c.sym) do s
        has_guess(c, s) ? get_guess(c, s) : nothing
    end
end
function _pdef(c::ComponentModel)::Vector{Union{Nothing,Float64}}
    map(c.psym) do s
        has_default_or_init(c, s) ? get_default_or_init(c, s) : nothing
    end
end
function _pguess(c::ComponentModel)::Vector{Union{Nothing,Float64}}
    map(c.psym) do s
        has_guess(c, s) ? get_guess(c, s) : nothing
    end
end

function stylesymbolarray(syms, defaults, guesses, symstyles=Dict{Int,Symbol}())
    @assert length(syms) == length(defaults)
    ret = "["
    for (i, sym, default, guess) in zip(1:length(syms), syms, defaults, guesses)
        style = get(symstyles, i, :default)
        ret = ret * styled"{$style:$(string(sym))}"
        if !isnothing(default)
            _str = str_significant(default; sigdigits=5)
            ret = ret * styled"{NetworkDynamics_defaultval:=$(_str)}"
        elseif !isnothing(guess)
            _str = str_significant(guess; sigdigits=5)
            ret = ret * styled"{NetworkDynamics_guessval:≈$(_str)}"
        end
        if i < length(syms)
            ret = ret * ", "
        end
    end
    ret = ret * "]"
end

function Base.show(io::IO, idx::VIndex)
    print(io, "VIndex(", repr_colon(idx.compidx), ", ", repr_colon(idx.subidx), ")")
end
function Base.show(io::IO, idx::EIndex)
    print(io, "EIndex(", repr_colon(idx.compidx), ", ", repr_colon(idx.subidx), ")")
end
function Base.show(io::IO, idx::VIndex{<:Any,Nothing})
    print(io, "VIndex(", repr_colon(idx.compidx), ")")
end
function Base.show(io::IO, idx::EIndex{<:Any,Nothing})
    print(io, "EIndex(", repr_colon(idx.compidx), ")")
end
Base.show(io::IO, idx::ParamIdx) = print(io, "ParamIdx(", repr_colon(idx.idx), ")")
Base.show(io::IO, idx::StateIdx) = print(io, "StateIdx(", repr_colon(idx.idx), ")")
repr_colon(::Colon) = ":"
repr_colon(x) = repr(x)

function Base.show(io::IO, mime::MIME"text/plain", @nospecialize(s::NWState))
    ioc = IOContext(io, :compact => true)
    if get(io, :limit, true)::Bool
        dsize = get(io, :displaysize, displaysize(io))::Tuple{Int,Int}
        rowmax = dsize[1] - 8
    else
        rowmax = typemax(Int)
    end

    print(io, "NWState{$(typeof(uflat(s)))} of ")
    show(ioc, mime, s.nw)

    if isempty(s.uflat) > 0
       print(io, "\n u = ", s.uflat)
    else
        strvec = map(SII.variable_symbols(s.nw), eachindex(s.uflat)) do sym, i
            buf = AnnotatedIOBuffer()
            print(buf, "&")
            show(buf, mime, sym)
            print(buf, " &&=> ")
            if isassigned(s.uflat, i)
                show(buf, mime, s.uflat[i])
            else
                print(buf, "#undef")
            end
            read(seekstart(buf), AnnotatedString)
        end
        print_treelike(io, align_strings(strvec); prefix="  ", rowmax)
    end

    buf = IOContext(AnnotatedIOBuffer(), :compact=>true, :limit=>true)
    print(buf, "\n p = ")
    show(buf, mime, s.p)
    print(buf, "\n t = ")
    show(buf, mime, s.t)
    str = read(seekstart(buf.io), AnnotatedString)
    if !isnothing(dim)
        print(io, styled"{NetworkDynamics_inactive:$str}")
    else
        print(io, str)
    end
end

function Base.show(io::IO, mime::MIME"text/plain", @nospecialize(p::NWParameter))
    compact = get(io, :compact, false)::Bool
    if get(io, :limit, true)::Bool
        dsize = get(io, :displaysize, displaysize(io))::Tuple{Int,Int}
        rowmax = dsize[1] - 5
    else
        rowmax = typemax(Int)
    end

    if compact
        print(io, "NWParameter(")
        show(io, p.pflat)
        print(io, ")")
    else
        ioc = IOContext(io, :compact => true)
        print(io, "Parameter{$(typeof(pflat(p)))} of ")
        show(ioc, mime, p.nw)

        strvec = map(SII.parameter_symbols(p.nw), eachindex(p.pflat)) do sym, i
            buf = AnnotatedIOBuffer()
            print(buf, "&")
            show(buf, mime, sym)
            print(buf, " &&=> ")
            if isassigned(p.pflat, i)
                show(buf, mime, p.pflat[i])
            else
                print(buf, "#undef")
            end
            read(seekstart(buf), AnnotatedString)
        end
        print_treelike(io, align_strings(strvec); prefix="  ", rowmax)
    end
end

function Base.show(io::IO, mime::MIME"text/plain", fp::FilteringProxy)
    printstyled(io, "FilteringProxy"; bold=true)
    print(io, " for ")
    if fp.data isa NWState
        print(io, "NWState()")
    elseif fp.data isa NWParameter
        print(io, "NWParameter()")
    else
        print(io, repr(typeof(fp.data)))
    end

    print(io, "\n  Component filter: ")
    if isnothing(fp.compfilter)
        printstyled(io, "none"; color=:light_black)
        printstyled(io, " <- filter further by obj.v, obj.e, obj[VIndex(..)], ..."; color=:light_black)
    else
        show(io, mime, fp.compfilter)
        if fp.compfilter isa Union{AllVertices, AllEdges}
            printstyled(io, " <- filter further by obj[1], obj[\"name\"], ..."; color=:light_black)
        end
    end
    print(io, "\n  State filter:     ")
    if isnothing(fp.varfilter)
        printstyled(io, "none"; color=:light_black)
        if fp.compfilter isa SymbolicIndex
            printstyled(io, " <- filter states by obj[\"δ\"], obj[:x], ..."; color=:light_black)
        end
    else
        show(io, mime, fp.varfilter)
    end
    print(io, "\n  Types:")

    # for c in [:blue, :cyan, :green, :hidden, :light_black, :light_blue, :light_cyan, :light_green, :light_magenta, :light_red, :light_white, :light_yellow, :magenta, :nothing, :red, :reverse, :underline, :white, :yellow]
    #    printstyled("\n",c, color=c)
    # end
    color_map = Dict(
        :states => :light_green,
        :parameters => :light_yellow,
        :inputs => :light_magenta,
        :outputs => :light_blue,
        :observables => :light_cyan,
    )

    for s in [:states, :parameters, :inputs, :outputs, :observables]
        if getfield(fp, s)
            printstyled(io, "  $s ✓", bold=true, color=color_map[s])
        else
            printstyled(io, "  $s ✗", color=:light_black)
        end
    end

    indices, types = generate_indices(fp; return_types=true)
    _val_to_str(x) = isnan(x) ? "NaN (undefined?)" : repr(x)
    value_str = try
        # try to resolve all at once (faster)
        values = fp.data[indices]
        map(_val_to_str, indices)
    catch
        # alternatively, go through them one by one
        map(indices) do _idx
            try
                _val_to_str(fp.data[_idx])
            catch
                "not accessible"
            end
        end
    end
    nw = extract_nw(fp)
    compnames = [string(getcomp(nw,idx).name) for idx in indices]

    if !isempty(indices)
        printstyled(io, "\nMatching Indices:", bold=true)
        # align all the printouts
        strvec = Ref("&") .* repr.(indices) .* Ref(" &&=> ") .* value_str .* " &:" .* compnames
        aligned_strvec =  align_strings(strvec)
        for i in eachindex(indices)
            isfirst = i == 1 ||
                      idxtype(indices[i]) != idxtype(indices[i-1]) ||
                      indices[i].compidx != indices[i-1].compidx
            islast  = i == length(indices) ||
                      idxtype(indices[i]) != idxtype(indices[i+1]) ||
                      indices[i].compidx != indices[i+1].compidx
            if isfirst && islast
                # ·  ∙  •  ●  ◦  Ø  ø
                # □  ■  ▫  ▪  ◆  ◊
                print("\n  ● ")
            elseif isfirst
                print("\n  ╭ ")
            elseif islast
                print("\n  ╰ ")
            else
                print("\n  │ ")
            end
            regex = r"(\s*[VE]Index\(.*?,\s*)(.*)(\)\s*=>.*)(:.*)\s*$"

            # t1 ="  VIndex(1, :Pm) => 0.9645731407002397 (kuramoto_second)"
            # match(regex, t1)[1] == "  VIndex(1, "
            # match(regex, t1)[2] == ":Pm"
            # t1 == match(regex, t1)[1]*match(regex, t1)[2]*match(regex, t1)[3]*match(regex, t1)[4]
            # t2 =" VIndex(4, :δ)  => 1.2211532806955903 (kuramoto_second)"
            # match(regex, t2)[1] == " VIndex(4, "
            # match(regex, t2)[2] == ":δ"
            # t2 == match(regex, t2)[1]*match(regex, t2)[2]*match(regex, t2)[3]*match(regex, t2)[4]
            m = match(regex, aligned_strvec[i])
            if !isnothing(m)
                print(m[1])
                printstyled(m[2], color=color_map[types[i]])
                print(m[3])
                isfirst && print(m[4])
            else
                print(aligned_strvec[i]) # fallback if regex fails
            end
        end
    else
        printstyled(io, "\nNo indices matching filter!", bold=true)
    end
end
Base.show(io::IO, ::MIME"text/plain", ::AllVertices) = print(io, "AllVertices()")
Base.show(io::IO, ::MIME"text/plain", ::AllEdges) = print(io, "AllEdges()")

function print_treelike(io, vec; prefix=" ", infix=" ", rowmax=typemax(Int))
    Base.require_one_based_indexing(vec)
    if length(vec) > rowmax
        upper = ceil(Int, rowmax / 2)
        lower = floor(Int, rowmax / 2) - 1
        for s in @views vec[1:upper]
            _print_tree_element(io, "├─", "| ", s, prefix, infix)
        end
        print(io, "\n", prefix, "⋮")
        for s in @views vec[end-lower:end-1]
            _print_tree_element(io, "├─", "| ", s, prefix, infix)
        end
    else
        for s in @views vec[begin:end-1]
            _print_tree_element(io, "├─", "| ", s, prefix, infix)
        end
    end
    _print_tree_element(io, "└─", "  ", vec[end], prefix, infix)
end

function _print_tree_element(io, treesym, sublinesym, s, prefix, infix)
    sublines = split(s, "\n")
    print(io, "\n", prefix, treesym, infix, sublines[begin])
    for sub in sublines[begin+1:end]
        print(io, "\n", prefix, sublinesym, infix, sub)
    end
end

function align_strings(_vecofstr::AbstractVector{<:AbstractString}; padding=:alternating)
    # FIXME: workaround for github.com/JuliaLang/StyledStrings.jl/issues/64
    @static if VERSION > v"1.10"
        vecofstr = _vecofstr
    else
        vecofstr = String.(_vecofstr)
    end
    splitted = Vector{AbstractString}[]
    sizehint!(splitted, length(vecofstr))
    i = 1
    newlines = Int[]
    for str in vecofstr
        for (ln, linesplit) in enumerate(eachsplit(str, '\n'))
            push!(splitted, split(linesplit, '&'))
            if ln > 1
                push!(newlines, i-1)
            end
            i += 1
        end
    end
    alig = align_strings(splitted; padding)

    @assert issorted(newlines)
    while !isempty(newlines)
        i = pop!(newlines)
        alig[i] = alig[i] * "\n" * alig[i+1]
        deleteat!(alig, i+1)
    end
    alig
end
function align_strings(vecofvec::AbstractVector{<:AbstractVector}; padding=:alternating)
    depth = maximum(length.(vecofvec))
    maxlength = zeros(Int,depth)

    for i in eachindex(vecofvec)
        for j in 1:depth
            j ≥ length(vecofvec[i]) && continue
            maxlength[j] = max(maxlength[j], length(vecofvec[i][j]))
        end
    end
    map(vecofvec) do strvec
        mapreduce(*,zip(1:depth, strvec,maxlength)) do (d, str, l)
            if padding == :alternating
                pad = isodd(d) ? lpad : rpad
            elseif padding == :left
                pad = lpad
            else
                pad = rpad
            end
            pad(str, l)
        end
    end
end

function maybe_plural(num, word, substitution=s"\1s")
    if num > 1 || num == 0
        word = replace(word, r"^(.*)$" => substitution)
    end
    num, word
end

function str_significant(x; sigdigits, phantom_minus=false)
    # isinteger(x) && return string(Int(x))
    formatted = @sprintf("%.*g", sigdigits, x)
    if occursin(r"e\+0*", formatted)
        formatted = replace(formatted, r"e\+0*" => "e")
    elseif occursin(r"e\-0+", formatted)
        formatted = replace(formatted, r"e\-0+" => "e-")
    end
    if phantom_minus && x >= 0
        formatted = " " * formatted
    end
    formatted
end

function Base.show(io::IO, ::MIME"text/plain", @nospecialize(cb::ComponentCallback))
    basename = readuntil(IOBuffer(repr(cb)), '{')
    print(io, basename)
    print(io, "(")
    print(io, shortrepr(cb))
    if cb isa VectorContinuousComponentCallback
        print(io, ", len=", cb.len)
    end
    print(io, ")")
end

function shortrepr(cb::ComponentCallback)
    io = IOBuffer()
    _print_condsyms(io, cb)
    print(io, " affecting ")
    _print_affectsyms(io, cb)
    String(take!(io))
end
function shortrepr(cb::PresetTimeComponentCallback)
    io = IOBuffer()
    _print_affectsyms(io, cb)
    print(io, " affected at t=")
    print(io, cb.ts)
    String(take!(io))
end

function _print_condsyms(io, @nospecialize(cb::ComponentCallback))
    print(io, "(")
    _print_syms(io, cb.condition.sym, true)
    _print_syms(io, cb.condition.psym, isempty(cb.condition.sym))
    print(io, ")")
end
function _print_affectsyms(io, @nospecialize(cb::ComponentCallback))
    print(io, "(")
    _print_syms(io, cb.affect.sym, true)
    _print_syms(io, cb.affect.psym, isempty(cb.affect.sym))
    print(io, ")")
end
function _print_syms(io, syms::Tuple, isfirst)
    for s in syms
        if !isfirst
            print(io, ", ")
        end
        print(io, ':', s)
        isfirst = false
    end
end
