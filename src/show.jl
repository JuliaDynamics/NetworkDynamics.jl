# AbstractTrees.TreeCharSet("├", "└","┌" "│", "─", "⋮", " ⇒ ")

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
    end
end

Base.show(io::IO, s::NaiveAggregator) = print(io, "NaiveAggregator($(repr(s.f)))")
Base.show(io::IO, s::KAAggregator) = print(io, "KAAggregator($(repr(s.f)))")
Base.show(io::IO, s::SequentialAggregator) = print(io, "SequentialAggregator($(repr(s.f)))")
Base.show(io::IO, s::PolyesterAggregator) = print(io, "PolyesterAggregator($(repr(s.f)))")

function Base.show(io::IO, ::MIME"text/plain", @nospecialize(c::ComponentFunction))
    type = match(r"^(.*?)\{", string(typeof(c)))[1]
    print(io, type, styled" {NetworkDynamics_name::$(c.name)}")
    print(io, styled" {NetworkDynamics_fftype:$(fftype(c))}")
    if has_graphelement(c)
        ge = get_graphelement(c)
        if c isa VertexFunction
            print(io, " @ Vertex $ge")
        else
            print(io, " @ Edge $(ge.src)=>$(ge.dst)")
        end
    end

    styling = Dict{Int,Symbol}()

    print_states_params(io, c, styling)
end

function print_states_params(io, @nospecialize(c::ComponentFunction), styling)
    info = AnnotatedString{String}[]

    if hasinsym(c)
        push!(info, _inout_string(c, insym, "input"))
    end

    num, word = maybe_plural(dim(c), "state")
    push!(info, styled"$num &$word: &&$(stylesymbolarray(c.sym, def(c), guess(c), styling))")

    if hasproperty(c, :mass_matrix) && c.mass_matrix != LinearAlgebra.I
        if LinearAlgebra.isdiag(c.mass_matrix) && !(c.mass_matrix isa UniformScaling)
            info[end] *= "\n&with diagonal mass matrix $(LinearAlgebra.diag(c.mass_matrix))"
        else
            info[end] *= "\n&with mass matrix $(c.mass_matrix)"
        end
    end

    push!(info, _inout_string(c, outsym, "output"))

    num, word = maybe_plural(pdim(c), "param")
    pdim(c) > 0 && push!(info, styled"$num &$word: &&$(stylesymbolarray(c.psym, pdef(c), pguess(c)))")

    print_treelike(io, align_strings(info))
end
function _inout_string(@nospecialize(c::VertexFunction), f, name)
    sym = f(c)
    num, word = maybe_plural(length(sym), name)
    defs = get_defaults_or_inits(c, sym)
    guesses = get_guesses(c, sym)
    styled"$num &$word: &&$(stylesymbolarray(sym, defs, guesses))"
end
function _inout_string(@nospecialize(c::EdgeFunction), f, name)
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

function stylesymbolarray(syms, defaults, guesses, symstyles=Dict{Int,Symbol}())
    @assert length(syms) == length(defaults)
    ret = "["
    for (i, sym, default, guess) in zip(1:length(syms), syms, defaults, guesses)
        style = get(symstyles, i, :default)
        ret = ret * styled"{$style:$(string(sym))}"
        if !isnothing(default)
            _str = str_significant(default; sigdigits=2)
            ret = ret * styled"{NetworkDynamics_defaultval:=$(_str)}"
        elseif !isnothing(guess)
            _str = str_significant(guess; sigdigits=2)
            ret = ret * styled"{NetworkDynamics_guessval:≈$(_str)}"
        end
        if i < length(syms)
            ret = ret * ", "
        end
    end
    ret = ret * "]"
end

function Base.show(io::IO, idx::VIndex)
    print(io, "VIndex(", repr(idx.compidx), ", ", repr(idx.subidx), ")")
end
function Base.show(io::IO, idx::EIndex)
    print(io, "EIndex(", repr(idx.compidx), ", ", repr(idx.subidx), ")")
end
function Base.show(io::IO, idx::VPIndex)
    print(io, "VPIndex(", repr(idx.compidx), ", ", repr(idx.subidx), ")")
end
function Base.show(io::IO, idx::EPIndex)
    print(io, "EPIndex(", repr(idx.compidx), ", ", repr(idx.subidx), ")")
end

function Base.show(io::IO, mime::MIME"text/plain", @nospecialize(s::NWState); dim=nothing)
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
            str = read(seekstart(buf), AnnotatedString)
            if !isnothing(dim) && dim(sym)
                str = styled"{NetworkDynamics_inactive:$str}"
            end
            str
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

function Base.show(io::IO, mime::MIME"text/plain", @nospecialize(p::NWParameter); dim=nothing)
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
            str = read(seekstart(buf), AnnotatedString)
            if !isnothing(dim) && dim(sym)
                str = styled"{NetworkDynamics_inactive:$str}"
            end
            str
        end
        print_treelike(io, align_strings(strvec); prefix="  ", rowmax)
    end
end

function Base.show(io::IO, mime::MIME"text/plain", p::Union{VProxy,EProxy})
    name = _proxyname(p)
    print(io, styled"{bright_blue:$name} for ")
    show(io, mime, p.s; dim = _dimcondition(p))
end
_proxyname(::VProxy) = "Vertex indexer"
_proxyname(::EProxy) = "Edge indexer"
_dimcondition(::VProxy) = sym -> sym isa EIndex || sym isa EPIndex
_dimcondition(::EProxy) = sym -> sym isa VIndex || sym isa VPIndex

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
    if !isempty(newlines)
        for i in newlines
            alig[i] = alig[i] * "\n" * alig[i+1]
        end
        deleteat!(alig, newlines .+1 )
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

function str_significant(x; sigdigits)
    (x == 0) && (return "0")
    x = round(x; sigdigits)
    n = length(@sprintf("%d", abs(x)))              # length of the integer part
    if (x ≤ -1 || x ≥ 1)
        decimals = max(sigdigits - n, 0)               # 'sig - n' decimals needed
    else
        Nzeros = ceil(Int, -log10(abs(x))) - 1      # No. zeros after decimal point before first number
        decimals = sigdigits + Nzeros
    end
    return @sprintf("%.*f", decimals, x)
end
