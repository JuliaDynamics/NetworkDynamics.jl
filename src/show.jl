# AbstractTrees.TreeCharSet("├", "└", "│", "─", "⋮", " ⇒ ")

function Base.show(io::IO, ::MIME"text/plain", nw::Network)
    println(io, "Dynamic network with:")
    num, word = maybe_plural(length(nw.vertexbatches), "type")
    println(io, " ├─ $(nv(nw.im.g)) vertices ($num unique $word)")
    num, word = maybe_plural(length(nw.layer.edgebatches), "type")
    println(io, " └─ $(ne(nw.im.g)) edges ($num unique $word)")
    print(io, "Aggregation using ")
    print(io, nw.layer.aggregator)
    println(io, " along")
    num, word = maybe_plural(nw.layer.edepth, "dimension")
    println(io, " ├─ $(num) edge $(word)")
    num, word = maybe_plural(nw.layer.vdepth, "dimension")
    println(io, " └─ $(num) vertex $(word)")
end

Base.show(io::IO, s::NNlibScatter) = print(io, "NNlibScatter($(repr(s.f)))")
Base.show(io::IO, s::NaiveAggregator) = print(io, "NaiveAggregator($(repr(s.f)))")
Base.show(io::IO, s::KAAggregator) = print(io, "KAAggregator($(repr(s.f)))")
Base.show(io::IO, s::SequentialAggregator) = print(io, "SequentialAggregator($(repr(s.f)))")
Base.show(io::IO, s::PolyesterAggregator) = print(io, "PolyesterAggregator($(repr(s.f)))")

function Base.show(io::IO, ::MIME"text/plain", c::StaticEdge)
    println(io, "StaticEdge :$(c.name) with $(coupling(c)) coupling")
    info = String[]
    push!(info, "$(dim(c)) states:\t$(c.sym)")

    pdim(c) > 0 && push!(info, "$(pdim(c)) params:\t$(c.psym)")

    print_treelike(io, info)
end

StyledStrings.addface!(:defaultval => StyledStrings.Face(color = :blue, inherit = :emphasis))

function stylesymbolarray(syms, defaults, symstyles=Dict{Int,Symbol}())
    @assert length(syms) == length(defaults)
    ret = "["
    for (i, sym, default) in zip(1:length(syms), syms,defaults)
        style = get(symstyles, i, nothing)
        if isnothing(style)
            ret = ret*string(sym)
        else
            ret = ret*styled"{$style:$(string(sym))}"
        end
        if !isnothing(default)
            ret = ret*styled"{defaultval:=$default}"
        end
        ret = ret*", "
    end
    ret = ret[1:end-2]*"]"
end

function print_treelike(io, vec, prefix=" ", infix=" ")
    for s in @views vec[begin:end-1]
        println(io, prefix, "├─", infix, s)
    end
    println(io, prefix, "└─", infix, vec[end])
end

function maybe_plural(num, word, substitution=s"\1s")
    if num > 1
        word = replace(word, r"^(.*)$"=>substitution)
    end
    num, word
end

# function Base.show(io::IO, ::MIME"text/plain", m::MyType)
#     print(io, "Examplary instance of MyType\n", m.x, " ± ", m.y)
# end

# Base.show(io::IO, m::MyType) = print(io, m.x, '(', m.y, ')')
