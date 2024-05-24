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

function Base.show(io::IO, ::MIME"text/plain", c::EdgeFunction)
    type = match(r"^(.*?)\{", string(typeof(c)))[1]
    println(io, type, styled" :$(c.name) with $(_styled_coupling(coupling(c))) coupling")

    styling = Dict{Int, Symbol}()
    if coupling(c) == Fiducial()
        for i in 1:depth(c)
            styling[i] = :NetworkDynamics_fordst
        end
        for i in depth(c)+1:2*depth(c)
            styling[i] = :NetworkDynamics_forsrc
        end
    elseif coupling(c) == Directed()
        for i in 1:depth(c)
            styling[i] = :NetworkDynamics_fordst
        end
    else
        for i in 1:depth(c)
            styling[i] = :NetworkDynamics_fordstsrc
        end
    end

    print_states_params(io, c, styling)
end
_styled_coupling(::Fiducial) = styled"{NetworkDynamics_fordst:Fidu}{NetworkDynamics_forsrc:cial}"
_styled_coupling(::Directed) = styled"{NetworkDynamics_fordst:Directed}"
_styled_coupling(::AntiSymmetric) = styled"{NetworkDynamics_fordstsrc:AntiSymmetric}"
_styled_coupling(::Symmetric) = styled"{NetworkDynamics_fordstsrc:Symmetric}"

function Base.show(io::IO, ::MIME"text/plain", c::VertexFunction)
    type = match(r"^(.*?)\{", string(typeof(c)))[1]
    println(io, type, styled" :$(c.name)")

    styling = Dict{Int, Symbol}()
    for i in 1:depth(c)
        styling[i] = :NetworkDynamics_forlayer
    end

    print_states_params(io, c, styling)
end

function print_states_params(io, c::ComponentFunction, styling)
    info = Base.AnnotatedString{String}[]
    num, word = maybe_plural(dim(c), "state")
    push!(info, styled"$num $word:\t$(stylesymbolarray(c.sym, c.def, styling))")

    if hasproperty(c, :mass_matrix) && c.mass_matrix != LinearAlgebra.I
        info[end] *= "\nmass m:\t$(c.mass_matrix)"
    end

    num, word = maybe_plural(pdim(c), "param")
    pdim(c) > 0 && push!(info, styled"$num $word:\t$(stylesymbolarray(c.psym, c.pdef))")

    print_treelike(io, info)
end

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
            ret = ret*styled"{NetworkDynamics_defaultval:=$default}"
        end
        if i < length(syms)
            ret = ret*", "
        end
    end
    ret = ret*"]"
end

function print_treelike(io, vec, prefix=" ", infix=" ")
    for s in @views vec[begin:end-1]
        sublines = split(s,"\n")
        println(io, prefix, "├─", infix, sublines[begin])
        for sub in sublines[begin+1:end]
            println(io, prefix, "│ ", infix, sub)
        end
    end
    sublines = split(vec[end],"\n")
    println(io, prefix, "└─", infix, sublines[begin])
    for sub in sublines[begin+1:end]
        println(io, prefix, "  ", infix, sub)
    end
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
