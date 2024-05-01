using Graphs.SimpleGraphs

function isvalid(g::SimpleGraph, colors)
    degrees = degree(g)
    usedc_nodes = [Vector{Int}(undef, d) for d in degrees]
    empty!.(usedc_nodes)

    for (c, e) in zip(colors, edges(g))
        usedc_src = usedc_nodes[e.src]
        usedc_dst = usedc_nodes[e.dst]
        if c ∈ usedc_src || c ∈ usedc_dst
            return false
        else
            push!(usedc_src, c)
            push!(usedc_dst, c)
        end
    end
    return true
end

function color_edges_greedy(g::SimpleGraph)
    degrees = degree(g)
    colors = zeros(Int, ne(g))
    usedc_nodes = [Vector{Int}(undef, d) for d in degrees]
    empty!.(usedc_nodes)

    maxd = maximum(degrees)
    cstat = zeros(Int, maxd)
    corder = collect(1:maxd+1)

    for (i, e) in enumerate(edges(g))
        usedc_src = usedc_nodes[e.src]
        usedc_dst = usedc_nodes[e.dst]

        # somewhat arbitrary: for big exclude lists it is beneficial to sort them first
        sortfirst = length(usedc_src) + length(usedc_dst) > 1000
        c = pickfrom_without(corder, usedc_src, usedc_dst; sortfirst)

        colors[i] = c
        push!(usedc_src, c)
        push!(usedc_dst, c)

        if c > length(cstat)
            @assert c == length(cstat) + 1
            resize!(cstat, c)
            cstat[c] = 0
            resize!(corder, c + 1)
            corder[c+1] = c + 1
        end
        cstat[c] += 1
        sortperm!(view(corder, 1:length(cstat)), cstat)
    end
    @assert isvalid(g, colors) "Could not find valid edge coloring 😱"
    return colors
end

@inline function pickfrom_without(X, exA, exB; sortfirst=false)
    if sortfirst
        ex = sort!(vcat(exA, exB))
        for x in X
            insorted(x, ex) || return x
        end
    else
        for x in X
            if x ∉ exA && x ∉ exB
                return x
            end
        end
    end
    error("did not find any element")
end
