OEdge(src, dst) = Edge(min(src, dst), max(src, dst))

LightGraphs.edges(g, u) = (OEdge(u, v) for v in neighbors(g, u))

function fan(g, colors, X, f; last = nothing)
    vec = [f]
    colored_neighbors = filter(i -> haskey(colors, OEdge(X, i)), neighbors(g, X))
    if last !== nothing && last != f
        idx = findfirst(isequal(last), colored_neighbors)
        deleteat!(colored_neighbors, idx)
        push!(colored_neighbors, last)
    end

    @show colored_neighbors
    for i in colored_neighbors
        c = colors[OEdge(X, i)]
        if isfree(g, colors, vec[end], c)
            push!(vec, i)
        else
            println("could not append to fan")
        end
    end

    @assert !haskey(colors, OEdge(X, vec[1])) "first element needs to be colorless"
    for i in 2:lastindex(vec)
        @assert haskey(colors, OEdge(X, i)) "element $i ($X->$(vec[i])) needs to have color"
    end

    return vec
end

function isfree(g, colors, u, c)
    for e in edges(g, u)
        haskey(colors, e) || continue
        c == colors[e] && return false
    end
    return true
end

function pickfree(g, colors, u; skip=-1)
    c = 1
    while c==skip || !isfree(g, colors, u, c)
        c = c + 1
    end
    return c
end

function invert_ucdpath!(g, colors, u, c, d)
    g2 = SimpleGraph(nv(g))
    for (e, col) in colors
        if col ∈ (c, d)
            add_edge!(g2, e)
        end
    end
    for component in connected_components(g2)
        if u ∈ component
            for e in edges(g2[component])
                oldc = colors[e]
                colors[e] = oldc == c ? d : c
            end
            return nothing
        end
    end
end

function find_w(g, colors, X, F, d)
    idx = 1
    w = F[idx]
    F′ = fan(g, colors, X, F[1], last=w)
    while !isfree(g, colors, w, d) || F′[end] !== w
        idx = idx + 1
        idx > length(F) && error("Could not find w")
        w = F[idx]
        F′ = fan(g, colors, X, F[1], last=w)
    end
    return w, F′
end

function rotate!(F, X, colors)
    for i in 1:lastindex(F)-1
        colors[OEdge(X, i)] = colors[OEdge(X, i+1)]
    end
    delete!(colors, OEdge(X, lastindex(F)))
end

function color_edges(g::SimpleGraph)
    colors = Dict{Edge, Int}()

    for e in edges(g)
        println("process edge $e")
        X, f = e.src, e.dst
        F = fan(g, colors, X, f)
        c = pickfree(g, colors, X)
        d = pickfree(g, colors, F[end])
        @show F c d
        invert_ucdpath!(g, colors, X, c, d)
        w, F′ = find_w(g, colors, X, F, d)
        @assert F′[end] == w
        rotate!(F′, X, colors)
        colors[e] = d
    end
    return colors
end


#=
algorithm Misra & Gries edge coloring algorithm is
    input: A graph G.
    output: A proper coloring c of the edges of G.

    Let U := E(G)

    while U ≠ ∅ do
        Let (u, v) be any edge in U.
        Let F[1:k] be a maximal fan of u starting at F[1] = v.
        Let c be a color that is free on u and d be a color that is free on F[k].
        Invert the cdu path
        Let w ∈ V(G) be such that w ∈ F, F' = [F[1]...w] is a fan and d is free on w.
        Rotate F' and set c(u, w) = d.
        U := U − {(u, v)}
    end while
=#
