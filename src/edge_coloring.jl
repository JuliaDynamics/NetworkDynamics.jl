using LightGraphs.SimpleGraphs

struct ColoredGraph
    g::SimpleGraph{Int}
    colors::Dict{SimpleEdge{Int}, Int}
end

ColoredGraph(g) = ColoredGraph(g, Dict{SimpleEdge{Int}, Int}())

OEdge(p::Pair) = OEdge(p.first, p.second)
OEdge(src::Int, dst::Int) = SimpleEdge(min(src, dst), max(src, dst))

LightGraphs.neighbors(g::ColoredGraph, u) = neighbors(g.g, u)
LightGraphs.nv(g::ColoredGraph) = nv(g.g)
LightGraphs.edges(g::ColoredGraph) = edges(g.g)

connected_edges(g::ColoredGraph, u) = (OEdge(u, v) for v in neighbors(g, u))

hascolor(g::ColoredGraph, e::SimpleEdge) = haskey(g.colors, e)
hascolor(g::ColoredGraph, p::Pair) = hascolor(g, OEdge(p))
color(g::ColoredGraph, e::SimpleEdge) = g.colors[e]
color(g::ColoredGraph, p::Pair) = color(g, OEdge(p))
setcolor!(g::ColoredGraph, e::SimpleEdge, c::Int) = g.colors[e] = c
setcolor!(g::ColoredGraph, p::Pair, c::Int) = setcolor!(g, OEdge(p), c)

colors(g::ColoredGraph) = g.colors

function isvalid(g::ColoredGraph)
    valid = true
    for v in 1:nv(g)
        colors = [color(g, e) for e in connected_edges(g, v) if hascolor(g, e)]
        unique = allunique(colors)
        if !unique
            println("Found colors $colors at vertex $v")
        end
        valid = valid && unique
    end
    return valid
end

function isfree(g::ColoredGraph, u::Int, c::Int)
    for e in connected_edges(g, u)
        hascolor(g, e) || continue
        c == color(g, e) && return false
    end
    return true
end

function pickfree(g::ColoredGraph, u)
    c = 1
    while !isfree(g, u, c)
        c = c + 1
    end
    return c
end

function pickfree(g::ColoredGraph, u1, u2)
    c = 1
    while !isfree(g, u1, c) || !isfree(g, u2, c)
        c = c + 1
    end
    return c
end

struct Fan
    g::ColoredGraph
    X::Int
    vec::Vector{Int}
end

Base.getindex(f::Fan, i) = f.vec[i]
Base.firstindex(f::Fan) = firstindex(f.vec)
Base.lastindex(f::Fan) = lastindex(f.vec)
Base.eachindex(f::Fan) = eachindex(f.vec)
Base.length(f::Fan) = length(f.vec)
Base.iterate(f::Fan) = iterate(f.vec)
Base.iterate(f::Fan, state) = iterate(f.vec, state)

function rotate!(F::Fan)
    for i in 1:lastindex(F)-1
        e = OEdge(F.X, F[i])
        newc = color(F.g, F.X => F[i+1])
        setcolor!(F.g, e, newc)
    end
    delete!(F.g.colors, OEdge(F.X, F[end]))
end

function maximalfan(g::ColoredGraph, X::Int, f::Int)
    vec = [f]
    colored_neighbors = filter(i -> hascolor(g, X=>i), neighbors(g, X))

    # isempty(colored_neighbors) && return Fan(g, X, vec)

    # if the last keyword is given then move element to the last place
    # if last !== nothing
    #     idx = findfirst(isequal(last), colored_neighbors)
    #     oldend = colored_neighbors[end]
    #     colored_neighbors[end] = last
    #     colored_neighbors[idx] = oldend
    # end

    # for i in colored_neighbors
    #     c = color(g, X=>i)
    #     if isfree(g, vec[end], c)
    #         push!(vec, i)
    #     else
    #         error("Could not append $X=>$i of color $c to fan (color not free on $(vec[end]))")
    #     end
    # end

    @assert !hascolor(g, X=>vec[1]) "first element needs to be colorless"
    for i in 2:lastindex(vec)
        @assert hascolor(g, X=>vec[i]) "element $i ($X->$(vec[i])) needs to have color"
    end

    return Fan(g, X, vec)
end

function traverse(g::ColoredGraph, pos, c)
    for n in neighbors(g, pos)
        hascolor(g, pos=>n) && color(g, pos=>n) == c && return n
    end
    return nothing
end

function invert_cdpath!(g::ColoredGraph, X, c, d)
    # we've chosen c in such a way
    @assert isfree(g, X, c)

    invert(i) = i==c ? d : c

    pos = X
    lookfor = d
    next = traverse(g, pos, lookfor)

    while next !== nothing
        edge = OEdge(pos, next)

        pos = next
        lookfor = invert(lookfor)
        next = traverse(g, pos, lookfor)

        setcolor!(g, edge, lookfor)
    end
end

"""
- d free at X
w satisfies:
  - w ∈ f..l (in the fan F)
  - d is free at w
  - f..w is Fan of X
"""
function find_w(g::ColoredGraph, X, F, d)
    @assert isfree(g, X, d) "d=$d not free on X=$X"

    # for i in lastindex(F):-1:firstindex(F)
    #     w = F[i]
    for (i, w) in enumerate(F)
        isfree(g, w, d) || continue
        # F′ = Fan(g, X, F[1], last=w)
        F′ = Fan(g, X, F[1:i])
        @assert F′[end] == w
        return w, F′
    end
    @warn "$d is not free at any of $(F.vec)"
end

function color_edges(g::SimpleGraph)
    gc = ColoredGraph(g)

    for e in edges(gc)
        @info "process edge $e" gc.colors
        X, f = e.src, e.dst
        F = Fan(gc, X, f)
        c = pickfree(gc, X)
        d = pickfree(gc, F[end])
        println("  Fan = $X->$(F.vec)")
        println("  c = $(c) free at $X")
        println("  d = $(d) free at $(F[end])")
        @info "Invert path c=$c d=$d around $X" gc.colors
        invert_cdpath!(gc, X, c, d)
        @info "after invert" gc.colors
        w, F′ = find_w(gc, X, F, d)
        println("  w = $(w) with Fan X->$(F′.vec)")
        rotate!(F′)
        println("  set color $X=>$w to $d")
        setcolor!(gc, OEdge(X, w), d)
        isvalid(gc) || @warn "oh no"
    end
    return gc
end

function color_edges_greedy(g::SimpleGraph)
    gc = ColoredGraph(g)
    for e in edges(gc)
        c = pickfree(gc, e.src, e.dst)
        setcolor!(gc, e, c)
    end
    @assert isvalid(gc) "Did not find valid edge coloring 😱"
    return gc
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
