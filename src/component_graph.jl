"""
    ComponentGraph(nv, edges::Vector{SimpleEdge{Int}})

Directed multigraph holding the topology of a [`Network`](@ref). Unlike `SimpleGraph` and
`SimpleDiGraph` it keeps parallel edges and stores edges in the given order, so the position in
`edges(g)` is the edge index in the network.

Degrees count edges with multiplicity, neighbor lists contain every neighbor only once.
"""
struct ComponentGraph <: AbstractGraph{Int}
    nv::Int
    edges::Vector{SimpleEdge{Int}}
    # lookup tables for the per-vertex queries, derived from `edges` once
    fadjlist::Vector{Vector{Int}} # sorted unique out-neighbors
    badjlist::Vector{Vector{Int}} # sorted unique in-neighbors
    outdeg::Vector{Int}           # with multiplicity
    indeg::Vector{Int}
    function ComponentGraph(nv::Integer, edges::AbstractVector{<:SimpleEdge})
        nv ≥ 0 || throw(ArgumentError("Number of vertices must be non-negative, got $nv"))
        edges = collect(SimpleEdge{Int}, edges)
        fadjlist = [Int[] for _ in 1:nv]
        badjlist = [Int[] for _ in 1:nv]
        outdeg = zeros(Int, nv)
        indeg = zeros(Int, nv)
        for e in edges
            if !(1 ≤ e.src ≤ nv && 1 ≤ e.dst ≤ nv)
                throw(ArgumentError("Edge $(e.src) => $(e.dst) references a vertex outside of 1:$nv"))
            end
            push!(fadjlist[e.src], e.dst)
            push!(badjlist[e.dst], e.src)
            outdeg[e.src] += 1
            indeg[e.dst] += 1
        end
        foreach(l -> unique!(sort!(l)), fadjlist)
        foreach(l -> unique!(sort!(l)), badjlist)
        new(nv, edges, fadjlist, badjlist, outdeg, indeg)
    end
end

Graphs.nv(g::ComponentGraph) = g.nv
Graphs.ne(g::ComponentGraph) = length(g.edges)
Graphs.vertices(g::ComponentGraph) = Base.OneTo(g.nv)
Graphs.edges(g::ComponentGraph) = g.edges
Graphs.edgetype(::ComponentGraph) = SimpleEdge{Int}
Base.eltype(::ComponentGraph) = Int
Base.eltype(::Type{ComponentGraph}) = Int
Base.zero(::Type{ComponentGraph}) = ComponentGraph(0, SimpleEdge{Int}[])
Graphs.has_vertex(g::ComponentGraph, v::Integer) = 1 ≤ v ≤ g.nv
Graphs.is_directed(::Type{ComponentGraph}) = true
Graphs.is_directed(::ComponentGraph) = true

function Graphs.has_edge(g::ComponentGraph, s::Integer, d::Integer)
    Graphs.has_vertex(g, s) && Graphs.has_vertex(g, d) && insorted(d, g.fadjlist[s])
end

Graphs.outneighbors(g::ComponentGraph, v::Integer) = g.fadjlist[v]
Graphs.inneighbors(g::ComponentGraph, v::Integer) = g.badjlist[v]

# degrees count parallel edges, so they come from their own tables and not the neighbor lists
Graphs.outdegree(g::ComponentGraph, v::Integer) = g.outdeg[v]
Graphs.indegree(g::ComponentGraph, v::Integer) = g.indeg[v]
Graphs.outdegree(g::ComponentGraph) = copy(g.outdeg)
Graphs.indegree(g::ComponentGraph) = copy(g.indeg)
Graphs.degree(g::ComponentGraph, v::Integer) = g.indeg[v] + g.outdeg[v]
Graphs.degree(g::ComponentGraph) = g.indeg .+ g.outdeg

Base.:(==)(a::ComponentGraph, b::ComponentGraph) = a.nv == b.nv && a.edges == b.edges
Base.hash(g::ComponentGraph, h::UInt) = hash(g.edges, hash(g.nv, hash(ComponentGraph, h)))

Base.show(io::IO, g::ComponentGraph) = print(io, "ComponentGraph($(g.nv) vertices, $(ne(g)) edges)")

"""
    legacy_graph_type(g::ComponentGraph) -> :simple | :digraph | :none

Which of `SimpleGraph` or `SimpleDiGraph` can represent `g` without losing an edge. Returns
`:none` if some ordered endpoint pair appears more than once.
"""
function legacy_graph_type(g::ComponentGraph)
    allunique(g.edges) || return :none
    all(e -> e.src < e.dst, g.edges) ? :simple : :digraph
end

"""
    edge_multiplicity(edgevec) -> Vector{Tuple{Int,Int}}

For every edge return `(k, n)`: it is the k-th of n edges with the same ordered endpoints.
"""
function edge_multiplicity(edgevec)
    total = Dict{eltype(edgevec),Int}()
    for e in edgevec
        total[e] = get(total, e, 0) + 1
    end
    seen = Dict{eltype(edgevec),Int}()
    map(edgevec) do e
        seen[e] = get(seen, e, 0) + 1
        (seen[e], total[e])
    end
end
