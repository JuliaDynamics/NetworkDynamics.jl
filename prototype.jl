using NetworkDynamics

@Base.kwdef struct StaticInteraction{T}
    f!::T # (e, v_s, v_d, p, t) -> nothing
    dim::Int # number of dimensions of x
    sym=[:e for i in 1:dim] # Symbols for the dimensions
    symmetric=:u
end

function interact!(inter, v_in, v_out, p ,t)
    inter = v_in - v_out
    nothing
end

si! = StaticInteraction(f! = interact!, dim = 1, symmetric = :a)

function aggregate!(agg, interactions, p, t) # actually we can let this depend on p as well and maybe come up with a new array type if we have to
  agg = sum(interactions)
  nothing
end

function node!(dv, v, p, t, inputs)
    dv .= -inputs
end

vertex! = ODEVertex(f! = node!, dim = 1)

## InteractionData

import Base.getindex, Base.setindex!, Base.length, Base.IndexStyle, Base.size, Base.eltype

struct InteractionData{T} <: AbstractArray{T, 1}
    data::AbstractArray{T,1}
    len::Int
    dim::Int
end
@inline Base.@propagate_inbounds function getindex(i_dat::InteractionData, idx<:Integer) #think about how slicing can be achieved
    # i_dat[j] should return a n_neighbors_j x i_dat dim array
    i_dat.data[range(idx-1 *i_dat.dim + 1, length=i_dat.dim)]
end

@inline Base.@propagate_inbounds function getindex(i_dat::InteractionData, idx) #think about how slicing can be achieved
    # i_dat[j] should return a n_neighbors_j x i_dat dim array

    i_dat.data[range.((idx.-1).*i_dat.dim .+ 1, length=i_dat.dim)]
end

@inline Base.@propagate_inbounds function setindex!(i_dat::InteractionData, x, idx)
    i_dat.data[(idx.-1).*i_dat.dim .+ 1:idx.*i_dat.dim] = x
    nothing
end

@inline function Base.length(i_dat::InteractionData)
    i_dat.len
end

@inline function Base.size(i_dat::InteractionData)
    (i_dat.len, )
end

@inline function Base.eltype(i_dat::InteractionData{T}) where T
    T
end

Base.IndexStyle(::Type{<:InteractionData}) = IndexLinear()

## This holds the layer dependent structural information

# constant and mutable information should be separate!!

# StaticInteractions should subtype AbstractInteractions


struct LayerData # {T2, T3, G, Ti} figure out typing later
    interactions # ::Ti
    aggregates
end
function LayerData(NL::NetworkLayer) # save symmetry = :u in interact!?
    # Let j be the index of an edge. For multidimensional edges interactions[j] should be an array. This seems to favour introducing a new data type with custom indexing.
    interactions =
    aggregates = Array{Float64}(undef, NL.num_v) # TBD types, should aggregations be able to return vectors?
    LayerData(interactions, aggregates)
end

struct NetworkLayer
    # Functions
    interact!
    aggregate!

    num_v::Int
    num_e::Int
    dim::Int # constant in layer
    len::Int # e_dim * num_e MAYBE
    # syms::Array{Symbol, 1} add later

    s_e::Array{Int, 1}
    d_e::Array{Int, 1}

    v_idx_s::Array{Array{Int, 1},1}
    v_idx_d::Array{Array{Int, 1},1}
    v_idx::Array{Array{Int, 1},1}

    symmetric::Symbol # :a, :asymmetric, :s, :symmetric, :u, :unsymmetric
end
function NetworkLayer(graph, interact!, aggregate!)
    num_v = nv(graph)
    num_e = ne(graph)
    len   = interact!.dim * num_e

    s_e = [src(e) for e in edges(graph)]
    d_e = [dst(e) for e in edges(graph)]

    v_idx_s = create_v_idxs(s_e, num_v)
    v_idx_d = create_v_idxs(d_e, num_v)

    if interact!.symmetric == :u # make sure this complies with the data
        v_idx_d .+= len # we will duplicate each edge variable
    end

    v_idx = Vector{Array{Int64,1}}(undef, num_v)
    for i in 1:num_v
        v_idx[i] = [v_idx_s; v_idx_d] # problems with self-loops?
    end

    LayerStruct(
    interact!,
    aggregate!,
    num_v,
    num_e,
    interact!.dim,
    len,
    #e_syms, # This should probably be an array
    s_e,
    d_e,
    v_idx_s,
    v_idx_d,
    v_idx)
end

function create_v_idxs(edge_idx, num_v)::Vector{Array{Int64,1}}
    v_idx = Vector{Array{Int64,1}}(undef, num_v)
    for i in 1:num_v
        v_idx[i] = findall( x -> (x==i), edge_idx)
    end
    v_idx
end
