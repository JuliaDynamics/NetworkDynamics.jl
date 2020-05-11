@inline function prep_gd(dx::T, x::T, gd::GraphData{T, T}, gs) where T
    gd.v_array = x # Does not allocate
    gd
end

@inline function prep_gd(dx, x, gd, gs)
    e_array = similar(dx, gs.dim_e)
    GraphData(x, e_array, gs)
end



@Base.kwdef struct nd_ODE_Static # {G, Tv, Ti, T1, T2, T3} figure out typing later
    vertices!
    graph_structure::GraphStruct
    graph_data::GraphData{Tv, Te}
    layers::Array{Layer, 1}
    parallel::Bool # enables multithreading for the core loop
end




function (d::nd_ODE_Static)(dx, x, p, t)
    gd = prep_gd(dx, x, d.graph_data, d.graph_structure)

    for layer in d.layers

    @nd_threads d.parallel for (i, j) in 1:layer.graph_structure.edges
        layer.interact(gd.inter[i, j], gd.vertex[i], gd.vertex[j], p, t)
    end
    @nd_threads d.parallel for i in 1:d.graph_structure.vertices
        maybe_idx(d.vertices!,i).f!(view(dx,d.graph_structure.v_idx[i]), gd.v[i], p_v_idx(p, i), t, layer.aggregate(gd.inter[i, :]) )
    end

    nothing
end

## This holds the global information (nodes)
"""
Create offsets for stacked array of dimensions dims
"""
function create_offsets(dims; counter=0)::Array{Int, 1}
    offs = [1 for dim in dims]
    for (i, dim) in enumerate(dims)
        offs[i] = counter
        counter += dim
    end
    offs
end

"""
Create indexes for stacked array of dimensions dims using the offsets offs
"""
function create_idxs(offs, dims)::Array{Idx, 1}
    idxs = [1+off:off+dim for (off, dim) in zip(offs, dims)]
end

struct GraphStruct
    num_v::Int
    num_layers::Int
    v_dims::Array{Int, 1}

    v_syms::Array{Symbol, 1}

    dim_v::Int # sum(v_dims)

    v_offs::Array{Int, 1}
    v_idx::Array{Idx, 1}

    # Add some information about layers here? Maybe the layer structs?
end
function GraphStruct(g, v_dims, v_syms)
    num_v = nv(g)
    v_offs = create_offsets(v_dims)

    v_idx = create_idxs(v_offs, v_dims)

    GraphStruct(
    num_v,
    v_dims,
    v_syms,
    sum(v_dims),
    v_offs,
    v_idx)
end


## This holds the layer dependent structural information

# constant and mutable information should be separate!!


struct NetworkLayer # {T2, T3, G, Ti} figure out typing later
    interact! # ::T2
    aggregate # ::T3
    layer_structure::LayerStruct
    interactions # ::Ti
end
function NetworkLayer(interact!, aggregate, graph) # save symmetry = :u in interact!?
    layer_structure = LayerStruct(graph, interact!.dim, interact!.sym)
    if interact!.symmetry in (:u, :unsymmetric)
        interactions = zeros{T, 2}{interact!.dim, 2 * layer_structure.num_e)
    else
        interactions = zeros{T, 2}{interact!.dim, layer_structure.num_e)
    end
    Layer(
      interact!,
      aggregate,
      layer_structure,
      interactions)
end

struct LayerStruct
    num_e::Int
    e_dim::Int # constant in layer
    len::Int # e_dim * num_e
    e_syms::Array{Symbol, 1}

    s_e::Array{Int, 1}
    d_e::Array{Int, 1}

    v_idx_s::Array{Array{Int, 1},1}
    v_idx_d::Array{Array{Int, 1},1}
    v_idx::Array{Array{Int, 1},1}

    symmetric::Symbol # :a, :asymmetric, :s, :symmetric, :u, :unsymmetric
end
function LayerStruct(g, e_dim, e_syms)
    num_e = ne(g)

    s_e = [src(e) for e in edges(g)]
    d_e = [dst(e) for e in edges(g)]

    v_idx_s = create_v_idxs(s_e, num_v)
    v_idx_d = create_v_idxs(d_e, num_v)
    v_idx = Vector{Array{Int64,1}}(undef, num_v)
    for i in 1:num_v
        v_idx[i] = [v_idx_s; v_idx_d]
    end

    LayerStruct(
    num_e,
    e_dim,
    e_dim * num_e,
    e_syms,
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

## Layer independent data

mutable struct GraphData{Tv}
    v_array::Tv
    v::Array{VertexData{GraphData{Tv}, Tv}, 1}
    function GraphData{Tv, Te}(v_array::Tv, gs::GraphStruct) where {Tv}
        gd = new{Tv, Te}(v_array, )
        gd.v = [VertexData{GraphData{Tv, Te}, Tv}(gd, offset, dim) for (offset,dim) in zip(gs.v_offs, gs.v_dims)]
    end
end
