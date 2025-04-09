module NetworkDynamicsDataFramesExt

using NetworkDynamics: NetworkDynamics, Network
using DataFrames: DataFrames, DataFrame
using OrderedCollections: OrderedDict

function _df_from_list(list, f)
    df = DataFrame()
    for (i, v) in list
        row = OrderedDict{Symbol,Any}()
        row[:idx] = i
        for sym in f(v)
            if NetworkDynamics.has_default_or_init(v, sym)
                row[sym] = NetworkDynamics.get_default_or_init(v, sym)
            else
                row[sym] = missing
            end
        end
        push!(df, NamedTuple(row); cols=:union)
    end
    df
end
function _df_from_pair(list, pair)
    key, f = pair
    df = DataFrame()
    for (i, v) in list
        row = OrderedDict{Symbol,Any}(:idx => i, key => f(v))
        push!(df, NamedTuple(row); cols=:union)
    end
    df
end

"""
    describe_vertices(nw::Network, extras...; parameters=true, states=true, batch=nothing)

Creates a DataFrame containing information about the vertices in a Network.

# Arguments
- `nw::Network`: The network to describe
- `extras...`: Additional pairs of (key, function) to include as columns,
   where the function gets the [`VertexModel`](@ref NetworkDynamics.VertexModel) as its only parameter
   to extract a custom metadata field for example..
- `parameters=true`: Whether to include parameter values
- `states=true`: Whether to include state values
- `batch=nothing`: Optionally filter by specific batches

# Returns
A DataFrame with columns for vertex indices, names, batch numbers, and any parameter/state values.
"""
function NetworkDynamics.describe_vertices(nw::Network, extras...; parameters=true, states=true, batch=nothing)
    pairs = enumerate(nw.im.vertexm);
    batches = map(idx -> findfirst(batch -> idx ∈ batch.indices, nw.vertexbatches), first.(pairs))

    if !isnothing(batch)
        idxs = findall(b -> b ∈ batch, batches)
        pairs = collect(pairs)[idxs]
        batch = collect(batches)[idxs]
    end
    isempty(pairs) && return DataFrame()

    basedf = DataFrame(
        idx = first.(pairs),
        name = map(v->last(v).name, pairs),
        batch = map(idx -> findfirst(batch -> idx ∈ batch.indices, nw.vertexbatches), first.(pairs)),
    )

    dfs = [basedf,]
    if parameters
        push!(dfs, _df_from_list(pairs, NetworkDynamics.psym))
    end
    if states
        push!(dfs, _df_from_list(pairs, NetworkDynamics.sym))
    end
    for p in extras
        push!(dfs, _df_from_pair(pairs, p))
    end

    foldl((a,b) -> DataFrames.leftjoin(a,b; on=:idx), dfs)
end

"""
    describe_edges(nw::Network, extras...; parameters=true, states=true, batch=nothing)

Creates a DataFrame containing information about the edges in a Network.

# Arguments
- `nw::Network`: The network to describe
- `extras...`: Additional pairs of (key, function) to include as columns,
   where the function gets the [`EdgeModel`](@ref NetworkDynamics.EdgeModel) as its only parameter
   to extract a custom metadata field for example..
- `parameters=true`: Whether to include parameter values
- `states=true`: Whether to include state values
- `batch=nothing`: Optionally filter by specific batches

# Returns
A DataFrame with columns for edge indices, source-destination pairs, names, batch numbers, and any parameter/state values.
"""
function NetworkDynamics.describe_edges(nw::Network, extras...; parameters=true, states=true, batch=nothing)
    pairs = enumerate(nw.im.edgem);
    batches = map(idx -> findfirst(batch -> idx ∈ batch.indices, nw.layer.edgebatches), first.(pairs))

    if !isnothing(batch)
        idxs = findall(b -> b ∈ batch, batches)
        pairs = collect(pairs)[idxs]
        batch = collect(batches)[idxs]
    end
    isempty(pairs) && return DataFrame()

    basedf = DataFrame(
        idx = first.(pairs),
        srcdst = map(idx -> nw.im.edgevec[idx].src => nw.im.edgevec[idx].dst, first.(pairs)),
        name = map(v->last(v).name, pairs),
        batch = map(idx -> findfirst(batch -> idx ∈ batch.indices, nw.layer.edgebatches), first.(pairs)),
    )

    dfs = [basedf,]
    if parameters
        push!(dfs, _df_from_list(pairs, NetworkDynamics.psym))
    end
    if states
        push!(dfs, _df_from_list(pairs, NetworkDynamics.sym))
    end
    for p in extras
        push!(dfs, _df_from_pair(pairs, p))
    end

    foldl((a,b) -> DataFrames.leftjoin(a,b; on=:idx), dfs)
end

end
