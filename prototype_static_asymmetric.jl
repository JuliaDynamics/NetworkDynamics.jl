using NetworkDynamics
using LightGraphs



## This holds the layer dependent structural information

# constant and mutable information should be separate!!

# StaticInteractions should subtype AbstractInteractions

struct LayerStruct
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

    #symmetric::Symbol # :a, :asymmetric, :s, :symmetric, :u, :unsymmetric
end
function LayerStruct(graph, interact!, aggregate!)
    num_v = nv(graph)
    num_e = ne(graph)
    len   = interact!.dim * num_e

    s_e = [src(e) for e in edges(graph)]
    d_e = [dst(e) for e in edges(graph)]

    v_idx_s = create_v_idxs(s_e, num_v) # edges that have their source at vertex i
    v_idx_d = create_v_idxs(d_e, num_v) # edges that have their destination at vertex i
    #
    # if interact!.symmetric == :u # make sure this complies with the data
    #     v_idx_d .+= len # we will duplicate each edge variable
    # end

    v_idx = Vector{Array{Int64,1}}(undef, num_v)
    for i in 1:num_v
        v_idx[i] = [v_idx_s[i]; v_idx_d[i]] # problems with self-loops?
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

mutable struct LayerData # {T2, T3, G, Ti} figure out typing later
    interactions::Array{Float64,1} # ::Ti
    aggregates::Array{Float64,1}
    v::Array{Float64,1}
end
function LayerData(NL::LayerStruct) # save symmetry = :u in interact!?
    interactions = zeros(NL.num_e * 2)
    if NL.dim > 1
        aggregates = zeros(NL.dim, NL.num_v) # TBD types, should aggregations be able to return vectors?
    else
        aggregates = zeros(NL.num_v)
    end
    v = zeros(NL.num_v)
    LayerData(interactions, aggregates, v)
end
function create_v_idxs(edge_idx, num_v)::Vector{Array{Int64,1}}
    v_idx = Vector{Array{Int64,1}}(undef, num_v)
    for i in 1:num_v
        v_idx[i] = findall( x -> (x==i), edge_idx)
    end
    v_idx
end



### GraphStructure

const Idx = UnitRange{Int}
struct NetworkStruct
    num_v::Int
    v_dims::Array{Int, 1}
    # v_syms::Array{Symbol, 1} tbd
    dim_v::Int
    v_offs::Array{Int, 1}
    v_idx::Array{Idx, 1}
    layers::Array{LayerStruct,1}
    # add layer offsets (once we get to ODEEdges)
end
function NetworkStruct(g, v_dims, layers)# syms later, v_syms)
    num_v = nv(g)
    v_offs = create_offsets(v_dims)
    v_idx = create_idxs(v_offs, v_dims)

    NetworkStruct(
    num_v,
    v_dims,
    #v_syms,
    sum(v_dims),
    v_offs,
    v_idx,
    layers)
end



### nd_ODE_Static

@Base.kwdef struct proto_ODE_Static{G, T1}
    vertices!::T1 # might be an array
    graph::G
    network_structure::NetworkStruct
    layer_data::Array{LayerData, 1}
    parallel::Bool # enables multithreading for the core loop
end


@Base.kwdef struct StaticInteraction{T}
    f!::T # (e, v_s, v_d, p, t) -> nothing
    dim::Int # number of dimensions of x
    sym=[:e for i in 1:dim] # Symbols for the dimensions
    symmetric=:u
end


@inline function interact!(iarr, v_in, v_out, p ,t)
    @inbounds iarr[1] = v_in - v_out
    nothing
end

si! = StaticInteraction(f! = interact!, dim = 1, symmetric = :a)

function aggregate!(interactions; x = nothing, p = nothing, t = nothing) # actually we can let this depend on x, p ,t as well and maybe come up with a new array type if we have to
  #println(agg, interaction)

end

function node!(dv, v, p, t, inputs)
    dv[1] = -inputs
    nothing
end

using OrdinaryDiffEq, Plots, BenchmarkTools

vertex! = ODEVertex(f! = node!, dim = 1)
N = 1000 # number of nodes
k = 20 # average degree
graph = barabasi_albert(N, k)

NL = LayerStruct(graph, si!, aggregate!)
LD = LayerData(NL)
NS = NetworkStruct(graph, [1 for v in 1:nv(graph)], [NL])

function (d::proto_ODE_Static)(dx, x, p, t)
    for (ls, ld) in zip(d.network_structure.layers, d.layer_data)
        for i in 1:ls.num_e
            @inbounds @views  ld.interactions[i] = x[ls.s_e[i]] - x[ls.d_e[i]]          #@inbounds ls.interact!.f!(view(ld.interactions,i),1., x[ls.d_e[i]], p, t) ## will need vertexdata and aggregationdata for convenient indexing,
            @inbounds @views ld.interactions[i + ls.num_e] = x[ls.d_e[i]] - x[ls.s_e[i]]
            #ls.interact!.f!(view(ld.interactions,i + ls.num_e), x[ls.d_e[i]], x[ls.s_e[i]], p, t)
        end
        for i in 1:d.network_structure.num_v

            @inbounds @views d.layer_data[1].aggregates[i] = sum(ld.interactions[ls.v_idx_s[i]]) + sum( ld.interactions[ls.v_idx_d[i] .+ ls.num_e])

            #@inbounds @views musum!(d.layer_data[1].aggregates[i], ld.interactions[ls.v_idx_s[i]]) #+ sum( ld.interactions[ls.v_idx_d[i] .+ ls.num_e])

            @inbounds @views  maybe_idx(d.vertices!,i).f!(
                view(dx, i),
                ld.v[i],
                p_v_idx(p, i),
                t,
                d.layer_data[1].aggregates[i])
       end
    end
    nothing
end

@btime @inbounds @views  sum($LD.interactions[$NL.v_idx_s[1]])

ode = proto_ODE_Static(vertex!, graph, NS, [LD], false)

x0 = randn(N)
dx0 = randn(N)
ode(x0, dx0, nothing, 0.)
@btime ode($x0, $dx0, nothing, 0.)

x0 = randn(N)

ode_prob = ODEProblem(ode, x0, (0., 5.))
sol = solve(ode_prob, Tsit5())

plot(sol)

display(@benchmark $ode($x0, $dx0, nothing, 0.))

#@btime solve(ode_prob, Tsit5())
