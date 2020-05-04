using NetworkDynamics
using LightGraphs

@Base.kwdef struct StaticInteraction{T}
    f!::T # (e, v_s, v_d, p, t) -> nothing
    dim::Int # number of dimensions of x
    sym=[:e for i in 1:dim] # Symbols for the dimensions
    symmetric=:u
end

function interact!(v_in, v_out, p ,t)
    sin.(v_in .- v_out)
end

si! = StaticInteraction(f! = interact!, dim = 1, symmetric = :a)

function aggregate!(agg, interaction; x = nothing, p = nothing, t = nothing) # actually we can let this depend on x, p ,t as well and maybe come up with a new array type if we have to
  #println(agg, interaction)
  agg .+= interaction[1]
  interaction
end

function node!(dv, v, p, t, inputs)
    dv .= -inputs
end

vertex! = ODEVertex(f! = node!, dim = 1)
N = 100 # number of nodes
k = 20  # average degree
graph = barabasi_albert(N, k)


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

    v_idx_s = create_v_idxs(s_e, num_v)
    v_idx_d = create_v_idxs(d_e, num_v)

    if interact!.symmetric == :u # make sure this complies with the data
        v_idx_d .+= len # we will duplicate each edge variable
    end

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

struct LayerData # {T2, T3, G, Ti} figure out typing later
    interactions # ::Ti
    aggregates
    temp # intermediate saving, this should really be temporary
end
function LayerData(NL::LayerStruct) # save symmetry = :u in interact!?
    interactions = nothing
    if !(typeof(NL.interact!) <: StaticInteraction)
        # Let j be the index of an edge. For multidimensional edges interactions[j] should be an array. This seems to favour introducing a new data type with custom indexing.
        interactions = something
    end
    if NL.dim > 1
        aggregates = zeros(NL.dim, NL.num_v) # TBD types, should aggregations be able to return vectors?
    else
        aggregates = zeros(NL.num_v)
    end
    temp = zeros(NL.dim)
    LayerData(interactions, aggregates,temp)
end
function create_v_idxs(edge_idx, num_v)::Vector{Array{Int64,1}}
    v_idx = Vector{Array{Int64,1}}(undef, num_v)
    for i in 1:num_v
        v_idx[i] = findall( x -> (x==i), edge_idx)
    end
    v_idx
end

NL = LayerStruct(graph, si!, aggregate!)
LD = LayerData(NL)

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

NS = NetworkStruct(graph, [1 for v in 1:nv(graph)], [NL])

### nd_ODE_Static

@Base.kwdef struct proto_ODE_Static{G, T1}
    vertices!::T1 # might be an array
    graph::G
    network_structure::NetworkStruct
    layer_data::Array{LayerData, 1}
    parallel::Bool # enables multithreading for the core loop
end


function (d::proto_ODE_Static)(dx, x, p, t)
    for (ls, ld) in zip(d.network_structure.layers, d.layer_data)
        ld.aggregates .*= 0 # this presumes aggregation is addition
        if ls.interact!.symmetric == :s
            for i in 1:ls.num_e #adjust syntax
                ld.temp .= ls.interact!.f!(x[ls.s_e[i]], x[ls.d_e[i]], p, t) ## will need vertexdata and aggregationdata for convenient indexing,
                ls.aggregate!(ld.aggregates[ls.s_e[i]], ld.temp)
                ls.aggregate!(ld.aggregates[ls.d_e[i]], ld.temp)
            end
        elseif ls.interact!.symmetric == :a
          for i in 1:ls.num_e
         #  will need vertexdata and aggregationdata for convenient indexing,
         # confusing functional tricks?
              ls.aggregate!(view(ld.aggregates, ls.d_e[i]), -1. .*
                ls.aggregate!(view(ld.aggregates, ls.s_e[i]),
                  ls.interact!.f!(x[ls.s_e[i]], x[ls.d_e[i]], p, t)))
         # easier syntax, more allocations
               # ld.temp .= ls.interact!.f!(x[ls.s_e[i]], x[ls.d_e[i]], p, t) ## will need vertexdata and aggregationdata for convenient indexing,
                # ls.aggregate!(view(ld.aggregates, ls.s_e[i]), ld.temp)
                # ls.aggregate!(view(ld.aggregates, ls.d_e[i]), -1. .* ld''.temp)
            end
       else
           for i in 1:ls.num_e #adjust syntax
               ls.aggregate!(ld.aggregates[ls.s_e[i]],
                 ls.interact!.f!(x[ls.s_e[i]], x[ls.d_e[i]], p, t))
               ls.aggregate!(ld.aggregates[ls.d_e[i]],
                 ls.interact!.f!(x[ls.d_e[i]], x[ls.s_e[i]], p, t))
            end
       end
       for i in 1:d.network_structure.num_v
           maybe_idx(d.vertices!,i).f!(view(dx, d.network_structure.v_idx[i]), x[d.network_structure.v_idx[i]], p_v_idx(p, i), t, ld.aggregates[i])
       end
    end
    nothing
end

ode = proto_ODE_Static(vertex!, graph, NS, [LD], false)

x0 = randn(N)
x = copy(x0)
dx0 = randn(N)
dx = copy(dx0)
ode(x0, dx0, nothing, 0.)
x.-x0

dx .- dx0

x0 = randn(N)
using OrdinaryDiffEq, Plots
ode_prob = ODEProblem(ode, x0, (0., 4.))
sol = solve(ode_prob, Tsit5())

plot(sol)

using BenchmarkTools



@benchmark ode(x0, randn(100), nothing, 0.)

@btime solve(ode_prob, Tsit5())
