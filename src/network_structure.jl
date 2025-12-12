mutable struct IndexManager{G}
    g::G
    edgevec::Vector{SimpleEdge{Int64}}
    # positions of vertex data
    v_data::Vector{UnitRange{Int}}  # v data in flat states
    v_out::Vector{UnitRange{Int}}   # v range in output buf
    v_para::Vector{UnitRange{Int}}  # v para in flat para
    v_aggr::Vector{UnitRange{Int}}  # v input in aggbuf
    v_ext::Vector{UnitRange{Int}}   # v external inputs in ext buf
    # positions of edge data
    e_data::Vector{UnitRange{Int}}  # e data in flat states
    e_out::Vector{@NamedTuple{src::UnitRange{Int},dst::UnitRange{Int}}}   # e range in output buf
    e_para::Vector{UnitRange{Int}}  # e para in flat para
    e_gbufr::Vector{@NamedTuple{src::UnitRange{Int},dst::UnitRange{Int}}}   # e range in input buf
    e_ext::Vector{UnitRange{Int}}   # v external inputs in ext buf
    # metadata
    edepth::Int
    vdepth::Int
    lastidx_dynamic::Int
    lastidx_out::Int
    lastidx_p::Int
    lastidx_aggr::Int
    lastidx_gbuf::Int
    lastidx_extbuf::Int
    vertexm::Vector{VertexModel}
    edgem::Vector{EdgeModel}
    aliased_vertexms::IdDict{VertexModel, @NamedTuple{idxs::Vector{Int}, hash::UInt}}
    aliased_edgems::IdDict{EdgeModel, @NamedTuple{idxs::Vector{Int}, hash::UInt}}
    unique_vnames::Dict{Symbol,Int}
    unique_enames::Dict{Symbol,Int}
    function IndexManager(g, dyn_states, edepth, vdepth, vertexm, edgem; valias, ealias)
        aliased_vertexm_hashes = _aliased_hashes(VertexModel, vertexm, valias)
        aliased_edgem_hashes = _aliased_hashes(EdgeModel, edgem, ealias)
        unique_vnames = unique_mappings(getproperty.(vertexm, :name), 1:nv(g))
        unique_enames = unique_mappings(getproperty.(edgem, :name), 1:ne(g))
        new{typeof(g)}(g, collect(edges(g)),
                       (Vector{UnitRange{Int}}(undef, nv(g)) for i in 1:5)...,
                       Vector{UnitRange{Int}}(undef, ne(g)),
                       Vector{@NamedTuple{src::UnitRange{Int},dst::UnitRange{Int}}}(undef, ne(g)),
                       Vector{UnitRange{Int}}(undef, ne(g)),
                       Vector{@NamedTuple{src::UnitRange{Int},dst::UnitRange{Int}}}(undef, ne(g)),
                       Vector{UnitRange{Int}}(undef, ne(g)),
                       edepth, vdepth,
                       0, 0, 0, 0, 0, 0,
                       vertexm, edgem,
                       aliased_vertexm_hashes,
                       aliased_edgem_hashes,
                       unique_vnames,
                       unique_enames)
    end
end
function _aliased_hashes(T, cfs, aliastype)
    hashdict = IdDict{T, @NamedTuple{idxs::Vector{Int}, hash::UInt}}()
    if aliastype == :some
        ag = aliasgroups(cfs)
        for (c, idxs) in ag
            h = hash(c)
            hashdict[c] = (; idxs=idxs, hash=h)
        end
    elseif aliastype == :all
        c = first(cfs)
        hashdict[c] =(; idxs=collect(eachindex(cfs)), hash=hash(c))
    end
    hashdict
end

dim(im::IndexManager) = im.lastidx_dynamic
pdim(im::IndexManager) = im.lastidx_p


struct Network{EX<:ExecutionStyle,G,NL,VTup,MM,CT,GBT,LM,EM}
    "vertex batches of same function"
    vertexbatches::VTup
    "network layer"
    layer::NL
    "index manager"
    im::IndexManager{G}
    "lazy cache pool"
    caches::@NamedTuple{output::CT,aggregation::CT,external::CT}
    "mass matrix"
    mass_matrix::MM
    "Gather buffer provider (lazy or eager)"
    gbufprovider::GBT
    "map for loopback edge gather"
    loopbackmap::LM
    "map to gather external inputs"
    extmap::EM
    "sparsity pattern"
    jac_prototype::Ref{Union{Nothing,SparseMatrixCSC{Bool,Int}}}
    function Network(ex, vb, nl, im, caches, mm, gbufp, loopmap, extmap, jac_prototype)
        new{
            ex,typeof(im.g),typeof(nl), typeof(vb),
            typeof(mm),eltype(caches),typeof(gbufp),
            typeof(loopmap),typeof(extmap)
        }(
            vb, nl, im, caches, mm, gbufp, loopmap, extmap, jac_prototype
        )
    end
end
function ConstructionBase.constructorof(::Type{<:Network{EX}}) where {EX}
    return (args...) -> Network(EX, args...)
end

executionstyle(::Network{ex}) where {ex} = ex()
nvbatches(::Network) = length(vertexbatches)

"""
    dim(nw::Network)

Returns the number of dynamic states in the network,
corresponts to the length of the flat state vector.
"""
dim(nw::Network) = dim(nw.im)

"""
    pdim(nw::Network)

Returns the number of parameters in the network,
corresponts to the length of the flat parameter vector.
"""
pdim(nw::Network) = pdim(nw.im)
Graphs.nv(nw::Network) = nv(nw.im.g)
Graphs.ne(nw::Network) = ne(nw.im.g)
Base.broadcastable(nw::Network) = Ref(nw)

"""
    get_graph(nw::Network)

Extracts the underlying graph of the network.
"""
get_graph(nw::Network) = nw.im.g

function get_output_cache(nw::Network, T)
    if eltype(T) <: AbstractFloat && eltype(nw.caches.output.du) != eltype(T)
        throw(ArgumentError("Network caches are initialized with $(eltype(nw.caches.output.du)) \
            but is used for $(eltype(T)) data!"))
    end
    o = get_tmp(nw.caches.output, T)
    fill!(o, convert(eltype(o), NaN))
end
get_aggregation_cache(nw::Network, T) = get_tmp(nw.caches.aggregation, T)
function get_extinput_cache(nw::Network, T)
    ext = get_tmp(nw.caches.external, T)
    fill!(ext, convert(eltype(ext), NaN))
end

iscudacompatible(nw::Network) = iscudacompatible(executionstyle(nw)) && iscudacompatible(nw.layer.aggregator)

struct NetworkLayer{GT,ETup,AF}
    "graph/toplogy of layer"
    g::GT
    "edge batches with same function"
    edgebatches::ETup
    "aggregator object"
    aggregator::AF
    "depth of edge aggregation"
    edepth::Int # potential becomes range for multilayer
    "vertex dimensions visible to edges"
    vdepth::Int # potential becomes range for multilayer
end

struct ComponentBatch{T,F,G,FFT,DIM,PDIM,INDIMS,OUTDIMS,EXTDIM,IV}
    "indices contained in batch"
    indices::IV
    "internal function"
    compf::F
    "output function"
    compg::G
    ff::FFT
    "state: dimension and first index"
    statestride::BatchStride{DIM}
    "parameter: dimension and first index"
    pstride::BatchStride{PDIM}
    "inputbuf: dimension(s) and first index"
    inbufstride::BatchStride{INDIMS}
    "outputbuf: dimension(s) and first index"
    outbufstride::BatchStride{OUTDIMS}
    "external inputs: dimension and first index"
    extbufstride::BatchStride{EXTDIM}
    function ComponentBatch(dT, i, f, g, ff, ss, ps, is, os, es)
        new{dT,typeof.((f,g,ff))...,stridesT.((ss, ps, is, os, es))...,typeof(i)}(
            i, f, g, ff, ss, ps, is, os, es)
    end
end
function ConstructionBase.constructorof(::Type{<:ComponentBatch{T}}) where {T}
    return (args...) -> ComponentBatch(T, args...)
end

@inline Base.length(cb::ComponentBatch) = Base.length(cb.indices)
@inline dispatchT(::ComponentBatch{T}) where {T} = T
@inline compf(b::ComponentBatch) = b.compf
@inline compg(b::ComponentBatch) = b.compg
@inline fftype(b::ComponentBatch) = b.ff
@inline pdim(b::ComponentBatch) = b.pstride.strides
@inline extdim(b::ComponentBatch) = b.extbufstride.strides

@inline state_range(batch) = _fullrange(batch.statestride, length(batch))

@inline state_range(batch, i)     = _range(batch.statestride, i)
@inline parameter_range(batch, i) = _range(batch.pstride, i)
@inline out_range(batch, i)       = _range(batch.outbufstride, i)
@inline out_range(batch, i, j)    = _range(batch.outbufstride, i, j)
@inline in_range(batch, i)        = _range(batch.inbufstride, i)
@inline in_range(batch, i, j)     = _range(batch.inbufstride, i, j)
@inline extbuf_range(batch, i)    = _range(batch.extbufstride, i)

function register_vertices!(im::IndexManager, idxs, dim, outdim, pdim, extdim)
    for i in idxs
        im.v_data[i] = _nexturange!(im, dim)
        im.v_out[i]  = _nextoutrange!(im, outdim)
        im.v_para[i] = _nextprange!(im, pdim)
        im.v_aggr[i] = _nextaggrrange!(im, im.edepth)
        im.v_ext[i]  = _nextextrange!(im, extdim)
    end
    (;
     state = BatchStride(first(im.v_data[first(idxs)]), dim),
     p     = BatchStride(first(im.v_para[first(idxs)]), pdim),
     in    = BatchStride(first(im.v_aggr[first(idxs)]), im.edepth),
     out   = BatchStride(first(im.v_out[first(idxs)]), outdim),
     ext   = BatchStride(first(im.v_ext[first(idxs)]), extdim),
    )
end
function register_edges!(im::IndexManager, idxs, dim, outdim, pdim, extdim)
    for i in idxs
        e = im.edgevec[i]
        im.e_data[i]  = _nexturange!(im, dim)
        im.e_out[i]   = (src = _nextoutrange!(im, outdim.src),
                         dst = _nextoutrange!(im, outdim.dst))
        im.e_para[i]  = _nextprange!(im, pdim)
        im.e_gbufr[i] = (src = _nextgbufrange!(im, im.vdepth),
                         dst = _nextgbufrange!(im, im.vdepth))
        im.e_ext[i]  = _nextextrange!(im, extdim)
    end
    (;
     state = BatchStride(first(im.e_data[first(idxs)]), dim),
     p     = BatchStride(first(im.e_para[first(idxs)]), pdim),
     in    = BatchStride(first(flatrange(im.e_gbufr[first(idxs)])), (;src=im.vdepth, dst=im.vdepth)),
     out   = BatchStride(first(flatrange(im.e_out[first(idxs)])), (;src=outdim.src, dst=outdim.dst)),
     ext   = BatchStride(first(im.e_ext[first(idxs)]), extdim),
    )
end
function _nexturange!(im::IndexManager, N)
    newlast, range = _nextrange(im.lastidx_dynamic, N)
    im.lastidx_dynamic = newlast
    return range
end
function _nextoutrange!(im::IndexManager, N)
    newlast, range = _nextrange(im.lastidx_out, N)
    im.lastidx_out = newlast
    return range
end
function _nextprange!(im::IndexManager, N)
    newlast, range = _nextrange(im.lastidx_p, N)
    im.lastidx_p = newlast
    return range
end
function _nextaggrrange!(im::IndexManager, N)
    newlast, range = _nextrange(im.lastidx_aggr, N)
    im.lastidx_aggr = newlast
    return range
end
function _nextgbufrange!(im::IndexManager, N)
    newlast, range = _nextrange(im.lastidx_gbuf, N)
    im.lastidx_gbuf = newlast
    return range
end
function _nextextrange!(im::IndexManager, N)
    newlast, range = _nextrange(im.lastidx_extbuf, N)
    im.lastidx_extbuf = newlast
    return range
end
_nextrange(last, N) = last + N, last+1:last+N

function isdense(im::IndexManager)
    stateidxs = sizehint!(Int[], im.lastidx_dynamic)
    pidxs     = sizehint!(Int[], im.lastidx_p)
    outidxs   = sizehint!(Int[], im.lastidx_out)
    extidxs   = sizehint!(Int[], im.lastidx_extbuf)
    for dataranges in (im.v_data, im.e_data)
        for range in dataranges
            append!(stateidxs, range)
        end
    end
    for pararanges in (im.v_para, im.e_para)
        for range in pararanges
            append!(pidxs, range)
        end
    end
    for outranges in (im.v_out, im.e_out)
        for range in outranges
            append!(outidxs, flatrange(range))
        end
    end
    for extranges in (im.v_ext, im.e_ext)
        for range in extranges
            append!(extidxs, flatrange(range))
        end
    end
    sort!(pidxs)
    sort!(stateidxs)
    sort!(outidxs)
    sort!(outidxs)
    @assert pidxs == 1:im.lastidx_p
    @assert stateidxs == 1:im.lastidx_dynamic
    @assert outidxs == 1:im.lastidx_out
    @assert extidxs == 1:im.lastidx_extbuf
    return true
end

#=
SciMLBase gets the index provider from ODEFunction.sys which defaults to f.sys so we provide it...
SSI Maintainer assured that f.sys is really only used for symbolic indexig so method seems legit
=#
SciMLBase.__has_sys(nw::Network) = true
SciMLBase.__has_jac_prototype(nw::Network) = !isnothing(nw.jac_prototype)
function Base.getproperty(nw::Network, s::Symbol)
    if s===:sys
        nw
    elseif s===:jac_prototype
        getfield(nw, :jac_prototype)[]
    else
        getfield(nw, s)
    end
end

"""
    set_jac_prototype!(nw::Network, jac::SparseMatrixCSC{Bool,Int})

Set the Jacobian prototype for a NetworkDynamics network.

This function stores a pre-computed Jacobian sparsity pattern in the network object,
which can be used by ODE solvers to improve performance during integration.

# Arguments
- `nw::Network`: The NetworkDynamics network to modify
- `jac::SparseMatrixCSC{Bool,Int}`: A sparse matrix representing the Jacobian sparsity pattern
"""
function set_jac_prototype!(nw::Network, jac::SparseMatrixCSC{Bool,Int})
    getfield(nw,:jac_prototype)[] = jac
    nw
end

"""
    set_jac_prototype!(nw::Network; kwargs...)

Compute and set the Jacobian prototype for a NetworkDynamics network.

This is a convenience function that automatically computes the Jacobian sparsity pattern
using `get_jac_prototype` and stores it in the network object.
Needs `SparseConnectivityTracer` to be loaded!

# Arguments
- `nw::Network`: The NetworkDynamics network to modify
- `kwargs...`: Keyword arguments passed to `get_jac_prototype` (e.g., `dense`, `remove_conditions`)

# Example Usage
```julia
nw = Network(...)
set_jac_prototype!(nw) # computs sparsity pattern and stores in network
prob = ODEProblem(nw, x0, (0.0, 1.0), p0)
sol = solve(prob, Rodas5P())
```

See also: [`get_jac_prototype`](@ref)
"""
function set_jac_prototype!(nw::Network; kwargs...)
    jac = get_jac_prototype(nw; kwargs...)
    set_jac_prototype!(nw, jac)
end
