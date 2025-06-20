struct BatchStride{T}
    first::Int
    strides::T # either StaticInt, Tuple or NamedTuple of StaticInt
end
stridesT(::BatchStride{T}) where {T} = T
function BatchStride(first::Int, strides::Union{Tuple,NamedTuple})
    strides = map(Static.static, strides)
    BatchStride{typeof(strides)}(first, strides)
end
BatchStride(first::Int, stride::Int) = BatchStride(first, Static.static(stride))

# get full stride length (all substrides)
_fullstride(bs::BatchStride) = sum(bs.strides)

# full range for N elements with this stride
@inline function _fullrange(bs::BatchStride, N)
    (bs.first):(bs.first+N*_fullstride(bs)-1)
end

# range of i-th element
@inline function _range(bs::BatchStride{<:StaticInt}, i)
    start = bs.first + (i - 1) * _fullstride(bs)
    start:start+_fullstride(bs)-1
end
# subrange j of i-th element
@inline function _range(bs::BatchStride{<:Union{Tuple, NamedTuple}}, i, j)
    start = bs.first + (i - 1) * _fullstride(bs) + _sum_upto_dim(bs.strides, j)
    start:start+bs.strides[j]-1
end
@inline _sum_upto_dim(t::Tuple, idx, init=0) = sum(t[1:idx-1], init=0)
@inline _sum_upto_dim(t::NamedTuple, idx::Int) = _sum_upto_dim(values(t), idx)
@inline function _sum_upto_dim(t::NamedTuple, name::Symbol)
    idx = findfirst(isequal(name), keys(t))
    _sum_upto_dim(values(t), idx)
end

flatrange(r::AbstractRange) = r
function flatrange(rs::@NamedTuple{src::UnitRange{Int}, dst::UnitRange{Int}})
    @assert last(rs.src) + step(rs.src) == first(rs.dst)
    @assert step(rs.src) == step(rs.dst)
    UnitRange(first(rs.src), last(rs.dst))
end

subscript(N) = String(_subscript.(reverse(digits(N))))
_subscript(i) = Char(0x02080 + i)

@inline function unrolled_foreach(f::F1, filter::F2, t::Tuple) where {F1,F2}
    filter(first(t)) && f(first(t))
    @inline unrolled_foreach(f, filter, Base.tail(t))
end
@inline unrolled_foreach(f::F1, filter::F2, t::Tuple{}) where {F1,F2} = nothing
# Abstract Vector, no unrolling
@inline function unrolled_foreach(f::F1, filter::F2, t::AbstractVector) where {F1,F2}
    for el in t
        filter(el) && @noinline f(el) #noinline as function barrier for unknown batch type
    end
    nothing
end
# no filter
@inline unrolled_foreach(f, t) = unrolled_foreach(f, nofilt, t)

nofilt(_) = true

"""
    unique_mappings([f=identity], from, to)

Given two vectors `from` and `to`, find all keys in `from` which exist only once.
For those unique keys, return a dict maping `from_unique => f(to)`
"""
unique_mappings(from, to) = unique_mappings(identity, from, to)
function unique_mappings(f, from, to)
    counts = Dict{eltype(from),Int}()
    for k in from
        counts[k] = get(counts, k, 0) + 1
    end
    unique = Dict{eltype(from),eltype(to)}()
    for (k, v) in zip(from, to)
        if get(counts, k, 0) == 1
            unique[k] = f(v)
        end
    end
    unique
end


"""
    hash_fields(obj::T, h)

This is @generated helper function which unrolls all fields of a struct `obj` and
recursively hashes them.
"""
@generated function hash_fields(obj::T, h::UInt) where {T}
    fields = fieldnames(obj)
    subhashes = Expr(:block, (:(h = hash(obj.$field, h)) for field in fields)...)

    quote
        h = hash(T, h)
        $subhashes
        h
    end
end

"""
    equal_fields(a::T, b::T) where {T}

Thise @generated helper function unrolls all fields of two structs `a` and `b` and
compares them.
"""
@generated function equal_fields(a::T, b::T) where {T}
    fields = fieldnames(T)
    subequals = Expr(:block, (:(a.$field == b.$field || return false) for field in fields)...)

    quote
        $subequals
        return true
    end
end

function rand_inputs_fg(rng, cf)
    @argcheck hasindim(cf) "ComponentModel has no specified input dimensions/syms"
    du = [NaN for _ in 1:dim(cf)]
    u = rand(rng, dim(cf))
    p = rand(rng, pdim(cf))
    ins = Tuple(rand(rng, l) for l in values(indim(cf)))
    if has_external_input(cf)
        ext = rand(rng, extdim(cf))
        ins = (ins..., ext)
    end
    outs = Tuple([NaN for _ in 1:l] for l in values(outdim(cf)))
    t = NaN
    (outs, du, u, ins, p, t) # fg extrancs outs and ints internally!
end
rand_inputs_fg(cf) = rand_inputs_fg(Random.default_rng(), cf)

function rand_inputs_obsf(rng, cf)
    (_, u, p, ins, p, t) = rand_inputs_fg(rng, cf)
    N = length(cf.obssym)
    outs = [NaN for _ in 1:N]
    (outs, u, ins..., p, t)
end
rand_inputs_obsf(cf) = rand_inputs_obsf(Random.default_rng(), cf)


# abstract symbolic index types
abstract type SymbolicIndex{C,S} end
abstract type SymbolicStateIndex{C,S} <: SymbolicIndex{C,S} end
abstract type SymbolicParameterIndex{C,S} <: SymbolicIndex{C,S} end

flatten_sym(v::NamedTuple) = reduce(vcat, values(v))
flatten_sym(v::AbstractVector{Symbol}) = v

"""
    find_identical(v::Vector;; equality)

Find identical elements in a vector `v` using the `equality` function.
Returns a vector of vectors where each vector contains the indices of identical elements.
"""
function find_identical(v::Vector{T}; equality=isequal) where {T}
    indices = eachindex(v)
    idxs_per_type = Vector{Int}[]
    unique_comp = T[]
    for i in eachindex(v)
        found = false
        for j in eachindex(unique_comp)
            if equality(v[i], unique_comp[j])
                found = true
                push!(idxs_per_type[j], indices[i])
                break
            end
        end
        if !found
            push!(unique_comp, v[i])
            push!(idxs_per_type, [indices[i]])
        end
    end
    @assert length(unique_comp) == length(idxs_per_type)
    return idxs_per_type
end

####
#### SymbolicView helper type
####
"""
    SymbolicView{N,VT} <: AbstractVetor{VT}

Is a (smallish) fixed size vector type with named dimensions.
Its main purpose is to allow named acces to variables in
[`ComponentCondition`](@ref) and [`ComponentAffect`](@ref) functions.

I.e. when the `ComponentAffect` declared `sym=[:x, :y]`, you can
acces `u[:x]` and `u[:y]` inside the condition function.
"""
struct SymbolicView{N,VT} <: AbstractVector{VT}
    v::VT
    syms::NTuple{N,Symbol}
    allow_int_indexing::Bool
end

struct IllegalIntIndexingError <: Exception end

SymbolicView(v, syms) = SymbolicView(v, syms, true)
SymbolicView(v, syms::Vector, allow) = SymbolicView(v, Tuple(syms), allow)
Base.IteratorSize(::Type{SymbolicView}) = IteratorSize(x.v)
Base.IteratorEltype(::Type{SymbolicView}) = IteratorEltype(x.v)
Base.eltype(::Type{SymbolicView}) = eltype(x.v)
Base.size(x::SymbolicView) = size(x.v)
Base.firstindex(x::SymbolicView) = firstindex(x.v)
Base.lastindex(x::SymbolicView) = lastindex(x.v)
Base.iterate(x::SymbolicView) = iterate(x.v)
Base.iterate(x::SymbolicView, state) = iterate(x.v, state)
Base.length(x::SymbolicView{N}) where {N} = N
Base.IndexStyle(::Type{SymbolicView}) = IndexLinear()
Base.getindex(x::SymbolicView, index) = x.v[_sym_to_int(x, index)]
Base.setindex!(x::SymbolicView, value, index) = x.v[_sym_to_int(x, index)] = value
function _sym_to_int(x::SymbolicView, sym::Symbol)
    idx = findfirst(isequal(sym), x.syms)
    isnothing(idx) && throw(ArgumentError("SymbolError: try to access SymbolicView($(x.syms)) with symbol $sym"))
    idx
end
function _sym_to_int(x::SymbolicView, idx::Int)
    x.allow_int_indexing || throw(IllegalIntIndexingError())
    idx
end
_sym_to_int(x::SymbolicView, idx) = _sym_to_int.(Ref(x), idx)



####
#### Type definition for InitConstraint
####
"""
    struct InitConstraint{F}
    InitConstraint(f, sym, dim)

A representation of an additionalconstraint that is applied during the initialization phase of a component.
It contains a function `f` that defines the constraint, a vector of symbols `sym` that are involved in the constraint,
and the dimension `dim` of the constraint.

    InitConstraint([:x, :y], 2) do res, u
        res[1] = u[:x]^2 + u[:y]^2 - 1
    end

See also [`@initconstraint`](@ref) for a macro to create such constraints.
"""
struct InitConstraint{F}
    f::F
    sym::Vector{Symbol}
    dim::Int
    prettyprint::Union{Nothing,String}
end
InitConstraint(f, sym, dim) = InitConstraint(f, sym, dim, nothing)

"""
    InitConstraint(subconstraints::InitConstraint...)

Combine multiple `InitConstraint` objects into a single constraint.

The resulting constraint will have:
- Combined symbols from all subconstraints (deduplicated)
- Total dimension equal to sum of individual dimensions
- Function that evaluates all subconstraints sequentially

Limitation: All subconstraints need to use `Symbol` indices internally, i.e. u[:x] rather than u[2]!
"""
function InitConstraint(subconstraints::InitConstraint...)
    isempty(subconstraints) && throw(ArgumentError("At least one subconstraint must be provided"))
    if length(subconstraints) == 1
        return only(subconstraints)
    end

    all_syms = unique!(reduce(vcat, c.sym for c in subconstraints))
    total_dim = sum(dim(c) for c in subconstraints)
    f_range = []
    offset = 0
    for c in subconstraints
        tup = (f=c.f, range=(offset+1):(offset+dim(c)))
        offset += dim(c)
        push!(f_range, tup)
    end

    f_range_tup = Tuple(f_range)
    function combined_f(res, u)
        offset = 0
        unrolled_foreach(f_range_tup) do (cf, cr)
            res_view = view(res, cr) 
            cf(res_view, u)
        end
        nothing
    end

    try
        su = SymbolicView(zeros(length(all_syms)), all_syms, false)
        res = zeros(total_dim)
        combined_f(res, su)
    catch e
        if e isa IllegalIntIndexingError
            throw(ArgumentError(
                "Cannot combine multiple init constraints, because at least \
                 on of them uses `u[::Int]` indexing internally. Use `u[::Symbol]` \
                 within the `InitConstraint` to fix this error."
            ))
        else
            rethrow(e)
        end
    end
    
    InitConstraint(combined_f, all_syms, total_dim, nothing)
end

dim(c::InitConstraint) = c.dim

(c::InitConstraint)(res, u) = c.f(res, SymbolicView(u, c.sym))

function Base.show(io::IO, ::MIME"text/plain", @nospecialize(c::InitConstraint))
    if c.prettyprint === nothing
        r = repr(c)
        s = replace(r, ", nothing)"=>")")
        print(io, s)
    else
        print(io, c.prettyprint)
    end
end
