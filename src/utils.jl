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
