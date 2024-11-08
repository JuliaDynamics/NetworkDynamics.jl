struct BatchStride{N}
    first::Int
    strides::NTuple{N,Int}
end
BatchStride(first::Int, stride::Int) = BatchStride(first, (stride,))

# get full stride length (all substrides)
_fullstride(bs::BatchStride{1}) = @inbounds bs.strides[1]
_fullstride(bs::BatchStride) = sum(bs.strides)

# full range for N elements with this stride
@inline function _fullrange(bs::BatchStride, N)
    (bs.first):(bs.first+N*_fullstride(bs)-1)
end

# range of i-th element
@inline function _range(bs::BatchStride{1}, i)
    start = bs.first + (i - 1) * _fullstride(bs)
    start:start+_fullstride(bs)-1
end
# subrange j of i-th element
@inline function _range(bs::BatchStride{N}, i, j) where {N}
    start = bs.first + (i - 1) * _fullstride(bs) + sum(bs.strides[1:j-1], init=0)
    start:start+bs.strides[j]-1
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
    @inline foreach(f, Iterators.filter(filter, t))
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
    @argcheck hasindim(cf) "ComponentFunction has no specified input dimensions/syms"
    du = rand(rng, dim(cf))
    u = rand(rng, dim(cf))
    p = rand(rng, pdim(cf))
    ins = Tuple(rand(rng, l) for l in values(indim(cf)))
    outs = Tuple(rand(rng, l) for l in values(outdim(cf)))
    t = NaN
    (outs..., du, u, ins..., p, t)
end
rand_inputs_fg(cf) = rand_inputs_fg(Random.default_rng(), cf)
