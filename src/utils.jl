struct BatchStride
    first::Int
    stride::Int
end
@inline function _range(bs::BatchStride, i)
    start = bs.first + (i - 1) * bs.stride
    start:start+bs.stride-1
end
@inline function _fullrange(bs::BatchStride, N)
    (bs.first):(bs.first+N*bs.stride-1)
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
