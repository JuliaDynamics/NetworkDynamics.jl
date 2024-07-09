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

@inline function unrolled_foreach(f::F, t::Tuple) where {F}
    f(first(t))
    @inline unrolled_foreach(f, Base.tail(t))
end
@inline unrolled_foreach(f::F, t::Tuple{}) where {F} = nothing
@inline unrolled_foreach(f::F, t::AbstractVector) where {F} = @inline foreach(f, t)
