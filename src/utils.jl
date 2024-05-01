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
