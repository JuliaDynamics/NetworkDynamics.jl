struct CachePool
    caches::Dict{Any, AbstractArray}
    CachePool() = new(Dict{Any, AbstractArray}())
end

function getcache(c::CachePool, T, size...)::T
    key = (T, size)
    return get!(c.caches, key) do
	   T(undef, size...)
    end::T
end
