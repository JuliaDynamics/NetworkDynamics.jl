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

"""
   check(cond, msg)

If `cond` evaluates false throw `ArgumentError` and print evaluation of `cond`.
"""
macro check(cond::Expr, msg)
    head = lstrip(repr(cond), ':')
    head = head * " == false"
    args = ()
    for (i,a) in enumerate(cond.args[2:end])
        lhs = lstrip(repr(a), ':')
        symbol = (i == length(cond.args)-1) ? "└ " : "├ "
        args  = (args..., :("\n   " * $symbol * $lhs * " = " * repr($(esc(a)))))
    end
    return :($(esc(cond)) ||
             throw(ArgumentError($msg * "\n  " * $head * $(args...))))
end
