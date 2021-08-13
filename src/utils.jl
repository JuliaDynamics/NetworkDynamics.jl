struct CachePool
    caches::Dict{Any, AbstractArray}
    CachePool() = new(Dict{Any, AbstractArray}())
end

function getcache(c::CachePool, T, length)::T
    key = (T, length)
    return get!(c.caches, key) do
        T(undef, length)
    end::T
end

"""
   check(cond, msg)

If `cond` evaluates false throw `ArgumentError` and print evaluation of `cond`.
"""
macro check(cond::Expr, msg)
    head = lstrip(repr(cond), ':')
    head = head * " evaluated false"
    args = ()
    for (i,a) in enumerate(cond.args[2:end])
        lhs = lstrip(repr(a), ':')
        symbol = (i == length(cond.args)-1) ? "└ " : "├ "
        args  = (args..., :("\n   " * $symbol * $lhs * " = " * repr($(esc(a)))))
    end
    return :($(esc(cond)) ||
             throw(ArgumentError($(esc(msg)) * "\n  " * $head * $(args...))))
end

"""
    @cond_threads

Mix @inbounds and @Threads.threads
Allows control over threading by the 1st argument, a boolean (likely from the
network object). This allows you to control for multithreading at runtime
without code duplication.
"""
macro cond_threads(cond,args...)
    na = length(args)
    if na != 1
        throw(ArgumentError("wrong number of arguments in @cond_threads"))
    end
    ex = args[1]
    if ex.head == :for
        # Need to escape the bounds in the loop, but not the whole expression,
        # in order to work with the base Threads.@thread macro.
        # That macro does it's own escaping, but also it's own expression
        # checking (the head of the first expression must be :for. So we will
        # escape the bounds here, keeping the expression formed in a suitable
        # way for @threads

        # make sure we have the form :for, iterator, block
        @assert ex.args[1].head === :(=)
        @assert ex.args[2].head === :block
        @assert length(ex.args) == 2

        # extract the loop iterator
        loop_iter = Expr(:(=), esc.(ex.args[1].args)...)
        # extract the loop body and annotate @inbounds
        loop_body_ib = :(@inbounds $(esc(ex.args[2])))
        # build new loop
        loop = Expr(:for, loop_iter, loop_body_ib)

        return quote
            if $(esc(cond))
                Threads.@threads $loop
            else
                @inbounds $(esc(ex))
            end
        end
    else
        throw(ArgumentError("unrecognized argument to @cond_threads"))
    end
end
