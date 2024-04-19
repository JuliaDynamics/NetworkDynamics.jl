function edgebyidx(g::AbstractGraph, i)
    iter = edges(g)
    i ≤ length(iter) || throw(BoundsError(collect(edges(g)), i))
    el, state = iterate(iter)
    for _ ∈ 1:(i-1)
        el, state = iterate(iter, state)
    end
    el
end
