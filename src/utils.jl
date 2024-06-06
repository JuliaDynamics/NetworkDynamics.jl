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

# function idxs_containing(nd, ex)
#     allidx = SII.variable_symbols(nd)
# end
function _extract_nw(inpr)
    sc = SII.symbolic_container(inpr)
    if sc isa SciMLBase.ODEFunction
        sc.sys
    else
        sc
    end
end
function vertex_idxs(inpr; static=true, filter=nothing)
    nw = _extract_nw(inpr)
    syms = []
    for (i,cf) in enumerate(nw.im.vertexf)
        static || NetworkDynamics.isdynamic(cf) || continue
        append!(syms, collect(VIndex(i, sym(cf))))
    end
    if isnothing(filter)
        syms
    else
        Base.filter(filter, syms)
    end
end
function edge_idxs(inpr; static=true, filter=nothing)
    nw = _extract_nw(inpr)
    syms = []
    for (i,cf) in enumerate(nw.im.edgef)
        static || NetworkDynamics.isdynamic(cf) || continue
        append!(syms, collect(EIndex(i, sym(cf))))
    end
    if isnothing(filter)
        syms
    else
        Base.filter(filter, syms)
    end
end
Base.contains(s::SymbolicIndex, ex) = contains(string(s.subidx), ex)
Base.contains(s::SymbolicIndex, ex::Symbol) = contains(string(s.subidx), string(ex))
