module SymbolicsExt

using Symbolics: Symbolics, @variables, Num
using MacroTools: postwalk, @capture
using NetworkDynamics: NetworkDynamics, ObservableExpression, SymbolicVertexIndex, SymbolicEdgeIndex, SII

function NetworkDynamics.generate_observable_expression(ex::Expr)
    if @capture(ex, capname_ = capcontent_)
        name = capname
        content = capcontent
    else
        name = nothing
        content = ex
    end
    # manually hygen on "mapping", otherwise we cannot escape x
    mapping_sym = gensym("mapping")
    symbolic_expr = postwalk(content) do x
        :(NetworkDynamics.collect_symbol!($mapping_sym, $x))
    end
    quote
        $(esc(mapping_sym)) = Dict()
        expr = $(esc(symbolic_expr))
        ObservableExpression($(esc(mapping_sym)), expr, $(Meta.quot(name)))
    end
end

function NetworkDynamics.ObservableExpression(mapping, expr, name)
    nwidx = collect(keys(mapping))
    syms = collect(values(mapping))
    f = Symbolics.build_function(expr, syms; expression=Val(false))
    ObservableExpression(nwidx, f, expr, name)
end

NetworkDynamics.collect_symbol!(_, x) = x
function NetworkDynamics.collect_symbol!(mapping, x::AbstractVector)
    if length(x) == 1 && only(x) isa Union{SymbolicVertexIndex,SymbolicEdgeIndex}
        return NetworkDynamics.collect_symbol!(mapping, only(x))
    elseif any(el -> el isa Union{SymbolicVertexIndex,SymbolicEdgeIndex}, x)
        throw(ArgumentError("Cannot handle vector expressions in @obsex! Encountered $x."))
    else
        x
    end
end
function NetworkDynamics.collect_symbol!(mapping, nwindex::Union{SymbolicVertexIndex,SymbolicEdgeIndex})
    if haskey(mapping, nwindex)
        return mapping[nwindex]
    else
        name = SII.getname(nwindex)
        sym = only(@variables $name)
        mapping[nwindex] = sym
        return sym
    end
end

end # module
