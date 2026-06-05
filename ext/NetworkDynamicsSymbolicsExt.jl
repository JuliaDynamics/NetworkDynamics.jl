module NetworkDynamicsSymbolicsExt

using Symbolics: Symbolics, @variables
using MacroTools: postwalk, @capture
using NetworkDynamics: NetworkDynamics, SymbolicIndex, SII, ParameterScope, ParameterDefaultFrom

# Register `scope` as symbolic variable metadata, so that it can be used in
# `@variables`/`@parameters`/`@mtkmodel`/`@component` style definitions, e.g.
# `@parameters p [scope=:global]`.
Symbolics.option_to_metadata_type(::Val{:scope}) = ParameterScope

# Register `default_from` as symbolic variable metadata, e.g.
# `@parameters p [default_from=:busbar]` or `@parameters p [default_from=(:src, :busbar)]`.
Symbolics.option_to_metadata_type(::Val{:default_from}) = ParameterDefaultFrom

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
        if x ∈ (:+, :-, :^, :/, :.+, :.-, :.^, :./, Symbol(":")) || x isa QuoteNode
            x
        else
            :(NetworkDynamics.collect_symbol!($mapping_sym, $x))
        end
    end
    quote
        $(esc(mapping_sym)) = Dict()
        expr = $(esc(symbolic_expr))
        NetworkDynamics.ObservableExpression($(esc(mapping_sym)), expr, $(Meta.quot(name)))
    end
end

function NetworkDynamics.ObservableExpression(mapping, expr::Vector, names)
    if names isa Symbol
        names = [Symbol(names, NetworkDynamics.subscript(i)) for i in 1:length(expr)]
    elseif isnothing(names)
        names = Iterators.repeated(nothing, length(expr))
    end
    [NetworkDynamics.ObservableExpression(mapping, e, n) for (e, n) in zip(expr, names)]
end

function NetworkDynamics.ObservableExpression(mapping, expr, name)
    nwidx = collect(keys(mapping))
    syms = collect(values(mapping))
    f = Symbolics.build_function(expr, syms; expression=Val(false))
    NetworkDynamics.ObservableExpression(nwidx, f, expr, name)
end

NetworkDynamics.collect_symbol!(_, x) = x
function NetworkDynamics.collect_symbol!(mapping, x::Union{AbstractVector, Tuple})
    map(el -> NetworkDynamics.collect_symbol!(mapping, el), x)
end
function NetworkDynamics.collect_symbol!(mapping, nwindex::SymbolicIndex)
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
