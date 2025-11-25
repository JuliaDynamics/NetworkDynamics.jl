const callback_keyword_docs = """
## Callback Keywords

The callback system supports three keyword arguments that control how callbacks are managed:

- **`add_comp_cb`**: Additional component callbacks: A `Dict` mapping component
  indices (e.g., `VIndex(1)` or `EIndex(2)`) to component callbacks. These are
  forwarded to [`get_callbacks`](@ref) and combined with callbacks stored in the
  network's component metadata. Use this to inject temporary component callbacks
  without modifying the network structure.

- **`add_nw_cb`**: Additional network/system callbacks: A network-level callback
  or `CallbackSet` (e.g., `PeriodicCallback`, `PresetTimeCallback`) that is
  combined with the network's component callbacks. Use this for callbacks that
  don't fit the component-based pattern, such as periodic saving or global
  termination conditions.

- **`override_cb`**: A callback or `CallbackSet` that completely replaces all network callbacks.
  When set, both `add_comp_cb` and `add_nw_cb` must be empty/nothing (enforced by ArgumentError).
  Use this for complete control over the callback system.
"""

"""
    SciMLBase.ODEProblem(nw::Network, args...;
        add_comp_cb=Dict(),
        add_nw_cb=nothing,
        override_cb=nothing,
        kwargs...
    )

Custom cosntructor for creating ODEProblem base of a `Network`-Object.
Its main purpose is to automaticially handle callback construction from the component level callbacks.

$callback_keyword_docs
"""
function SciMLBase.ODEProblem(
    nw::Network, args...;
    add_comp_cb=Dict(),
    add_nw_cb=nothing,
    override_cb=nothing,
    kwargs...
)

    if haskey(kwargs, :callback)
        if typeof(kwargs[:callback]) == typeof(get_callbacks(nw))
            @warn "Passing `callback=get_callbacks(nw)` to ODEProblem(nw, ...) is deprecated. The ODEConstructor will allways extract the Network callbacks automaticially."
            kwargs = filter(kv -> kv.first != :callback, kwargs)
            @assert !haskey(kwargs, :callback)
        else
            throw(ArgumentError("""
            Cannot pass `callback` keyword to ODEProblem(nw::Network, ...) constructor. Callbacks are always generated from the Network object using `get_callbacks(nw)`. You can either
            - pass additional component callbacks using `add_comp_cb=Dict(VIndex(1)=>comp_callback)`
            - pass additional network level callbacks using `add_nw_cb=callback/CallbackSet` or
            - override all callbacks using `override_cb=callback/CallbackSet`
            """))
        end
    end
    if !isnothing(override_cb) && (!isempty(add_comp_cb) || !isnothing(add_nw_cb))
        throw(ArgumentError("Cannot pass `override_cb` together with `add_comp_cb` or `add_nw_cb`. When overriding the default network callbacks, no additional callbacks are allowed."))
    end

    if !isnothing(override_cb)
        finalcallback = override_cb
    else
        nw_callback = get_callbacks(nw, add_comp_cb)
        if isnothing(add_nw_cb)
            finalcallback = nw_callback
        else
            finalcallback = SciMLBase.CallbackSet(nw_callback, add_nw_cb)
        end
    end

    SciMLBase.ODEProblem(SciMLBase.ODEFunction(nw), args...; callback=finalcallback, kwargs...)
end

"""
    SciMLBase.ODEProblem(nw::Network, s0::NWState, tspan; kwargs...)
    SciMLBase.ODEProblem(nw::Network, s0::NWState, tspan, p0::NWParameter; kwargs...)

This is a simple wrapper which:
- extracts the flat state and parameter vectors from `s0` (and `p0` if provided)
- makes a copy of the parameter vector to avoid side effects due to callbacks
- constructs the callbacks from the network and combines them with any additional callbacks

$callback_keyword_docs
"""
function SciMLBase.ODEProblem(nw::Network, s0::NWState, tspan, p::NWParameter=s0.p; kwargs...)
    u = uflat(s0)
    p = copy(pflat(p))

    SciMLBase.ODEProblem(nw, u, tspan, p; kwargs...)
end
