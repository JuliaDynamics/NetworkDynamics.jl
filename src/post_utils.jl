"""
    SciMLBase.ODEProblem(nw::Network, s0::NWState, tspan;
        additional_cb=Dict(),
        override_cb=nothing,
        kwargs...)
    SciMLBase.ODEProblem(nw::Network, s0::NWState, tspan, p0::NWParameter;
        additional_cb=Dict(),
        kwargs...)

This is a simple wrapper which:
- extracts the flat state and parameter vectors from `s0` (and `p0` if provided)
- makes a copy of the parameter vector to avoid side effects due to callbacks
- constructs the callbacks from the network and any additional callbacks provided.
  `additional_cb` is forwarded to [`get_callbacks`](@ref) and can be used to inject
  component callbacks which are not part of the metadata.

If you need to pass "normal" callbacks, please use the underlying `ODEProblem` constructor directly.
To avoid confusion, passing `callback` will throw an error.

This function is a convenience wrapper around it:

    ODEProblem(nw, uflat(s0), tspan, copy(pflat(p0));
               callback=get_callbacks(nw, additional_cb), kwargs...)

where `p0 = s0.p` if not provided.
"""
function SciMLBase.ODEProblem(
    nw::Network, args...;
    add_comp_cb=Dict(),
    add_nw_cb=nothing,
    override_cb=nothing,
    kwargs...
)

    if haskey(kwargs, :callback)
        throw(ArgumentError("Cannot pass `callback` keyword to ODEProblem(nw::Network, ...) constructor. \
                             Callbacks are allways generated from the Network object using `get_callbacks(nw)`. \
                             You can either\n
                                - pass additional component callbacks using `add_comp_cb=Dict(VIndex(1)=>comp_callback)`\n
                                - pass additional network level callbacks using `add_nw_cb=callback/CallbackSet` or\n
                                - override all callbacks using `override_cb=callback/CallbackSet`"))
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

    SciMLBase.ODEProblem(nw, args...; callback=finalcallback, kwargs...)
end


function SciMLBase.ODEProblem(nw::Network, s0::NWState, tspan, p::NWParameter=s0.p; kwargs...)
    u = uflat(s0)
    p = copy(pflat(p))

    SciMLBase.ODEProblem(nw, u, tspan, p; kwargs...)
end
