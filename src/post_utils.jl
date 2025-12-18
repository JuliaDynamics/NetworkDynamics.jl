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

Custom constructor for creating ODEProblem base of a `Network`-Object.
Its main purpose is to automatically handle callback construction from the component level callbacks.

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

####
#### Loopback Connections
####
function LOOPBACK_G(outdst, insrc, indst, p, t)
    outdst .= -1 .* insrc
    nothing
end

"""
    LoopbackConnection(; potential, flow, kwargs...)

A `LoopbackConnection` is a special `EdgeModel` that enables direct connection of
"injector nodes" to a "hub" node without requiring aggregation logic.
An injector node is an "inverted" `VertexModel`, which gets the networks potential as an
input and outputs a flow variable.

The LoopbackConnection allows a direct, star-like connection of injector nodes
to a single hub nodes.
The LoopbackConnection is a **directed edge model from injector to hub**!

```asciiart
       ┊
    ┄┄┄◯   ● injector 1
      ╱ ╲ ╱
   ┄┄◯╶─╴◯╶─╴● injector 2
     ┊    ╲
           ● injector 3
```

Injector nodes:
- have a flipped interface (potential in, flow out)
- must be leaf nodes (one neighbor only),
- must be connected through a LoopbackConnection EdgeModel and
- may have feed-forward (direct dependency of flow-output on potential-input).

!!! note "Sign Convention"
    For normal vertices, positive flow as an input means flow *into* the vertex.
    This convention is maintained for injector nodes (though it may seem counter-intuitive):
    - **Positive flow**: Draw from the hub (consumption)
    - **Negative flow**: Injection into the hub (production)

    When using ModelingToolkit models, you only need to flip the input/output variable
    declarations—the equations themselves remain unchanged. For example, a resistor with
    `p.i ~ p.v/R` keeps the same equation; only the interface changes from
    `VertexModel(..., [:p₊i], [:p₊v])` to `VertexModel(..., [:p₊v], [:p₊i])`.

```asciiart
                       △
      ╭────────────────┼────────────╮
      │      potential │ φ out      │
━━━━━━▽━━┓   ╔═════════△═════════╗  │  ┏━━━━━━━┓   ╔═══════════════════╗
 normal  ┃   ║ VertexModel (hub) ║  ╰──▷┄┄┄┄┄┄┄▷───▷ Injector Vertex   ║
EdgeModel┃   ║ ẋ = f(x, Φ, p, t) ║     ┃       ┃   ║ ẋ = f(x, φ, p, t) ║
         ┃   ║ φ = g(x, p, t)    ║  ╭──◁┄×(-1)┄◁───◁ Φ = g(x, φ, p, t) ║
━━━━━━▽━━┛   ╚═════════△═════════╝  │  ┗━━━━━━━┛   ╚═══════════════════╝
 flow │ Φ out        ╭─┴─╮          │  special      ⋅ flipped interface:
      ╰──────────────▷ + ◁──────────╯  "Loopback"     ▷ potential φ in
       (aggregation) ╰─△─╯             EdgeModel      ◁ flow Φ out
                       │               inj => hub   ⋅ feed forward allowed
```
For input-output naming you need to provide the `potential` and `flow` symbols.

```julia-repl
julia> LoopbackConnection(; potential=[:u_r, :u_i], flow=[:i_r, :i_i], src=1, dst=2)
EdgeModel :loopback PureFeedForward() @ Edge 1=>2
 ├─ 2/2 inputs:  src=[injector₊i_r, injector₊i_i] dst=[hub₊u_r, hub₊u_i]
 ├─   0 states:  []
 └─ 2/2 outputs: src=[injector₊u_r, injector₊u_i] dst=[hub₊i_r, hub₊i_i]
```

"""
function LoopbackConnection(; potential, flow, kwargs...)
    insym = (;
        src=[Symbol(:injector₊, s) for s in flow],
        dst=[Symbol(:hub₊, s) for s in potential],
    )
    outsym = (;
        src=[Symbol(:injector₊, s) for s in potential],
        dst=[Symbol(:hub₊, s) for s in flow],
    )

    g = Directed(NetworkDynamics.LOOPBACK_G)
    EdgeModel(; g, insym, outsym, check=false, name=:loopback, kwargs...)
end

is_loopback(eb::ComponentBatch) = isnothing(compf(eb)) && compg(eb) == NetworkDynamics.LOOPBACK_G
is_loopback(em::EdgeModel) = em.g isa Directed && em.g.g == NetworkDynamics.LOOPBACK_G
has_loopback_edges(im::IndexManager) = any(is_loopback, im.edgem)

function gen_loopback_map(im::IndexManager)
    outindex = Int[]
    aggindex = Int[]
    for i in 1:ne(im.g)
        is_loopback(im.edgem[i]) || continue
        e = im.edgevec[i]
        # we want to copy from dst output to src input
        append!(outindex, im.v_out[e.dst])  # output idx of dst vertex
        append!(aggindex, im.v_aggr[e.src]) # input idx of src vertex
    end
    outindex, aggindex
end

apply_loopback!(aggbuf, obuf, map) = _apply_loopback!(get_backend(aggbuf), aggbuf, obuf, map)

function _apply_loopback!(::KernelAbstractions.CPU, aggbuf, obuf, map)
    outindex, aggindex = map
    for (oi, ai) in zip(outindex, aggindex)
        aggbuf[ai] = obuf[oi]
    end
    nothing
end

function _apply_loopback!(backend::KernelAbstractions.GPU, aggbuf, obuf, map)
    outindex, aggindex = map
    kernel = _lb_kernel!(backend)
    kernel(aggbuf, obuf, outindex, aggindex; ndrange=length(outindex))
    nothing
end
@kernel function _lb_kernel!(aggbuf, @Const(obuf), @Const(outindex), @Const(aggindex))
    I = @index(Global)
    @inbounds if I ≤ length(outindex)
        aggbuf[aggindex[I]] = obuf[outindex[I]]
    end
    nothing
end
