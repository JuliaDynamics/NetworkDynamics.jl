# Accessing edge values

 In the following, we will present two different ways how the edge values of an `ODEFunction` with `StaticEdge`s  can be retrieved

## Accessing edge values via GetGD

`NetworkDynamics.jl` includes the types `GetGD` and `GetGS`, which can be used to access to the underlying GraphData and GraphStruct objects by multiple dispatching of the generated `ODEFunction`. Let `nd` be the ODEFunction that has been formed by calling `network_dynamics`. To get access to the edge values at any time $t$, `nd` can be called as follows:

```@example accessing_edge_values
gd_nd = nd(x, p, t, GetGD) # exposes underlying graph data struct
e_values = gd_nd.gdb.e_array
nothing # hide
```

By adding `GetGD` as an argument to the call of `nd`, a `GraphData` instance at the given time, parameters and states will be returned. `gdb` is the internal `GraphDataBuffer` struct. The argument `x` has to be the precomputed state of the system at time `t`, e.g. from a solution object.

## Accessing edge values using SavingCallback

Instead of recomputing the values of the edge variables they can be saved in parallel to an integration with `DiffEqCallbacks`.

```@example
using DiffEqCallbacks: SavingCallback, SavedValues

saved_values = SavedValues(Float64, Vector{Float64})
function saving_func(u, t, integrator)
    edgevals = Float64[]
    for i in 1:integrator.f.f.graph_structure.num_e
        push!(edgevals, integrator.f.f.graph_data.gdb.e[i]...)
    end
    edgevals
end
cb = SavingCallback(saving_func, saved_values)
sol = solve(prob, Tsit5(), callback=cb)
```

The variables `saved_values` will contain the stored edge values. At the moment this method requires detailed knowledge of the internal data structures. Our plan is to simplify this in a future release.
