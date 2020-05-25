# Parameter handling

Let `nd!` be an ODEFunction returned by `network_dynamics`, e.g.

```julia
nd! = network_dynamics(vertices!, edges!, graph)
```

 then the behaviour of `nd!` has the signature `(dx, x, t, p)` and its behaviour changes with the type of parameters `p` being passed.

  * When `p` is an Array, a Dict, a struct... then the entire object is passed to each `VertexFunction` and `EdgeFunction`.
  * When `p = (p_v, p_e)` is a Tuple of two values, then the first value will be passed to all vertices and the second to all edges.
  * If `p = (p_v_arr, p_e_arr)` is a Tuple of two Arrays with lengths corresponding to the number of nodes and number of edges respectively, then the edges or nodes receive only the parameter with the corresponding index.
  * If all nodes and/or edges have no internal parameters the value `nothing` may be passed. Using `nothing` instead of dummy parameters is usually faster, since then less data are copied.

Another option for specifying heterogeneous parameters is to make each `VertexFunction` a callable struct with the parameters hardcoded as fields. This approach is used in [PowerDynamics.jl](https://github.com/JuliaEnergy/PowerDynamics.jl). However it provides considerably less flexibility and interoperability with other packages.

For its greater speed and flexibility in modeling we recommend to use the tuple syntax.



## Compatability with DiffEqFlux

Most of the sensitivity algorithms that DiffEqFlux makes use of assume that the parameters `p` are a subtype of `AbstractArray`. Therefore they are not compatible with the tuple syntax.`TrackerAdjoint` does not have these limitations but works best on out-of-place problems. Unfortunately, `network_dynamics` returns an ODEProblem that is in-place (mutating its inputs) leading to slow performance with Tracker.

However, wrapping the function in such a way that it accepts arrays of parameters that are later pasted into the tuple syntax sidesteps these issues and enables the use of adjoint methods. Depending on the use case such a wrapper might look like this:

```julia
function nd_wrapper!(dx, x, p, t)
  nd!(dx, x, (p, nothing), t)
end
```

`nd_wrapper!` should now work with `BacksolveAdjoint` and `InterpolatingAdjoint`. At the moment we recommend this way for combining NetworDynamics and DiffEqFlux.

Forward mode (ForwardDiff.jl) and source-to-source (Zygote.jl)  automatic differentiation  is not fully-supported yet. For more detailed discussion see [this issue](https://github.com/FHell/NetworkDynamics.jl/issues/34).

Further resources:

* [DiffEqSensitivity algorithms](https://docs.sciml.ai/stable/analysis/sensitivity/#Sensitivity-Algorithms-1)
