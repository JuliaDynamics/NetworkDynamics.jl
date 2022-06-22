# [Parameter handling](@id parameters)

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



## Compatability with specific packages

Some other packages from the Julia ecosystem, e.g. DiffEqFlux, assume that the parameters `p` are a subtype of `AbstractArray`. Therefore they are not fully compatible with the tuple syntax.

However, wrapping the function in such a way that it accepts arrays of parameters that are later pasted into the tuple syntax sidesteps potential issues. Depending on the use case such a wrapper might look like this:

```julia
function nd_wrapper!(dx, x, p, t)
  nd!(dx, x, (p, nothing), t)
end
```
At the moment we recommend this way for combining NetworkDynamics and DiffEqFlux.