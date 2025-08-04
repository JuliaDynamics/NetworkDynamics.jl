# Symbolic Indexing

By using SciML's [`SymbolicIndexingInterface.jl`](https://github.com/SciML/SymbolicIndexingInterface.jl), `ND.jl` 
provides numerous methods to access and change variables and parameters.

## Provide Symbol Names
When constructing component models, you can pass symbolic names using the `sym` and `psym` keywords.
```@example si
using NetworkDynamics, Graphs, OrdinaryDiffEqTsit5, Plots
function _edgef!(e, v_s, v_d, (K,), t)
    e .= K * (v_s[1] .- v_d[1])
end
edgef = EdgeModel(;g=AntiSymmetric(_edgef!), outsym=[:flow], psym=[:K=>1])
```
Here we create a static diffusion edge with suitable variable and parameter names.
Similarly, we define the diffusion vertex with symbolic names.
```@example si
function _vertexf!(dv, v, esum, p, t)
    dv[1] = esum[1]
end
vertexf = VertexModel(f=_vertexf!, g=1, sym=[:storage])
```


## Fundamental Symbolic Indices
The default types for this access are the types [`VIndex`](@ref), [`EIndex`](@ref), [`VPIndex`](@ref) and [`EPIndex`](@ref).
Each of those symbolic indices consists of 2 elements: a reference to the network component and a reference to the 
symbol within that component.
As such: 
- `VIndex(2, :x)` refers to variable with symbolic name `:x` in vertex number 2.
- `EPIndex(4, 2)` refers to the *second* parameter of the edge component for the 4th edge.

!!! details "Setup code to make following examples work"
    ```@example si
    g = wheel_graph(5)
    nw = Network(g, vertexf, edgef)
    s = NWState(nw)
    s.v[:,:storage] .= randn(5)
    prob = ODEProblem(nw, uflat(s), (0,2), pflat(s))
    sol = solve(prob, Tsit5()) 
    nothing #hide
    ```

Those fundamental indices can be used in a lot of scenarios. Most importantly you can use them to
```@example si
sol(sol.t; idxs=VIndex(1, :storage))   # extract timeseries out of a solution object
plot(sol; idxs=[VIndex(1, :storage), VIndex(5,:storage)]) # plot storage of two nodes
```

## Generate Symbolic Indices
Often, you need many individual symbolic indices. To achieve this you can use the helper methods [`vidxs`](@ref), 
[`eidxs`](@ref), [`vpidxs`](@ref) and [`epidxs`](@ref). With their help you can generate arrays of symbolic indices:

```@example si
vidxs(nw, :, :storage) # get variable "storage" for all vertices
```
```@example si
plot(sol; idxs=vidxs(nw, :, :storage))
```

## `NWState` and `NWParameter` Objects
Internally, both state and parameters of a `Network` are represented using flat arrays.
To access the state or parameters of a network, you can use the [`NWState`](@ref) and [`NWParameter`](@ref) objects.
```@example si
p = NWParameter(nw)
```
creates a `NWParameter` object for the network `nw`.
It essentially creates a new flat parameter array and fills it with the default parameter values defined in the component.
The parameters in the `NWParameter` object can be accessed using symbolic indices.
```@example si
p[EPIndex(5, :K)] = 2.0 # change the parameter K of the 5th edge
nothing #hide
```
Similarly, you can create a `NWState` object for the network `nw` using
```@example si
s = NWState(nw)
```
No default values were provided in the network components, so the state array is filled with `NaN` values.
```@example si
s[VIndex(:, :storage)] .= randn(5) # set the (initial) storage for all vertices 
s #hide
```
For both `NWState` and `NWParameter` objects, there is a more convenient way to access the variables and parameters.
```@example si
@assert s.v[1, :storage] == s[VIndex(1, :storage)] # s.v -> access vertex states
@assert s.e[1, :flow]    == s[EIndex(1, :flow)]    # s.e -> access edge states
@assert s.p.e[1,:K]      == p[EPIndex(1, :K)]      # s.p -> access parameters
```

The `NWState` and `NWParameter` objects are mutable, thus changing them will also change the underlying wrapped flat arrays.
You can always access the flat representations by calling [`uflat`](@ref) and [`pflat`](@ref). The ordering of elements 
in these flat arrays corresponds exactly to the order returned by [`variable_symbols`](@ref) and 
[`parameter_symbols`](@ref) respectively.

!!! note
    The `NWState` and `NWParameter` wrappers can be constructed from various objects.
    For example, within a callback you might construct `p = NWParameter(integrator)` to then change the parameters of 
the network within the callback.


## Observables
Sometimes, the "states" you're interested in aren't really states in the DAE sense but rather
algebraic derivations from DAE states, parameters, and time -- in accordance with the naming in 
the `SciML` ecosystem, these values are called Observables.

A prime example of Observables are edge/vertex-outputs, such as the `flow` in the edge model defined above.
It is also possible to define additional Observables manually by using the `obssym` and `obsf` keyword
on the `EdgeModel`/`VertexModel` constructors.
When building models using ModelingToolkit, the reduced algebraic states will be preserved automatically as observables.

Observables can be accessed like any other state. For example, the flows in the network don't show up in the state array 
but can be accessed in all the ways discussed above. 
For example:
```@example si
plot(sol; idxs=eidxs(nw, :, :flow))
```

## Derived `ObservableExpressions` using `@obsex`

Sometimes it is useful to plot or observe simple derived quantities. For that,
one can use the [`@obsex`](@ref) macro to define simple derived quantities.

For example, we can directly plot the storage difference with respect to storage of node 1.

```@example si
plot(sol; idxs=@obsex(vidxs(nw,:,:storage) .- VIndex(1,:storage)))
```

Other examples include calculating the magnitude and argument of complex values that are modeled using real and 
imaginary parts.
```
@obsex mag = sqrt(VIndex(1, :u_r)^2 + VIndex(2, :u_i)^2)
```

## Low-level accessors for flat array indices
Sometimes, you want to know the indices of your states in the flat arrays.
For that, you can use the low-level methods defined in `SymbolicIndexingInterface.jl`:

```@example si
using NetworkDynamics: SII # SII = SymbolicIndexingInterface
idxs = SII.variable_index(nw, vidxs(1:2, :storage))
```
```@example si
uflat(s)[idxs] == s.v[1:2, :storage]
```
Analogous with parmeters:
```@example si
idxs = SII.parameter_index(nw, eidxs(1:2, :K))
pflat(s)[idxs] == s.p.e[1:2, :K]
```

If you need the symbols of all the states/parameters in order, you can use:
```@example si
SII.variable_symbols(nw)
```
and
```@example si
SII.parameter_symbols(nw)
```
These functions return the symbolic indices in the exact order they appear in the flat arrays
returned by [`uflat`](@ref) and [`pflat`](@ref), making them essential when you need to map
between flat array indices and symbolic representations.

All above examples also work on other "symbolic containers", e.g. `SII.variable_symbols(::NWState)`.
