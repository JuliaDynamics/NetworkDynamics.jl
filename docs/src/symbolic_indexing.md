# Symbolic Indexing

Using SciML's [`SymblicIndexingInterface.jl`](https://github.com/SciML/SymbolicIndexingInterface.jl), `ND.jl` provides lots of methods to access and change variables and Parameters.

!!! details "Setup code to make following examples work"
    ```@example si
    using NetworkDynamics
    using Graphs
    using OrdinaryDiffEqTsit5
    using Plots
    ```

## Provide Symbol Names
When construction component functions, you can pass symbolic names using the `sym` and `psym` keywords.
```@example si
function _edgef!(e, v_s, v_d, (K,), t)
    e .= K * (v_s[1] .- v_d[1])
end
edgef = StaticEdge(_edgef!; sym=[:flow], psym=[:K=>1], coupling=AntiSymmetric())
```
Here we created a static diffusion edge with suitable variable and parameter names.
Similarly, we define the diffusion vertex with symbolic names.
```@example si
function _vertexf!(dv, v, esum, p, t)
    dv[1] = esum[1]
end
vertexf = ODEVertex(_vertexf!; sym=[:storage])
```


## Fundamental Symblic Indices
The default types for this access are the types [`VIndex`](@ref), [`EIndex`](@ref), [`VPIndex`](@ref) and [`EPIndex`](@ref).
Each of those symbolic indices consists of 2 elements: a reference to the network componen and a reference to the symbol within that component.
As such, `VIndex(2, :x)` refers to variable with symbolic name `:x` in vertex number 2.
`EPIndex(4, 2)` would refer to the *second* parameter of the edge component for the 4th edge.

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
sol(sol.t; idxs=VIndex(1, :storage))   # extract timeseries out ouf solution object
plot(sol; idxs=[VIndex(1, :storage), VIndex(5,:storage)]) # plot storage of two nodes
```

## Generate Symbolic Indices
Often, you need many individual symbolic indices. For that there are the helper methods [`vidxs`](@ref), [`eidxs`](@ref), [`vpidxs`](@ref) and [`epidxs`](@ref).
With the help of those methods you can generate arrays of symbolic indices:

```@example si
vidxs(nw, :, :storage) # get variable "storage" for all nodes
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
It essentially creates a new flat parameter array and fills it with the default parameter values define in the component.
The parameters in the `NWParameter` object can be accessed using the symbolic indices.
```@example si
p[EPIndex(5, :K)] = 2.0 # change the parameter K of the 5th edge
nothing #hide
```
Similarly, you can create a `NWState` object for the network `nw` using
```@example si
s = NWState(nw)
```
No default values were provided in the network components, so the state array is filled with `NaN`s.
```@example si
s[VIndex(:, :storage)] .= randn(5) # set the (initial) storage for alle nodes 
s #hide
```
For both `NWState` and `NWParameter` objects, the there is a more convenient way to access the variables and parameters.
```@example si
@assert s.v[1, :storage] == s[VIndex(1, :storage)] # s.v -> access vertex states
@assert s.e[1, :flow]    == s[EIndex(1, :flow)]    # s.e -> access edge states
@assert s.p.e[1,:K]      == p[EPIndex(1, :K)]      # s.p -> access parameters
```

The `NWState` and `NWParameter` objects are mutable, thus changing them will also change the underlying wrapped flat arrays.
You can allways access the flat representations by calling [`uflat`](@ref) and [`pflat`](@ref).

!!! note
    The `NWState` and `NWParameter` wrappers can be constructed from various objects.
    Fore example, within a callback you might construct `p = NWParameter(integrator)` to then change the parameters of the network within the callback.


## Observables
Some "States" arn't states in the sense of the ODE, but rather observables.
That means, they can be calculated explicitly from the state and parameters.
Those observables can be also accessed using symbolic indexing.
For example, the flows in the network don't show up in the state array but can be accessed in all the ways discussed above, for example

```@example si
plot(sol; idxs=eidxs(nw, :, :flow))
```
