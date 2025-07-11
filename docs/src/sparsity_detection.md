# Sparsity Detection

NetworkDynamics.jl can automatically detect and exploit the sparsity structure of the Jacobian matrix to significantly improve the performance of ODE solvers. This feature uses [SparseConnectivityTracer.jl](https://github.com/adrhill/SparseConnectivityTracer.jl) to analyze the network's dynamics and create a sparse Jacobian prototype that modern solvers can use for more efficient linear algebra operations.

The sparsity detection is particularly beneficial for:
- Large networks where the Jacobian matrix is sparse
- Stiff systems that require implicit solvers
- Networks with complex component interactions
- Components with conditional statements that complicate automatic differentiation

## Core Function

The main interface is the [`get_jac_prototype`](@ref) function:

```@docs
get_jac_prototype
```

The sparsity pattern can be passed to ODE solvers to improve performance:

```julia
f_ode = ODEFunction(nw; jac_prototype=get_jac_prototype(nw))
prob = ODEProblem(f_ode, x0, (0.0, 1.0), p0)
sol = solve(prob, Rodas5P())
```

Alternatively, you can store the sparsity pattern directly in the network:

```julia
set_jac_prototype!(nw; kwargs_for_get_jac_prototype...)
prob = ODEProblem(nw, x0, (0.0, 1.0), p0)  # automatically uses stored prototype
```

## Example: Handling Conditional Statements

A key feature of NetworkDynamics.jl's sparsity detection is the ability to handle conditional statements in component functions. This is particularly useful for ModelingToolkit-based components that use `ifelse` statements.

The conditional statements will be resolved in favor of a "global" sparsity pattern by
replacing them temporarily with `trueblock + falseblock` which is then inferable by
SparseConnectivityTracer.jl.

!!! details "Setup code"
    ```@example sparsity
    using NetworkDynamics, ModelingToolkit, Graphs
    using SparseArrays, OrdinaryDiffEq
    using ModelingToolkit: D_nounits as Dt, t_nounits as t
    nothing #hide
    ```

```@example sparsity
# Define a component with conditional logic
@mtkmodel ValveModel begin
    @variables begin
        p_src(t), [description="source pressure"]
        p_dst(t), [description="destination pressure"]
        q(t), [description="flow through valve"]
    end
    @parameters begin
        K=1, [description="conductance"]
        active=1, [description="valve state"]
    end
    @equations begin
        q ~ ifelse(active > 0, K * (p_src - p_dst), 0)
    end
end

@mtkmodel NodeModel begin
    @variables begin
        p(t)=1, [description="pressure"]
        q_nw(t), [description="network flow"]
    end
    @parameters begin
        C=1, [description="capacitance"]
        q_ext, [description="external flow"]
    end
    @equations begin
        C*Dt(p) ~ q_ext + q_nw
    end
end
```

```@example sparsity
# Create network
@named valve = ValveModel()
@named node = NodeModel()

g = wheel_graph(10)
v = VertexModel(node, [:q_nw], [:p])
e = EdgeModel(valve, [:p_src], [:p_dst], AntiSymmetric([:q]))

nw = Network(g, v, e)
```

```@example sparsity
# This will fail due to conditional statements
try
    get_jac_prototype(nw)
catch
    println("Error: Sparsity detection failed due to conditional statements")
end
```

```@example sparsity
# This works by removing conditionals
jac_prototype = get_jac_prototype(nw; remove_conditions=true)

# Store the prototype directly in the network
set_jac_prototype!(nw, jac_prototype)
```

## Performance Benefits

Using sparsity detection can significantly improve solver performance, especially for large networks and stiff systems:

```@example sparsity
using OrdinaryDiffEqRosenbrock, Chairmarks

# Create a large sparse network for benchmarking
g_large = grid([20, 20])  # 400 nodes in a 2D grid (very sparse)
nw_large = Network(g_large, v, e)

# Setup initial conditions and parameters
using Random # hide
Random.seed!(42) # hide
s0 = NWState(nw_large)
s0.v[:, :p] .= randn(400)  # random initial pressures

p0 = NWParameter(nw_large)
p0.v[:, :q_ext] .= randn(400)  # small external flow

nothing #hide
```

The network is now ready for benchmarking. Let's first time the solution without sparsity detection:

```@example sparsity
# Without sparsity detection (dense Jacobian)
prob_dense = ODEProblem(nw_large, uflat(s0), (0.0, 1.0), pflat(p0))
@b solve($prob_dense, Rodas5P()) seconds=1
```

Now let's enable sparsity detection:
```@example sparsity
jac = get_jac_prototype(nw_large; remove_conditions=true)
```
The pattern already shows that the Jacobian is really sparse due to the sparse network connections.

```@example sparsity
set_jac_prototype!(nw_large, jac)
```

Now we can benchmark the sparse version:
```@example sparsity
# Solve with sparsity detection
prob_sparse = ODEProblem(nw_large, uflat(s0), (0.0, 1.0), pflat(p0))
@b solve($prob_sparse, Rodas5P()) seconds=1
```

For this network, we see a substantial speedup due to the sparse solver!

## Troubleshooting

**Sparsity detection fails with conditional statements:**
- Use `remove_conditions=true` to handle `ifelse` statements in MTK components
- For specific problematic components, pass a vector of indices: `remove_conditions=[EIndex(1), VIndex(2)]`

**Detection fails for complex components:**
- Use `dense=true` to treat all components as dense (fallback option)
- For specific components, use `dense=[EIndex(1)]` to treat only those components as dense

**Performance doesn't improve:**
- Sparsity detection is most beneficial for large networks (>50 nodes) with sparse connectivity
- Dense networks or small systems may not see significant speedup
- Ensure you're using a solver that can exploit sparsity (e.g., `Rodas5P`, `FBDF`)

The sparsity detection feature requires the `SparseConnectivityTracer.jl` package, which is automatically loaded as a conditional dependency when needed.
