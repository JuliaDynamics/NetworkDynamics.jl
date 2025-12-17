# Sparsity Detection

NetworkDynamics.jl can automatically detect and exploit the sparsity structure of the Jacobian matrix to significantly improve the performance of ODE solvers. This feature uses [SparseConnectivityTracer.jl](https://github.com/adrhill/SparseConnectivityTracer.jl) to analyze the network's dynamics and create a sparse Jacobian prototype that modern solvers can use for more efficient linear algebra operations.

The sparsity detection is particularly beneficial for:
- Large networks where the Jacobian matrix is sparse
- Stiff systems that require implicit solvers
- Networks with complex component interactions
- Components with conditional statements that complicate automatic differentiation

## Core Function

The main interface is the [`get_jac_prototype`](@ref) function, which takes a `Network` object as an argument and returns a sparse boolean matrix containing the sparsity pattern.

You can store the sparsity pattern directly in the network, which will then be
picked up by the `ODEProblem` and `ODEFunction` constructors.

```julia
set_jac_prototype!(nw; kwargs_for_get_jac_prototype...)
prob = ODEProblem(nw, x0, (0.0, 1.0), p0)  # automatically uses stored prototype
```

The `get_jac_prototype` function will operate on batches of identical components.
If the user provided component functions are **not SCT compatible**, it'll first try to
resolve `if..else..end` statements in MTK-generated code and fall back to dense component functions.
Even the dense component fallback can lead to a substantial speedup because most of the sparsity stems
from the Network sparsity rather than the component sparsity.

## Example: Handling Conditional Statements

A key feature of NetworkDynamics.jl's sparsity detection is the ability to automatically handle conditional statements in MTK component functions.

The conditional `if...else...end` statements generated in the codegen phase will be replaced by equivalent `ifelse(..,..,..)` statements which can be handled by SCT.

!!! details "Setup code"
    ```@example sparsity
    using NetworkDynamics, ModelingToolkit, Graphs
    using SparseArrays, OrdinaryDiffEqRosenbrock, OrdinaryDiffEqNonlinearSolve, NonlinearSolve
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
nothing # hide
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
# This works by removing conditionals
jac_prototype = get_jac_prototype(nw)

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
@b solve($prob_dense, Rodas5P()) evals=1
```

Now let's enable sparsity detection:
```@example sparsity
jac = get_jac_prototype(nw_large)
```
The pattern already shows that the Jacobian is really sparse due to the sparse network connections.

```@example sparsity
set_jac_prototype!(nw_large, jac)
```

Now we can benchmark the sparse version:
```@example sparsity
# Solve with sparsity detection
prob_sparse = ODEProblem(nw_large, uflat(s0), (0.0, 1.0), pflat(p0))
@b solve($prob_sparse, Rodas5P(linsolve=KLUFactorization())) evals=1
```

For this network, we see a substantial speedup due to the sparse solver!
