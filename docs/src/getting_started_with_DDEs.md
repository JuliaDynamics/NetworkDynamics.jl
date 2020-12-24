# Time-delayed network dynamics
 An `IJulia` [notebook](https://github.com/FHell/NetworkDynamics.jl/tree/master/examples) corresponding to this tutorial is available on GitHub.

 #### Topics covered in this tutorial include:
  * modeling time-delays with Delay Differential Equations (DDEs)
  * handling parameters for DDEs
  * time-delayed vertex variables
  * time-delayed coupling


## Delayed network diffusion

We revisit the introductory example on [Network diffusion](@ref getting_started), this time with time-delayed vertex dynamics. In this section our goal is not to model real physical phenomena but rather to learn how time-delays can be modeled in `NetworkDynamics.jl`.

Let $g$ be a graph with $N$ nodes and adjacency matrix $A$. Let $v = (v_1, \dots, v_n)$ be a vector of (abstract) temperatures or concentrations at each node $i = 1, \dots, N$. The dynamics of state $v_i$ at time $t$ are determined by its past value $v_i(t-\tau)$, where $\tau$ is a constant time lag and its difference with its neighbors with coupling strength $\sigma$.

```math
\begin{aligned}
\dot v_i(t) = - v_i(t-\tau) - \sigma * \sum_{i=1}^N A_{ij} (v_i(t) - v_j(t))
\end{aligned}
```

## Modeling diffusion with NetworkDynamics.jl

The sum on the right hand side plays the role of a (discrete) gradient. If the temperature at node $i$ is higher than at its neighboring node $j$, it will decrease along that edge.

The coupling function is the same as in the [previous example](@ref getting_started)
```@example DDEVertex
function diffusionedge!(e, v_s, v_d, p, t)
   e .= .1 * (v_s .- v_d)
   nothing
end
nothing # hide
```
However, the internal vertex dynamics are now determined by the time-delayed equation $\dot v_i(t) = - v_i(t-\tau)$ and are described in the vertex function with help of the history array $h_v$ containing the past values of the vertex variables.
```@example DDEVertex

function delaydiffusionvertex!(dv, v, edges, h_v, p, t)
   # h_v is the history array of the vertex variables
   dv .= -h_v
   sum_coupling!(dv, edges)
   nothing
end
nothing # hide
```
The Watts-Strogatz algorithm constructs a regular ring network with $N$ nodes connected to $k$ neighbors and a probability $p$ of rewiring links.  Since $p=0$ no rewiring happens and `g` is a simple ring network.

```@example DDEVertex
using LightGraphs

### Defining a graph

N = 10 # number of nodes
k = 4  # average degree
g = watts_strogatz(N, k, 0.) # ring-network
nothing # hide
```
While the `EdgeFunction` is constructed via `StaticEdge` just as before, the vertex functions is wrapped in `DDEVertex` to account for the time-lag in the internal dynamics. `network_dynamics` then returns a `DDEFunction` compatible with the solvers of DifferentialEquations.jl. Combining those into a `DDEFunction` via `network_dynamics` is then straightforward.

```@example DDEVertex
using NetworkDynamics

nd_diffusion_vertex = DDEVertex(f! = delaydiffusionvertex!, dim = 1)
nd_diffusion_edge = StaticEdge(f! = diffusionedge!, dim = 1)

nd = network_dynamics(nd_diffusion_vertex, nd_diffusion_edge, g)
nothing # hide
```


# Simulation
When simulating a time-delayed system the initial conditions have to be specified on the whole interval $[-\tau, 0 ]$. This is usually done by specifying initial conditions at `t=0` and a history function `h` for extrapolating these values to time $t= - \tau$. For simplicity `h` is often chosen constant. When solving a DDE with DifferentialEquations.jl `h` is later overloaded to interpolate between time-steps, for details [see the docs](https://diffeq.sciml.ai/stable/tutorials/dde_example/).

```@example DDEVertex

const x0 = randn(N) # random initial conditions
tspan = (0., 10.) # time interval

# history function keeps the initial conditions constant and is in-place to save allocations
h(out, p, t) = (out .= x0)
nothing # hide
```
We still have to specify the *constant* time-lag $\tau$. At the moment this is only possible by using the ND-specific [tuple syntax](@ref parameters). For DDEs these tuples have to contain a third argument specifying the delay time.

```@example DDEVertex
# parameters p = (vertex parameters, edge parameters, delay time)
p = (nothing, nothing, 1.)
nothing # hide
```
 We solve the problem on the time interval $[0, 10]$ with the `Tsit5()` algorithm, recommended for solving non-stiff problems. The `MethodOfSteps(..)` translates an `OrdinaryDiffEq.jl`  solver into an efficient method for delay differential equations. `DDEProblem` addtional accepts the keyword `constant_lags` that can be useful in some situation, see [their docs](https://diffeq.sciml.ai/stable/tutorials/dde_example/) for details.

```@example DDEVertex
using DelayDiffEq, Plots

dde_prob = DDEProblem(nd, x0, h, tspan, p)
sol = solve(dde_prob, MethodOfSteps(Tsit5()))
plot(sol, vars = syms_containing(nd, "v"))
```

## Bonus: Two independent diffusions

In this extension of the first example, there are two independent diffusions on the same network with variables $x$ and $\phi$ - hence the dimension is set to `dim=2`.

```@example DDEVertex
nd_diffusion_vertex_2 = DDEVertex(f! = delaydiffusionvertex!, dim = 2, sym = [:x, :ϕ])
nd_diffusion_edge_2 = StaticEdge(f! = diffusionedge!, dim = 2)
nd_2 = network_dynamics(nd_diffusion_vertex_2, nd_diffusion_edge_2, g)
nothing # hide
```
The initial conditions are sampled from (squared) normal distributinos such that the first $N$ values correspond to variable `x` and the values with indices from $N+1$ to $2N$ belong to variable `ϕ`, where $x \sim \mathcal{N}(0,1)$; $ϕ \sim \mathcal{N}(0,1)^2$.

```@example DDEVertex
const x0_2 = Array{Float64,1}(vec([randn(N).-10 randn(N).^2]')) # x ~ \mathcal{N}(0,1); ϕ ~ \mathcal{N}(0,1)^2

h_2(out, p, t) = (out .= x0_2)

p = (nothing, nothing, 1.) # p = (vertexparameters, edgeparameters, delaytime)
nothing # hide
```
Now we can define the `DDEProblem` and solve it.

```@example DDEVertex
dde_prob_2 = DDEProblem(nd_2, x0_2, h_2, tspan, p)
sol_2 = solve(dde_prob_2, MethodOfSteps(Tsit5()));
plot(sol_2, legend=false)
```

## Kuramoto model with delay

A common modification of the [Kuramoto model](https://en.wikipedia.org/wiki/Kuramoto_model) is to include a time-lag in the coupling function. In neuroscience such a coupling is used to account for transmission delays along synapses connecting different neurons.

Unlike in the diffusion example, the edges depend on past values of the vertex variables . For this reason the edgefunction has the history arrays of the destination and source vertices as arguments.

```@example DDEVertex
function kuramoto_delay_edge!(e, v_s, v_d, h_v_s, h_v_d, p, t)
    e[1] = p * sin(v_s[1] - h_v_d[1])
    nothing
end
nothing # hide
```

The dynamics of vertices in the Kuramoto model are determined by a constant rotation frequency `p`. Contributions of the edges are summed up according to their destinations.

```@example DDEVertex
function kuramoto_vertex!(dv, v, edges, p, t)
    dv[1] = p
    for e in edges
        dv .+= e
    end
    nothing
end
nothing # hide
```

In this case the vertex function does not depend on any past values and is simply constructed as an `ODEVertex`. Since the edges depend on the history of the vertex variables their calling signature has changed and they are handed to the delay-aware constructor `StaticDelayEdge`. *Static* refers to the fact that they don't have any internal dynamics.

```@example DDEVertex
kdedge! = StaticDelayEdge(f! = kuramoto_delay_edge!, dim = 2)
kdvertex! = ODEVertex(f! = kuramoto_vertex!, dim = 1)

nd! = network_dynamics(kdvertex!, kdedge!, g)
nothing # hide
```

The remaining steps for simulating the system are the same as above.

```@example DDEVertex
const x0_3 = randn(N) # random initial conditions
h_3(out, p, t) = (out .= x0_3)
# p = (vertexparameters, edgeparameters, delaytime)
ω = randn(N)
ω .-= sum(ω)/N # center at 0
p = (ω, 2., 1.)
tspan = (0.,10.)

dde_prob = DDEProblem(nd!, x0_3, h_3, tspan, p)

sol = solve(dde_prob, MethodOfSteps(Tsit5()));

### Plotting

plot(sol, vars = syms_containing(nd, "v"), legend=false)

```

## Delay differential equations with algebraic constraints

The vertex dynamics may in principle contain algebraic equations in mass matrix form. For an experimental test case have a look at the last section of this [example on github](https://github.com/FHell/NetworkDynamics.jl/blob/master/examples/getting-started-with-DDEs.jl).
