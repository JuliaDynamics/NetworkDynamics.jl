# SDE Tutorial

 An `IJulia` [notebook](https://github.com/FHell/NetworkDynamics.jl/tree/master/examples) corresponding to this tutorial will be available on GitHub soon.

#### Topics covered in this tutorial include:
 * the swing equation
 * fixpoint search of ODE systems
 * constructing an SDE problem with NetworkDynamics.jl

## Example: Fluctuations in Power Grids

This tutorial explains the use of stochastic differential equations (SDE's) in NetworkDynamics.jl. As an example system we will simulate fluctuations in a simple power grid model.

The phase and frequency dynamics in power grids can be modeled by the swing equation. In the complex systems community this is also known as the Kuramoto model with inertia.

```math
\begin{aligned}
\dot \theta_i &= \omega_i \\
M_i \dot \omega_i &= P_i(t) - D_i \omega_i - \sum_{i=1}^N K_{ij} \sin(\theta_i - \theta_j)
\end{aligned}
```

Here, $M_i$ and $D_i$ are inertia and damping parameters, $K_{ij}$ is the power capacitiy of the line and $P_i(t)$ is the power infeed at the node. We assume that this power can be separated into a constant power setpoint $P^\circ_i$ and a stochastic fluctuation.

```math
\begin{aligned}
P_i(t) = P^\circ_i + \sigma_i \cdot \xi_i(t)
\end{aligned}
```

In this tutorial we assume the fluctuations $\xi_i(t)$ to be white Gaussian noise. Separating the deterministic and stochastic part, we can write this problem as a stochastic differential equation.

```math
\begin{aligned}
d\omega_i = \left[P_i - D_i \omega_i - \sum_{i=1}^N K_{ij} \sin(\theta_i - \theta_j)\right] dt + \sigma_i dW_i
\end{aligned}
```

Here, $dW_i = \xi_i dt$ is the infinitesimal increment of the Wiener process. Problems of this type can be numerically solved in Julia as described in the [SDE Tutorial](https://diffeq.sciml.ai/stable/tutorials/sde_example/) of `DifferentialEquations`.

## Implementing the Swing Equation

First we will implement the node and edge functions for the deterministic case without the fluctuations. We assume the inertia and damping parameters to be $M_i = 1.0$ and $D_i = 0.1$.

```@example SDEVertex
function swing_equation!(dv, v, e_s, e_d, P, t)
    dv[1] = v[2]
    dv[2] = P - 0.1 * v[2]
    for e in e_s
        dv[2] -= e[1]
    end
    for e in e_d
        dv[2] += e[1]
    end
end

function powerflow!(e, v_s, v_d, K, t)
    e[1] = K * sin(v_s[1] - v_d[1])
end
nothing # hide
```

## Contructing the Deterministic Dynamics

For the graph structure we will use a simple 4 node ring network.

```@example SDEVertex
using LightGraphs
g = watts_strogatz(4, 2, 0.)
nothing # hide
```

Then we can construct the `ODEFunction` of the deterministic system by using `network_dynamics()`.

```@example SDEVertex
using NetworkDynamics

swing_vertex = ODEVertex(f! = swing_equation!, dim = 2, sym=[:θ, :ω])
powerflow_edge = StaticEdge(f! = powerflow!, dim = 1)

nd = network_dynamics(swing_vertex, powerflow_edge, g)
nothing # hide
```

## Fixpoint Search

Now we need to define the dynamic parameters of vertices and edges. For simplicity we assume homogeneous capacities on the lines. For the nodes we assume half of them to be net producers ($P = 1.0$) and half of them to be net consumers ($P = -1.0$) of power.

```@example SDEVertex
K = 6.0
P = [1.,-1.,1.,-1.]
p = (P,K)
nothing # hide
```

We want to simulate fluctuations around an equilibrium state of our model system. Therefore, we need to find a fixpoint of the determinitic system which can be done by using the utility function `find_fixpoint()`. As an initial guess we take all variables equal to zero.

```@example SDEVertex
u0 = find_fixpoint(nd, p, zeros(8))

using DifferentialEquations
ode_prob = ODEProblem(nd, u0, (0.,500.), p)
ode_sol = solve(ode_prob)

using Plots, LaTeXStrings
plot(ode_sol, vars = syms_containing(nd, "ω"), ylims = (-1.0, 1.0), ylabel = L"\omega", legend = false)
```

We see that this is in fact a fixpoint solution. We will later use this as an initial condition for the numerical integration of the SDE system.

## Adding a Stochastic Layer

For adding the stochastic part of the dynamics we have to define a second graph layer. In our example, the fluctuations at different nodes are independent of each other. Therefore, we define a second graph with the same number of vertices but without any edges.

```@example SDEVertex
h =  SimpleGraph(4, 0)
nothing # hide
```
The dynamics at the nodes has to have the same dimension as in the deterministic case. In our example we only have fluctuations in the second variable.

```@example SDEVertex
function fluctuation!(dx, x, e_s, e_d, p, t)
    dx[1] = 0.0
    dx[2] = 0.05
end
nothing # hide
```

Now we can construct the dynamics of the second layer by using `network_dynamics()`. Since the graph structure of the stochastic layer has no edges we can take the edge function of the deterministic case as a placeholder.

```@example SDEVertex
fluctuation_vertex = ODEVertex(f! = fluctuation!, dim = 2)
nd_noise = network_dynamics(fluctuation_vertex, powerflow_edge, h)
nothing # hide
```

## Simulating the SDE

Finally, we can create an `SDEProblem` and solve it with `DifferentialEquations`.

```@example SDEVertex
sde_prob = SDEProblem(nd, nd_noise, u0, (0., 500.), p)
sde_sol = solve(sde_prob)
plot(sde_sol, vars = syms_containing(nd, "ω"), ylims = (-1.0, 1.0), ylabel = L"\omega", legend = false)
```

More details on SDE problems, e.g. how to include correlations or how to define an `EnsembleProblem`, can be found in the [documentation](https://diffeq.sciml.ai/stable/types/sde_types/) of `DifferentialEquations`.
