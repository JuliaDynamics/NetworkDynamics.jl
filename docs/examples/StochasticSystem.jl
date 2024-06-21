#=
# SDE Tutorial

This example can be dowloaded as a normal Julia script [here](@__NAME__.jl). #md

#### Topics covered in this tutorial include:

  - the swing equation
  - fixpoint search of ODE systems
  - constructing an SDE problem with NetworkDynamics.jl

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

First we will implement the node and edge functions for the deterministic case
without the fluctuations. We set the defaults for the inertia and damping
parameters to be $M_i = 1.0$ and $D_i = 0.1$. The default coupling strength is $K=6$.
=#
using NetworkDynamics, Graphs
using StochasticDiffEq, OrdinaryDiffEq
using Plots, LaTeXStrings

function swing_equation!(dv, v, esum, (M, P, D), t)
    dv[1] = v[2]
    dv[2] = 1/M *(P - D * v[2] + esum[1])
    nothing
end
swing_vertex = ODEVertex(swing_equation!; sym=[:θ, :ω], psym=[:M=>1, :P, :D=>0.1])
#-

function powerflow!(e, v_s, v_d, (K,), t)
    e[1] = K * sin(v_s[1] - v_d[1])
end
powerflow_edge = StaticEdge(powerflow!; dim=1, psym=[:K=>6], coupling=AntiSymmetric())

#=
## Contructing the Deterministic Dynamics

For the graph structure we will use a simple 4 node ring network.
=#

g = watts_strogatz(4, 2, 0.0)
nothing #hide #md

#=
Then we can construct the `Network` of the deterministic system.
=#

nd = Network(g, swing_vertex, powerflow_edge)

#=
## Fixpoint Search

Now we need to define the dynamic parameters of vertices and edges.
For that, we start by creating the `NWParameter` object, pre filling it with the default values:
=#

p = NWParameter(nd)

#=
Most of the parameters are allready set
For the nodes we assume half of them to be net producers ($P = 1.0$) and half of them to be net consumers ($P = -1.0$) of power.
=#
p.v[1:4, :P] = [1, -1, 1, -1]
nothing #hide #md

#=
We want to simulate fluctuations around an equilibrium state of our model
system. Therefore, we need to find a fixpoint of the determinitic system which
can be done by using the utility function `find_fixpoint()`. As an initial guess
we take all variables equal to zero.
=#

u0 = find_fixpoint(nd, p)

ode_prob = ODEProblem(nd, uflat(u0), (0.0, 500.0), pflat(p))
Main.test_execution_styles(ode_prob) # testing all ex styles #src
ode_sol = solve(ode_prob, Tsit5())

plot(ode_sol; idxs=vidxs(nd,:,:ω), ylims=(-1.0, 1.0), ylabel=L"\omega", legend=false, fmt=:png)

#=
We see that this is in fact a fixpoint solution. We will later use this as an initial condition for the numerical integration of the SDE system.

## Adding a Stochastic Layer

For adding the stochastic part of the dynamics we have to define a second graph layer. In our example, the fluctuations at different nodes are independent of each other. Therefore, we define a second graph with the same number of vertices but without any edges.
=#

h = SimpleGraph(4, 0)
nothing #hide #md

#=
The dynamics at the nodes has to have the same dimension as in the deterministic case. In our example we only have fluctuations in the second variable.
=#

function fluctuation!(dx, x, edges, p, t)
    dx[1] = 0.0
    dx[2] = 0.05
end
nothing #hide #md

#=
Now we can construct the dynamics of the second layer by using `network_dynamics()`. Since the graph structure of the stochastic layer has no edges we can take the edge function of the deterministic case as a placeholder.
=#

fluctuation_vertex = ODEVertex(fluctuation!; dim=2)
nd_noise = Network(h, fluctuation_vertex, NetworkDynamics.EdgeFunction[])
nothing #hide #md

#=
## Simulating the SDE

Finally, we can create an `SDEProblem` and solve it with `DifferentialEquations`.
=#

sde_prob = SDEProblem(nd, nd_noise, uflat(u0), (0.0, 500.0), pflat(p))
sde_sol = solve(sde_prob, SOSRA())
plot(sde_sol; idxs=vidxs(nd,:,:ω), ylims=(-1.0, 1.0), ylabel=L"\omega", legend=false, fmt=:png)

#=
More details on SDE problems, e.g. how to include correlations or how to define an `EnsembleProblem`, can be found in the [documentation](https://diffeq.sciml.ai/stable/types/sde_types/) of `DifferentialEquations`.
=#
