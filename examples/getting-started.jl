# You will find a step-by-step guide to this example in the docs and the
# corresponding jupyter notebook on our github repository.

using NetworkDynamics
using LightGraphs
using OrdinaryDiffEq
using Plots

### Defining a graph

N = 20 # number of nodes
k = 4  # average degree
g = barabasi_albert(N, k) # a little more exciting than a bare random graph


### Functions for edges and vertices

# the edge and vertex functions are mutating, hence
# by convention `!` is appended to their names

# e=edges, v_s=source-vertices, v_d=destination-vertices, p=parameters, t=time
function diffusionedge!(e, v_s, v_d, p, t)
    # usually e, v_s, v_d are arrays, hence we use the broadcasting operator .
    e .= v_s .- v_d
    nothing
end

# dv=derivative of vertex variables, v=vertex,
# e_s=source edges, e_d=destination edges, p=parameters, t=time
function diffusionvertex!(dv, v, e_s, e_d, p, t)
    # usually dv, v, e_s, e_d are arrays, hence we use the broadcasting operator .
    dv .= 0.
    # edges for which v is the source
    for e in e_s
        dv .-= e
    end
    # edges for which v is the destination
    for e in e_d
        dv .+= e
    end
    nothing
end

### Constructing the network dynamics

# ODEVertex and StaticEdge are structs that contain additional information on the function f!, such as number of variables of the internal function (dim), the symbols of those variables, and if a mass_matrix should be used
# VertexFunction/EdgeFunction is an abstract supertype for all vertex/edge function structs in NetworkDynamics

# signature of ODEVertex: (f!, dim, mass_matrix, sym) (mass_matrix and sym are optional arguments)
nd_diffusion_vertex = ODEVertex(f! = diffusionvertex!, dim = 1)

# signature of StaticEdge: (f!, dim, sym) (sym is optional)
nd_diffusion_edge = StaticEdge(f! = diffusionedge!, dim = 1)

# setting up the key constructor network_dynamics
# signature of network_dynamics: (vertices!, edges!, g; parallel = false)
# parameter parallel of type bool enables a parallel environment
# returned object nd of type ODEFunction is compatible with the solvers of OrdinaryDiffEq
nd = network_dynamics(nd_diffusion_vertex, nd_diffusion_edge, g)

### Simulation

x0 = randn(N) # random initial conditions
# ODEProblem is a struct of OrdinaryDiffEq.jl
# signature:(f::ODEFunction,u0,tspan, ...)
ode_prob = ODEProblem(nd, x0, (0., 4.))
# solve has signature: (prob::PDEProblem,alg::DiffEqBase.DEAlgorithm,args,kwargs)
sol = solve(ode_prob, Tsit5());

### Plotting

# vars=list of variables we want to plot, in this case we want to plot variables with symbol "v"
plot(sol, vars = syms_containing(nd, "v"))


### Bonus: Two independent diffusions with fancy symbols


# we will have two independent diffusions on the network, hence dim = 2 (x,ϕ)
nd_diffusion_vertex_2 = ODEVertex(f! = diffusionvertex!, dim = 2, sym = [:x, :ϕ])
nd_diffusion_edge_2 = StaticEdge(f! = diffusionedge!, dim = 2)
nd_2 = network_dynamics(nd_diffusion_vertex_2, nd_diffusion_edge_2, g)

x0_2 = vec(transpose([randn(N) randn(N).^2])) # x ~ normal distribution, ϕ ~ squared normal distribution
ode_prob_2 = ODEProblem(nd_2, x0_2, (0., 4.))
sol_2 = solve(ode_prob_2, Tsit5());



# try plotting the variables ϕ_i yourself. [Type \phi and press TAB]
plot(sol_2, vars = syms_containing(nd_2, "x"))
