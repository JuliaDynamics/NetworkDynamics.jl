# using NetworkDynamics

using LightGraphs
using Reexport

include("/Users/lalalulu/Desktop/PIK/github/NetworkDynamics.jl/src/Utilities.jl")
@reexport using .Utilities

include("/Users/lalalulu/Desktop/PIK/github/NetworkDynamics.jl/Jacobians/jac_structs.jl")
@reexport using .jac_structs

include("/Users/lalalulu/Desktop/PIK/github/NetworkDynamics.jl/Jacobians/jac_functions.jl")
@reexport using .jac_functions


function diffusionedge!(e, v_s, v_d, p, t)
    # usually e, v_s, v_d are arrays, hence we use the broadcasting operator .
    e .= v_s .- v_d
    nothing
end

function diffusionvertex!(dv, v, edges, p, t)
    # usually dv, v, edges are arrays, hence we use the broadcasting operator .
    dv .= 0.
    for e in edges
        dv .+= e
    end
    nothing
end

nd_diffusion_vertex = ODEVertex(f! = diffusionvertex!, dim = 1)

nd_diffusion_edge = StaticEdge(f! = diffusionedge!, dim = 1)

nd = network_dynamics(nd_diffusion_vertex, nd_diffusion_edge, g)

### Simulation

x0 = randn(N) # random initial conditions
ode_prob = ODEProblem(nd, x0, (0., 4.))
sol = solve(ode_prob, Tsit5());

function VertexJacobian!(J, v, p, t)
   J = internal_jacobian(v, p, t)
   # get_vertex_jacobian(gd, i) = something(v,p,t)
   # gd.v_jac_array[i] = something(v,p,t)
end
