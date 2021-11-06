# using Flux, Optim
using DiffEqFlux
using Flux: ADAM
#using Optim: BFGS
using DiffEqSensitivity
using NetworkDynamics
using Graphs
using OrdinaryDiffEq
using Plots


### Defining a graph

N = 20 # number of nodes
k = 4  # average degree
g = barabasi_albert(N, k) # a little more exciting than a bare random graph


### Functions for edges and vertices

function diffusionedge!(e, v_s, v_d, σ, t)
    e .= σ .* (v_s .- v_d)
    nothing
end

function diffusionvertex!(dv, v, e_s, e_d, p, t)
    dv .= 0.0
    for e in e_s
        dv .-= e
    end
    for e in e_d
        dv .+= e
    end
    nothing
end

### Constructing the network dynamics

nd_diffusion_vertex = ODEVertex(; f=diffusionvertex!, dim=1)
nd_diffusion_edge = StaticEdge(; f=diffusionedge!, dim=1)
nd! = network_dynamics(nd_diffusion_vertex, nd_diffusion_edge, g)

### Simulation

x0 = randn(N) # random initial conditions
# These are the parameters for the edges, i.e. the heterogeneous diffusion constants
σ = rand(ne(g))
tspan = (0.0, 4.0)
ode_prob = ODEProblem(nd, x0, tspan, (nothing, σ))
sol = solve(ode_prob, Tsit5())

plot(sol)


# DiffEqFlux

# Some AD algorithms don't like parameters which are tuples.
# In order to improve compatibility with AD algorithms it would be desirable to have
# a version of network_dynamics that utilizes only arrays of parameters, without the
# tuple syntax. At the moment the workaround is to use a wrapper function.

function nd_wrapper!(dx, x, σ, t)
    nd!(dx, x, (nothing, σ), t)
end


probflux = ODEProblem(nd_wrapper!, x0, tspan, σ)


# #### The prediction function
#
# The function `predict` integrates the system forward in time.
# Sensitivity analysis refers to computing the sensitivity of the solution with respect to # changes in the parameters, i.e. partial derivatives.
#
# The choice of sensitiviy algorithm is crucial, because it determines which method will
# be used for AD. For details consult the DiffEqSensitivity manual.  I recommend using
# InterpolatingAdjoint or BacksolveAdjoint. Both algorithms work only on ODEs.
# Tracker is more generally applicable but is tailored towards out-of-place problems and
# is slow on in-place problems. Zygote fails.


function predict(p)
    ## default sensealg is InterpolatingAdjoint
    solve(probflux, Tsit5(); p=p, saveat=tspan[1]:0.01:tspan[end], sensealg=ForwardDiffSensitivity())

end


function loss(p)
    pred = predict(p)
    # converge to 0 as fast as possible
    loss = sqrt(sum(abs2, pred))# + sqrt(sum(abs2,p))
    loss, pred
end

cb = function (p, l, pred) # callback function to observe training
    println(" loss = ", l)
    display(plot(pred))
    return false
end

cb(σ, loss(σ)...)

# A crucial thing to realize is that sciml_train works best with Arrays of parameters

# We optimize for optimal local diffusion constants
res = DiffEqFlux.sciml_train(loss, σ, ADAM(0.5); cb=cb, maxiters=20)

# res = DiffEqFlux.sciml_train(loss, σ, BFGS(), cb = cb)
