# call from shell: env JULIA_NUM_THREADS=1 julia mutltithreading_benchmarks.jl

using NetworkDynamics
using LightGraphs
using BenchmarkTools

### Defining a graph

N = 1000 # number of nodes
k = 20  # average degree
g = barabasi_albert(N, k) # a little more exciting than a bare random graph


### Functions for edges and vertices

@inline Base.@propagate_inbounds function diffusionedge!(e, v_s, v_d, p, t)
    # usually e, v_s, v_d are arrays, hence we use the broadcasting operator .
    e[1] = p * (v_s[1] - v_d[1])
    nothing
end

@inline Base.@propagate_inbounds function diffusionvertex!(dv, v, e_s, e_d, p, t)
    # usually v, e_s, e_d are arrays, hence we use the broadcasting operator .
    dv[1] = p
    # edges for which v is the source
    for e in e_s
        dv[1] -= e[1]
    end
    # edges for which v is the destination
    for e in e_d
        dv[1] += e[1]
    end
    nothing
end

### Constructing the network dynamics

nd_diffusion_vertex = ODEVertex(f! = diffusionvertex!, dim = 1)
nd_diffusion_edge = StaticEdge(f! = diffusionedge!, dim = 1)


### Benchmarking
println("Number of Threads: ", haskey(ENV, "JULIA_NUM_THREADS") ? ENV["JULIA_NUM_THREADS"] : "1")
println("\nBenchmarking ODE_Static...")
p  = (collect(1:nv(g))./nv(g) .- 0.5, .5 .* ones(ne(g)))
println("\nsingle-threaded...")
nd = network_dynamics(nd_diffusion_vertex, nd_diffusion_edge, g, parallel=false)
x0 = randn(N)
display(@benchmark $nd($x0, $x0, $p , 0,))
println("\nparallel...")
nd = network_dynamics(nd_diffusion_vertex, nd_diffusion_edge, g, parallel=true)
display(@benchmark $nd($x0, $x0, $p , 0,))
#=
### Simulation
using OrdinaryDiffEq

x0 = randn(N)

ode_prob = ODEProblem(nd, x0, (0., 4.), p)
sol = solve(ode_prob, Tsit5());

### Plotting
using Plots
plot(sol, vars = syms_containing(nd, "v"))
=#
### Redefinition

@inline Base.@propagate_inbounds function diffusion_dedge!(de, e, v_s, v_d, p, t)
    de[1] = 100. * (p * sin(v_s[1] - v_d[1]) - e[1])
    nothing
end

nd_diffusion_dedge = ODEEdge(f! = diffusion_dedge!, dim = 1)



p  = (collect(1:nv(g))./nv(g) .- 0.5, 5 .* ones(ne(g)))

### Benchmarking
println("\nBenchmarking ODE_ODE...")
println("\nsingle-threaded...")
nd_ode = network_dynamics(nd_diffusion_vertex, nd_diffusion_dedge, g, parallel=false)
x0 = randn(nv(g) + ne(g))
dx0 = randn(nv(g) + ne(g))
display(@benchmark $nd_ode($x0, $dx0, $p , 0,))
println("\nparallel...")
nd_ode = network_dynamics(nd_diffusion_vertex, nd_diffusion_dedge, g, parallel=true)
display(@benchmark $nd_ode($x0, $dx0, $p , 0,))


println("")
#=
### Simulation

x0 = randn(nv(g) + ne(g))

ode_ode_prob = ODEProblem(nd_ode, x0, (0., 2.), p)
sol_ode = solve(ode_ode_prob, Tsit5());

### Plotting

plot(sol_ode, vars = syms_containing(nd, "v"))

=#
