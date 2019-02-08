include("src/NetworkDynamics.jl")

using LightGraphs


A = barabasi_albert(40,6)
dnd = NetworkDynamics.diffusive_network_dynamics(A, x -> -x^2)

x0 = rand(40)
dx0 = zeros(40)

dnd(dx0, x0, nothing, 0.)
println(dnd.L)
