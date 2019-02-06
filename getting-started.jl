include("NetworkDynamics/src/NetworkDynamics.jl")

using LightGraphs


A = barabasi_albert(40,6)
L = laplacian_matrix(A, Float64)
dnd = network-dynamics.diffusive_network_dynamics(L, x -> -x^2)

x0 = rand(40)
dx0 = zeros(40)

dnd(dx0, x0, nothing, 0.)
