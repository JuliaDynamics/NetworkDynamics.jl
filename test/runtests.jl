# basic tests for indivitual modules
include("NetworkStructures_test.jl")
include("ComponentFunctions_test.jl")
include("nd_ODE_Static_test.jl")
include("nd_ODE_ODE_test.jl")

# complex tests of networks
include("diffusion_test.jl")
include("inhomogeneous_test.jl")
include("autodiff_test.jl")
