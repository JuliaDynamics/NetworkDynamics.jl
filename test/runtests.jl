using Test
@testset "NetworkDynamics Tests" begin
    # basic tests for individual modules
    include("NetworkStructures_test.jl")
    include("ComponentFunctions_test.jl")
    include("NetworkDiffEq_test.jl")
    include("checkbounds_test.jl")

    # complex tests of networks
    include("diffusion_test.jl")
    include("inhomogeneous_test.jl")
    include("autodiff_test.jl")
    include("massmatrix_test.jl")
end
