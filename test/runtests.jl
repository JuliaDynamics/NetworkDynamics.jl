using Test
using SafeTestsets

@testset "NetworkDynamics Tests" begin

@safetestset "utils test" begin include("utils_test.jl") end
@safetestset "construction test" begin include("construction_test.jl") end
@safetestset "Aggregation Tests" begin include("aggregators_test.jl") end
@safetestset "Symbolic Indexing Tests" begin include("symbolicindexing_test.jl") end

@safetestset "Diffusion test" begin include("diffusion_test.jl") end
@safetestset "inhomogeneous test" begin include("inhomogeneous_test.jl") end

end
