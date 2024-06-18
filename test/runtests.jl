using Test
using SafeTestsets
using Pkg
using NetworkDynamics

@testset "NetworkDynamics Tests" begin
    @safetestset "utils test" begin include("utils_test.jl") end
    @safetestset "construction test" begin include("construction_test.jl") end
    @safetestset "Aggregation Tests" begin include("aggregators_test.jl") end
    @safetestset "Symbolic Indexing Tests" begin include("symbolicindexing_test.jl") end

    @safetestset "Diffusion test" begin include("diffusion_test.jl") end
    @safetestset "inhomogeneous test" begin include("inhomogeneous_test.jl") end
    @safetestset "massmatrix test" begin include("massmatrix_test.jl") end
    @safetestset "doctor test" begin include("doctor_test.jl") end
    @safetestset "initialization test" begin include("initialization_test.jl") end
end

@testset "Test Doc Examples" begin
    @info "Activate doc environment and test examples"
    Pkg.activate(joinpath(pkgdir(NetworkDynamics), "docs"))
    Pkg.develop(path=pkgdir(NetworkDynamics))
    Pkg.instantiate()

    examples = joinpath(pkgdir(NetworkDynamics), "docs", "examples")
    for file in readdir(examples; join=true)
        endswith(file, ".jl") || continue
        name = basename(file)
        @info "Test $name"
        eval(:(@safetestset $name begin include($file) end))
    end
end
