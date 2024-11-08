using Test
using SafeTestsets
using Pkg
using NetworkDynamics
using SciMLBase
using InteractiveUtils
using CUDA
using KernelAbstractions
using Adapt
using NetworkDynamics: iscudacompatible
using Aqua
using ExplicitImports

isinteractive() ? includet("testutils.jl") : include("testutils.jl")

@testset "NetworkDynamics Tests" begin
    @testset "Package Quality Tests" begin
        # print_explicit_imports(NetworkDynamics)
        @test check_no_implicit_imports(NetworkDynamics) === nothing
        @test check_no_stale_explicit_imports(NetworkDynamics) === nothing
        Aqua.test_all(NetworkDynamics;
            ambiguities=false,
            persistent_tasks=false)
        @test_broken isempty(Docs.undocumented_names(NetworkDynamics))
    end

    @safetestset "utils test" begin include("utils_test.jl") end

    NetworkDynamics.CHECK_COMPONENT[] = false
    @safetestset "construction test" begin include("construction_test.jl") end
    @safetestset "Aggregation Tests" begin include("aggregators_test.jl") end
    @safetestset "Symbolic Indexing Tests" begin include("symbolicindexing_test.jl") end
    @safetestset "massmatrix test" begin include("massmatrix_test.jl") end
    NetworkDynamics.CHECK_COMPONENT[] = true

    @safetestset "doctor test" begin include("doctor_test.jl") end

    @safetestset "Diffusion test" begin include("diffusion_test.jl") end
    @safetestset "inhomogeneous test" begin include("inhomogeneous_test.jl") end
    @safetestset "initialization test" begin include("initialization_test.jl") end

    @safetestset "AD test" begin include("AD_test.jl") end

    if CUDA.functional()
        @safetestset "GPU test" begin include("GPU_test.jl") end
    end

    @safetestset "MTK extension test" begin include("MTK_test.jl") end

    # check on the precompile files
    # @safetestset "Precompile workload" begin include("../src/precompile_workload.jl") end
    # @safetestset "MTK precompile workload" begin include("../ext/precompile_workload.jl") end
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
        if name == "kuramoto_delay.jl"
            @test_broken false
            continue
        end
        eval(:(@safetestset $name begin include($file) end))
    end
end

if !CUDA.functional()
    @test gethostname() != "hw-g14" # on this pc, CUDA *should* be available
    @warn "Skipped all CUDA tests because no device is available."
end
