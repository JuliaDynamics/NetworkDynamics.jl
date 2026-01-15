using Test
using SafeTestsets
using Pkg
using NetworkDynamics
using SciMLBase
using InteractiveUtils
using CUDA
using KernelAbstractions
using Adapt
using NetworkDynamics: iscudacompatible, NaiveAggregator
using Aqua
using ExplicitImports

(isinteractive() ? includet : include)(joinpath(pkgdir(NetworkDynamics, "test", "testutils.jl")))

BUILDKITE = haskey(ENV, "BUILDKITE")
BUILDKITE && @test CUDA.functional() # fail early in buildkite if cuda is not available

@testset "NetworkDynamics Tests" begin
    BUILDKITE || @testset "Package Quality Tests" begin
        # print_explicit_imports(NetworkDynamics)
        @test check_no_implicit_imports(NetworkDynamics) === nothing
        @test check_no_stale_explicit_imports(NetworkDynamics, ignore=(:Symbolics,)) === nothing
        Aqua.test_all(NetworkDynamics;
            ambiguities=false,
            stale_deps=true,
            persistent_tasks=false)
        @test_broken isempty(Docs.undocumented_names(NetworkDynamics))
    end

    # skip non-GPU tests on buildkite
    if !BUILDKITE
        @safetestset "utils test" begin include("utils_test.jl") end

        NetworkDynamics.CHECK_COMPONENT[] = false
        @safetestset "construction test" begin include("construction_test.jl") end
        @safetestset "Aggregation Tests" begin include("aggregators_test.jl") end
        @safetestset "massmatrix test" begin include("massmatrix_test.jl") end
        NetworkDynamics.CHECK_COMPONENT[] = true

        @safetestset "initialization test" begin include("initialization_test.jl") end
        @safetestset "Callbacks test" begin include("callbacks_test.jl") end
        @safetestset "Metadata test" begin include("metadata_test.jl") end
        @safetestset "Linear Stability test" begin include("linear_stability_test.jl") end
        @safetestset "Show-methods test" begin include("show_test.jl") end
        @safetestset "Spinners test" begin include("spinners_test.jl") end
        @safetestset "sparsity test" begin include("sparsity_test.jl") end

        @safetestset "MTK extension test" begin include("MTK_test.jl") end

        # check on the precompile files
        @safetestset "Precompile workload" begin include("../src/precompile_workload.jl") end
        @safetestset "MTK precompile workload" begin include("../ext/MTKExt_precomp_workload.jl") end

        @safetestset "AD test" begin include("AD_test.jl") end
        @safetestset "doctor test" begin include("doctor_test.jl") end
    end

    @safetestset "Symbolic Indexing Tests" begin include("symbolicindexing_test.jl") end
    @safetestset "external input test" begin include("external_inputs_test.jl") end
    @safetestset "Loopback Connection test" begin include("loopback_test.jl") end

    @safetestset "Diffusion test" begin include("diffusion_test.jl") end
    @safetestset "inhomogeneous test" begin include("inhomogeneous_test.jl") end

    if CUDA.functional()
        @safetestset "GPU test" begin include("GPU_test.jl") end
    end
end

@testset "Test Doc Examples" begin
    examples = joinpath(pkgdir(NetworkDynamics), "docs", "examples")
    for file in readdir(examples; join=true)
        endswith(file, ".jl") || continue
        name = basename(file)

        @info "Test $name"
        eval(:(@safetestset $name begin include($file) end))
    end
end

if !CUDA.functional()
    @test gethostname() != "hw-g14" # on this pc, CUDA *should* be available
    @warn "Skipped all CUDA tests because no device is available."
end
