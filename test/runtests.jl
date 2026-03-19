using Test
using Testfiles
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
(isinteractive() ? includet : include)(joinpath(pkgdir(NetworkDynamics, "test", "ComponentLibrary.jl")))

BUILDKITE = haskey(ENV, "BUILDKITE")
BUILDKITE && @test CUDA.functional() # fail early in buildkite if cuda is not available

@testset "NetworkDynamics Tests" begin
    # skip majority of tests under BUILDKITE env
    if !BUILDKITE
        @testset "Package Quality Tests" begin
            # print_explicit_imports(NetworkDynamics)
            @test check_no_implicit_imports(NetworkDynamics) === nothing
            @test check_no_stale_explicit_imports(NetworkDynamics, ignore=(:Symbolics,)) === nothing
            Aqua.test_all(NetworkDynamics;
                ambiguities=false,
                stale_deps=VERSION ≥ v"1.11", # don't check stale deps on LTS (we add Testfiles to main env)
                deps_compat=VERSION ≥ v"1.11", # don't check compat on LTS
                persistent_tasks=false)
            @test_broken isempty(Docs.undocumented_names(NetworkDynamics))
        end

        @testfile "utils_test.jl"

        NetworkDynamics.CHECK_COMPONENT[] = false
        @testfile "construction_test.jl"
        @testfile "aggregators_test.jl"
        @testfile "massmatrix_test.jl"
        NetworkDynamics.CHECK_COMPONENT[] = true

        @testfile "initialization_test.jl"
        @testfile "callbacks_test.jl"
        @testfile "metadata_test.jl"
        @testfile "linear_analysis_test.jl"
        @testfile "show_test.jl"
        @testfile "spinners_test.jl"
        @testfile "sparsity_test.jl"

        @testfile "MTK_test.jl"

        # check on the precompile files
        @testfile "../src/precompile_workload.jl"
        @testfile "../ext/MTKExt_precomp_workload.jl"

        @testfile "AD_test.jl"
        @testfile "doctor_test.jl"
    end

    @testfile "symbolicindexing_test.jl"
    @testfile "external_inputs_test.jl"
    @testfile "loopback_test.jl"

    @testfile "diffusion_test.jl"
    @testfile "inhomogeneous_test.jl"

    if CUDA.functional()
        @testfile "GPU_test.jl"
    end

    # the docs should work with MTK loaded
    @testset "Test Doc Examples" begin
        examples = joinpath(pkgdir(NetworkDynamics), "docs", "examples")
        for file in readdir(examples; join=true)
            endswith(file, ".jl") || continue
            eval(:(@testfile $file))
        end
    end

    if !CUDA.functional()
        @test gethostname() != "hw-g14" # on this pc, CUDA *should* be available
        @warn "Skipped all CUDA tests because no device is available."
    end
end
