using NetworkDynamics
using ParallelTestRunner
using CUDA
using Test

# `runtests` runs every entry of `testsuite` in its own module on a worker process. The files are
# listed explicitly instead of autodiscovered, so this file stays the single place which decides
# what runs where (normal CI vs. the CUDA machine on buildkite).

TESTDIR = pkgdir(NetworkDynamics, "test") # not @__DIR__, that breaks when sent to a REPL
testsuite = Dict{String,Expr}()

# init_worker_code -> workers Main, init_code -> every test module
init_worker_code = quote
    # on a worker this is always `include`, in the REPL Revise tracks the helpers instead
    (isinteractive() ? includet : include)($(joinpath(TESTDIR, "testutils.jl")))
    (isinteractive() ? includet : include)($(joinpath(TESTDIR, "ComponentLibrary.jl")))
end

"""
    includetest!(file)

Add a test file:
- plain repl: `include(file)` in Main
- ] test/include("runtests.jl"): add to paarallel test pool and run on worker
"""
function includetest!(file)
    path = joinpath(TESTDIR, file)
    if any(f -> f.func === :include || f.func === :_include, stacktrace())
        testsuite[file] = :(Main.quiet_include(@__MODULE__, $path))
    else
        Core.eval(Main, init_worker_code) # no worker to set up here
        @testset "$file" begin include(path) end
    end
end

# on buildkite (CUDA machine) we only run the tests which actually exercise the GPU paths
BUILDKITE = haskey(ENV, "BUILDKITE")
BUILDKITE && !CUDA.functional() && error("CUDA is not functional on the buildkite agent!")

if !BUILDKITE
    includetest!("quality_test.jl");
    includetest!("utils_test.jl");
    includetest!("construction_test.jl");
    includetest!("aggregators_test.jl");
    includetest!("massmatrix_test.jl");
    includetest!("initialization_test.jl");
    includetest!("callbacks_test.jl");
    includetest!("metadata_test.jl");
    includetest!("alias_normalization_test.jl");
    includetest!("init_resolution_test.jl");
    includetest!("obsrules_test.jl");
    includetest!("linear_analysis_test.jl");
    includetest!("show_test.jl");
    includetest!("spinners_test.jl");
    includetest!("sparsity_test.jl");
    includetest!("MTK_test.jl");
    includetest!("bound_to_test.jl");
    includetest!("weak_test.jl");
    includetest!("optional_formula_test.jl");
    includetest!("default_from_test.jl");
    includetest!("AD_test.jl");
    includetest!("doctor_test.jl");

    # # check on the precompile files
    includetest!("../src/precompile_workload.jl");
    includetest!("../ext/MTKExt_precomp_workload.jl");
end

# stuff with GPU tests
includetest!("symbolicindexing_test.jl");
includetest!("external_inputs_test.jl");
includetest!("loopback_test.jl");
includetest!("diffusion_test.jl");
includetest!("inhomogeneous_test.jl");

# the docs should work with MTK loaded
for file in readdir(joinpath(TESTDIR, "..", "docs", "examples"))
    endswith(file, ".jl") || continue
    includetest!("../docs/examples/" * file);
end;

if CUDA.functional()
    includetest!("GPU_test.jl")
else
    gethostname() == "hw-g14" && error("CUDA should be available on this machine!")
    @warn "Skipped all CUDA tests because no device is available."
end

# empty when the lines above were sent to a REPL one by one, those tests already ran
isempty(testsuite) || runtests(NetworkDynamics, ARGS; testsuite, init_worker_code)
