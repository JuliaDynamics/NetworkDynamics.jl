#!julia

baseline_id = isempty(ARGS) ? "main" : ARGS[1]
@info "Run benchmarks on current directory and compare to $baseline_id"

# set root path of the ND repo
ndpath = begin
    path = splitpath(abspath(@__DIR__))
    i = findlast(isequal("NetworkDynamics"), path)
    abspath(path[begin:i]...)
end
bmpath = joinpath(ndpath, "benchmark")
cd(ndpath)
import Pkg; Pkg.activate(ndpath)

# make sure PkgBenchamark is available in your standard environment
using PkgBenchmark
using PkgBenchmark.LibGit2
isdirty = LibGit2.with(LibGit2.isdirty, LibGit2.GitRepo(ndpath))
isdirty && throw(error("Git environment $ndpath is dirty! Please stash or commit changes."))

# setup the base config for the benchmarks
bmconfig(id=nothing) = BenchmarkConfig(;id,
                                       env = Dict("JULIA_NUM_THREADS" => 4),
                                       juliacmd = `julia`)

# Run benchmarks for current version
target = benchmarkpkg(ndpath, bmconfig())
writeresults(joinpath(bmpath, "target.data"), target)
export_markdown(joinpath(bmpath, "target.md"), target)

# Run benchmarks defined in current directory for specific version
tmpdir = tempname()
cp(bmpath, tmpdir)
script = joinpath(tmpdir, "benchmarks.jl")
baseline = benchmarkpkg(ndpath, bmconfig(baseline_id); script)
writeresults(joinpath(bmpath, "baseline.data"), baseline)
export_markdown(joinpath(bmpath, "baseline.md"), baseline)

# compare the two
judgment = judge(target, baseline)
export_markdown(joinpath(bmpath, "judgment.md"), judgment)
