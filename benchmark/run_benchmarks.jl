#!julia

# make sure PkgBenchamark is available in your standard environment
using PkgBenchmark
using LibGit2
using Dates

baseline_id = isempty(ARGS) ? "main" : ARGS[1]
@info "Run benchmarks on current directory and compare to $baseline_id"

# set root path of the ND repo
ndpath = begin
    path = splitpath(abspath(@__DIR__))
    i = findlast(isequal("NetworkDynamics"), path)
    abspath(path[begin:i]...)
end
bmpath = joinpath(ndpath, "benchmark")

ndpath_tmp = tempname()
@info "Copy repo to $ndpath_tmp"
cp(ndpath, ndpath_tmp)
cd(ndpath_tmp)

# if directory is dirty, commet everything
isdirty = with(LibGit2.isdirty, GitRepo(ndpath_tmp))
if isdirty
    @info "Dirty directory, add everything to new commit!"
    @assert pwd() == ndpath_tmp
    run(`git add -A`)
    run(`git commit -m "tmp for benchmarking"`)
    # assert that repo is clean now
    with(LibGit2.isdirty, GitRepo(ndpath_tmp)) && throw(error("Repository is still dirty!"))
end

import Pkg; Pkg.activate(ndpath_tmp)

# setup the base config for the benchmarks
bmconfig(id=nothing) = BenchmarkConfig(;id,
                                       env = Dict("JULIA_NUM_THREADS" => 4),
                                       juliacmd = `julia`)


ts = Dates.format(now(), "yyyy-mm-dd_HHMMSS_")

# Run benchmarks for current version
@info "Run benchmarks for current directory"
target = benchmarkpkg(ndpath_tmp, bmconfig())
writeresults(joinpath(bmpath, ts*"target.data"), target)
export_markdown(joinpath(bmpath, ts*"target.md"), target)

# Run benchmarks defined in current directory for specific version
@info "Run benchmarks for $baseline_id"
script_dir = tempname()
cp(bmpath, script_dir)
script = joinpath(script_dir, "benchmarks.jl")
baseline = benchmarkpkg(ndpath_tmp, bmconfig(baseline_id); script)
writeresults(joinpath(bmpath, ts*"baseline.data"), baseline)
export_markdown(joinpath(bmpath, ts*"baseline.md"), baseline)

# compare the two
judgment = judge(target, baseline)
export_markdown(joinpath(bmpath, ts*"judgment.md"), judgment)
