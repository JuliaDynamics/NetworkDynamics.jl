# set root path of the ND repo
ndpath = begin
    path = splitpath(abspath(@__DIR__))
    i = findlast(isequal("NetworkDynamics"), path)
    abspath(path[begin:i]...)
end
bmpath = joinpath(ndpath, "benchmark")
cd(bmpath)

# import Pkg; Pkg.activate(bmpath)
using PkgBenchmark

# setup the base config for the benchmarks
bmconfig(id=nothing) = BenchmarkConfig(;id,
                                       env = Dict("JULIA_NUM_THREADS" => 4),
                                       juliacmd = `julia`)

# Example 1 : run benchmarks for current version
current = benchmarkpkg(ndpath, bmconfig())

# Example 2 : run benchmarks for tagged version
baseline = benchmarkpkg(ndpath, bmconfig("v0.5.0"))

# Example 3: run benchmakrs from HEAD for specific version
tmpdir = tempname()
cp(bmpath, tmpdir)
script = joinpath(tmpdir, "benchmarks.jl")
baseline = benchmarkpkg(ndpath, bmconfig("v0.5.0"); script, retune=true)

judge(current, baseline)

judge(ndpath, bmconfig(), bmconfig("v0.5.0"))
