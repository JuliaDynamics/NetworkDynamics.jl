# set root path of the ND repo
ndpath = begin
    path = splitpath(abspath(@__DIR__))
    i = findlast(isequal("NetworkDynamics"), path)
    abspath(path[begin:i]...)
end
bmpath = joinpath(ndpath, "benchmark")
cd(bmpath)

import Pkg; Pkg.activate(bmpath)
using PkgBenchmark

# setup the base config for the benchmarks
bmconfig(; id="HEAD") = BenchmarkConfig(id = id,
                                      env = Dict("JULIA_NUM_THREADS" => 4),
                                      juliacmd = `julia`)

# Example 1 : run benchmarks for current version
current = benchmarkpkg(ndpath, bmconfig())

# Example 2 : run benchmarks for tagged version
baseline = benchmarkpkg(ndpath, bmconfig("v0.5.0"))

judge(current, baseline)
