#!julia

# make sure PkgBenchamark is available in your standard environment
using PkgBenchmark
using LibGit2
using Dates

tstart = time()

function getarg(default, flags...)
    idx = findfirst(x -> x ∈ flags, ARGS)
    return idx===nothing ? default : ARGS[idx+1]
end

# target defaults to nothing which means the current directory
target_id   = getarg("directory", "-t", "--target")
# baseline defaults to "main" branach. if "nothing" or "none" don't run second benchmark
baseline_id = getarg("main", "-b", "--baseline")
# number of threads, defaults to 4
num_threads = parse(Int, getarg("4", "--threads"))
# chose julia executable (the default one)
julia_cmd   = getarg("julia", "--command")
# chose julia executable for target specificially
julia_cmd_target = getarg("", "--tcommand")
# chose julia executable for baseline specificially
julia_cmd_baseline = getarg("", "--bcommand")
# verbose flag, print current benchmark
verbose = "-v" ∈ ARGS || "--verbose" ∈ ARGS
# force retune
retune = "--retune" ∈ ARGS
# export raw benchmark data
exportraw = "--exportraw" ∈ ARGS
# target defaults to nothing which means the current directory
prefix = getarg(Dates.format(now(), "yyyy-mm-dd_HHMMSS"), "-p", "--prefix") * "_"

@info "Run benchmarks on $target_id and compare to $baseline_id"

# set root path of the ND repo
ndpath = begin
    path = splitpath(abspath(@__DIR__))
    i = findlast(isequal("NetworkDynamics"), path)
    abspath(path[begin:i]...)
end
bmpath = joinpath(ndpath, "benchmark")

ndpath_tmp = tempname()*"_NetworkDynamics"
@info "Copy repo to $ndpath_tmp"
cp(ndpath, ndpath_tmp)
cd(ndpath_tmp)

# if directory is dirty, commet everything
isdirty = with(LibGit2.isdirty, GitRepo(ndpath_tmp))
if isdirty
    @info "Dirty directory, add everything to new commit!"
    @assert realpath(pwd()) == realpath(ndpath_tmp) "Julia is in $(pwd()) not it $ndpath_tmp"
    run(`git add -A`)
    run(`git commit -m "tmp for benchmarking"`)
    # assert that repo is clean now
    with(LibGit2.isdirty, GitRepo(ndpath_tmp)) && throw(error("Repository is still dirty!"))
end

import Pkg; Pkg.activate(ndpath_tmp)

# copy bmscripts to separate folder
script_dir = tempname()*"_benchmark"
cp(bmpath, script_dir)
script_path = joinpath(script_dir, "benchmarks.jl")

"""
Create a runnable `Cmd` object from a string like "juila-1.5".
Due to a bug in PkgBenchmark this also adds the default sysimge
to the command. See

https://github.com/JuliaCI/PkgBenchmark.jl/issues/136
"""
function string_to_command(cmd::AbstractString)
    @assert !contains(cmd, " ") "julia commands with arguments not allowed, got $cmd"
    out = IOBuffer()
    capturecommand = `$cmd -e "println(Base.julia_cmd()[3])"`
    run(pipeline(capturecommand; stdout=out))
    sysimg = String(take!(out))[1:end-1] # remove \n
    @assert contains(sysimg, r"^-J.*dylib$") "Captured sysimg looks weird? $sysimg"
    return Cmd([cmd, sysimg])
end

"""
    benchmark(; name, id, cmd)

- name determines the filenames
- id determines the git id, i.e. branch name. If id=="directory" use currend
directory even if dirty!
- cmd if empty use default julia command
"""
function benchmark(; name, id, cmd::AbstractString, retune)
    # choose specific command only if given
    cmd = isempty(cmd) ? julia_cmd : cmd
    @info "Run $name benchmarks on $id with $cmd"
    # magic id directory will benchmark dirty state
    id = id=="directory" ? nothing : id
    config = BenchmarkConfig(;id,
                             env = Dict("JULIA_NUM_THREADS" => num_threads),
                             juliacmd = string_to_command(cmd))
    results = benchmarkpkg(ndpath_tmp, config;
                           script=script_path,
                           retune,
                           verbose)
    exportraw && writeresults(joinpath(bmpath, prefix*name*".data"), results)
    export_markdown(joinpath(bmpath, prefix*name*".md"), results)
    return results
end

target = benchmark(; name="target", id=target_id, cmd=julia_cmd_target, retune)

if baseline_id ∉ ["nothing", "none"]
    baseline = benchmark(; name="baseline", id=baseline_id, cmd=julia_cmd_baseline, retune=false)

    # compare the two
    judgment = judge(target, baseline)
    export_markdown(joinpath(bmpath, prefix*"judgment.md"), judgment)
end

# copy tune file over to real repo if there is none yet or it was retuned
if !isfile(joinpath(bmpath, "tune.json")) || retune
    println("Update tune.json in $bmpath")
    cp(joinpath(ndpath_tmp, "benchmark", "tune.json"), joinpath(bmpath, "tune.json"); force=retune)
end

s = round(Int, time()-tstart)
m, s = s ÷ 60, s % 60
@info "Benchmarking endet after $m min $s seconds"
