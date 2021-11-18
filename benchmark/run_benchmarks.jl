#!/usr/bin/env julia

# record start time
tstart = time()

# set root path of the ND repo
NDPATH = begin
    path = splitpath(abspath(@__DIR__))
    i = findlast(contains("NetworkDynamics"), path)
    abspath(path[begin:i]...)
end
BMPATH = joinpath(NDPATH, "benchmark")

# activate the benchmark folder as current environment
io = IOBuffer() # hide output
import Pkg;
Pkg.activate(BMPATH; io);
Pkg.develop(; path=NDPATH, io);

using PkgBenchmark
using LibGit2
using Dates
using Random
using ArgParse

s = ArgParseSettings()
#! format: off
@add_arg_table! s begin
    "--target", "-t"
        help = "Specify branch, commit id, tag, … for target benchmark. If `directory` use the current state of the directory. If *.date file use prev. exporte raw results."
        default = "directory"
    "--baseline", "-b"
        help = "Same for baseline benchmark. If `none` or `nothing` skip the baseline benckmark and coparison step."
        default = "main"
    "--command"
        help = "Julia command for benchmarks"
        default = "julia"
    "--tcommand"
        help = "Julia command for target benchmarks. If nothing use `--command`"
    "--bcommand"
        help = "Julia command for baseline benchmarks. If nothing use `--command`"
    "--threads"
        help = "Set number of threads."
        arg_type = Int
        default = 4
    "--prefix", "-p"
        help = "Prefix for filenames, defaults to timestamp"
        default = Dates.format(now(), "yyyy-mm-dd_HHMMSS")
    "--verbose", "-v"
        help = "Print out the current benchmark."
        action = :store_true
    "--retune"
        help = "Force retuneing of parameters befor the first benchmark. (Second benchmark will use the tune file)"
        action = :store_true
    "--exportraw"
        help = "Export raw data of trials. I.e. to use benchmarks results again als baseline or target."
        action = :store_true
end
#! format: on
args = parse_args(s; as_symbols=true)
args[:prefix] *= '_'

@info "Run benchmarks on $(args[:target]) and compare to $(args[:baseline])"

#=
PkgBenchmark will do stuff to the repo (i.e. check out the desired branch),
therefore we copy the repo to tmp/NetworkDynamics
We want to have the same benchmarks for target and base therefore we copy
NetworkDynamics/benchmark to tmp/benchmark so it won't change if the repo
state changes.
=#
tmp = mkpath(tempname());
ndpath_tmp = joinpath(tmp, "NetworkDynamics")
@info "Copy repo to $ndpath_tmp"
cp(NDPATH, ndpath_tmp)
cd(ndpath_tmp) # the tmp repo is a good place to chill

# copy bm scripts to separate folder (so it won't change if the state of the nd_temp repo chagnes)
bmpath_tmp = joinpath(tmp, "benchmark")
cp(BMPATH, bmpath_tmp)
script_path = joinpath(bmpath_tmp, "benchmarks.jl")

# if directory is dirty, commit everything (otherwise we cant switch do a differen commit/branch)
isdirty = with(LibGit2.isdirty, GitRepo(ndpath_tmp))
if isdirty
    @info "Dirty directory, add everything to new commit!"
    @assert realpath(pwd()) == realpath(ndpath_tmp) "Julia is in $(pwd()) not it $ndpath_tmp"
    run(`git checkout -b $(randstring(15))`)
    run(`git add -A`)
    run(`git commit -m "tmp commit for benchmarking"`)
    # assert that repo is clean now
    with(LibGit2.isdirty, GitRepo(ndpath_tmp)) && throw(error("Repository is still dirty!"))
end

# activate the environment in the tmp copy of the scripts
# env will be used for the julia processes which do the benchmarking
Pkg.activate(bmpath_tmp)
Pkg.develop(; path=ndpath_tmp)

"""
Create a runnable `Cmd` object from a string like "juila-1.5".
Due to a bug in PkgBenchmark this also adds the default sysimge
to the command. See

https://github.com/JuliaCI/PkgBenchmark.jl/issues/136
"""
function string_to_command(cmd::AbstractString)
    @assert !contains(cmd, " ") "julia commands with arguments not allowed, got $cmd"
    out = IOBuffer()
    capturecommand = `$cmd --startup-file=no -e "println(Base.julia_cmd()[3])"`
    run(pipeline(capturecommand; stdout=out))
    sysimg = String(take!(out))[1:end-1] # remove \n
    @assert contains(sysimg, r"^-J.*[dylib|so]$") "Captured sysimg looks weird? $sysimg"
    return Cmd([cmd, sysimg])
end

"""
    benchmark(; name, id, cmd)

  - name determines the filenames
  - id determines the git id, i.e. branch name. If id=="directory" use currend
    directory even if dirty!
  - cmd if empty use default julia command
"""
function benchmark(; name, id, cmd, retune)
    # magic id string directory will benchmark dirty state
    id = id == "directory" ? nothing : id

    # PkgBenchmarks does not call Pkg.resolve after checking out a different
    # state of ndpath_tmp! Therefore we need to check out the desired commit in
    # the tmp directory and resolve the manifest of the currently activated
    # environment (-> bmpath_tmp)
    if !isnothing(id)
        @assert realpath(pwd()) == realpath(ndpath_tmp) "Julia is in $(pwd()) not it $ndpath_tmp"
        run(`git checkout $id`)
        Pkg.resolve()
    end

    # choose specific command only if given, elso the default command
    cmd = isnothing(cmd) ? args[:command] : cmd
    @info "Run $name benchmarks on $id with $cmd"

    config = BenchmarkConfig(; id,
                             env=Dict("JULIA_NUM_THREADS" => args[:threads]),
                             juliacmd=string_to_command(cmd))
    results = benchmarkpkg(ndpath_tmp, config;
                           script=script_path,
                           retune,
                           verbose=args[:verbose])

    args[:exportraw] && writeresults(joinpath(BMPATH, args[:prefix] * name * ".data"), results)
    export_markdown(joinpath(BMPATH, args[:prefix] * name * ".md"), results)
    return results
end

if contains(args[:target], r".data$")
    path = joinpath(BMPATH, args[:target])
    target = readresults(joinpath(BMPATH, args[:target]))
else
    target = benchmark(; name="target", id=args[:target], cmd=args[:tcommand], retune=args[:retune])
end

if args[:baseline] ∉ ["nothing", "none"]
    if contains(args[:baseline], r".data$")
        path = joinpath(BMPATH, args[:baseline])
        baseline = readresults(joinpath(BMPATH, args[:baseline]))
    elseif args[:baseline] ∉ ["nothing", "none"]
        baseline = benchmark(; name="baseline", id=args[:baseline], cmd=args[:bcommand], retune=false)
    end
    # compare the two
    judgment = judge(target, baseline)
    export_markdown(joinpath(BMPATH, args[:prefix] * "judgment.md"), judgment)
end

# copy tune file over to real repo if there is none yet or it was retuned
if !isfile(joinpath(BMPATH, "tune.json")) || args[:retune]
    println("Update tune.json in $BMPATH")
    cp(joinpath(ndpath_tmp, "benchmark", "tune.json"), joinpath(BMPATH, "tune.json"); force=args[:retune])
end

s = round(Int, time() - tstart)
m, s = s ÷ 60, s % 60
@info "Benchmarking endet after $m min $s seconds"
