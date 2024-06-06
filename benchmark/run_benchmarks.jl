#!julia --startup-file=no
@info "Start Benchmarking..."

original_path = pwd()

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
using Serialization
using StyledStrings
using CairoMakie

(isinteractive() ? includet : include)("benchmark_utils.jl")

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
    "--prefix", "-p"
        help = "Prefix for filenames, defaults to timestamp"
        default = Dates.format(now(), "yyyy-mm-dd_HHMMSS")
    "--verbose", "-v"
        help = "Print out the current benchmark."
        action = :store_true
    "--no-data-export"
        help = "Don't copy the .data files to working directory."
        action = :store_true
    "--no-plot"
        help = "Don't copy plot to the working directory."
        action = :store_true
    "--export-md"
        help = "Export Comparison table as markdown file.q"
        action = :store_true
end
#! format: on
args = parse_args(s; as_symbols=true)
if args[:prefix] != ""
    args[:prefix] *= '_'
end

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
cd(ndpath_tmp)

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

# # activate the environment in the tmp copy of the scripts
# # env will be used for the julia processes which do the benchmarking
# Pkg.activate(bmpath_tmp)
# Pkg.develop(; path=ndpath_tmp)


"""
    benchmark(; name, rev, cmd)

  - name determines the filenames
  - rev determines the git rev, i.e. branch name. If rev=="directory" use currend
    directory even if dirty!
  - cmd if empty use default julia command
"""
function benchmark(; name, rev, cmd)
    # magic rev string directory will benchmark dirty state
    if rev != "directory"
        @assert realpath(pwd()) == realpath(ndpath_tmp) "Julia is in $(pwd()) not it $ndpath_tmp"
        run(`git checkout $rev`)
    end

    # choose specific command only if given, elso the default command
    cmd_str = isnothing(cmd) ? args[:command] : cmd
    @assert !contains(cmd_str, r"--startup-file") "--startup-file is not implemeted"
    @assert !contains(cmd_str, r"--project") "--project is not implemeted"

    @info "Prepare environment"
    prepare_cmd = "import Pkg; Pkg.develop(path=\"$ndpath_tmp\"); Pkg.update()"
    run(`$cmd_str --startup-file=no --project=$bmpath_tmp -e $prepare_cmd`)

    exp_tmp = tempname(bmpath_tmp)
    cmd = `$cmd_str --startup-file=no --project=$bmpath_tmp $script_path $exp_tmp`
    @info "Run $name benchmarks on \"$rev\" with $cmd"
    println()
    println("------------ spwan benchmark process --------------")
    println()

    run(cmd)

    println()
    println("------------ benchmark process ended --------------")
    println()

    args[Symbol("no-data-export")] || cp(exp_tmp, joinpath(BMPATH, args[:prefix] * name * ".data"))
    result = deserialize(exp_tmp)
end

target = if contains(args[:target], r".data$")
    path = joinpath(original_path, args[:target])
    deserialize(path)
elseif args[:target] == "latest"
    file = sort(filter(contains("target.data"), readdir(original_path)))[end]
    @info "Use file $file as target"
    deserialize(joinpath(original_path, file))
else
    benchmark(; name="target", rev=args[:target], cmd=args[:tcommand])
end

baseline = if args[:baseline] ∉ ["nothing", "none"]
    if args[:baseline] == "latest"
        file = sort(filter(contains("baseline.data"), readdir(original_path)))[end]
        @info "Use file $file as baseline"
        deserialize(joinpath(original_path, file))
    elseif contains(args[:baseline], r".data$")
        path = joinpath(original_path, args[:baseline])
        deserialize(path)
    elseif args[:baseline] ∉ ["nothing", "none"]
        benchmark(; name="baseline", rev=args[:baseline], cmd=args[:bcommand])
    end
else
    nothing
end

if !isnothing(baseline)
    println(styled"{bright_red:Baseline}")
    println()
    display(baseline)
end

println()
println(styled"{bright_red:Target}")
println()
display(target)
println()

if !isnothing(baseline)
    println()
    println(styled"{bright_red:Comparison}")
    println()
    comp = compare(target, baseline)
    display(comp)

    res = test_return_values(comp)
    failed = res.anynonpass

    if !args[Symbol("no-plot")]
        @info "Save plot..."
        fig = plot_over_N(target, baseline)
        save(joinpath(original_path, args[:prefix] * "comparison.pdf"), fig)

    end
    if args[Symbol("export-md")]
        path = joinpath(original_path, args[:prefix] * "comparison.md")
        @info "Save markdown table to $path"
        open(path, "w") do io
            pretty_table(io, comp; backend=Val(:markdown))
        end
    end
else
    failed=false
end

s = round(Int, time() - tstart)
m, s = s ÷ 60, s % 60
@info "Benchmarking endet after $m min $s seconds"

exit(failed)
