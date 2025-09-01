#! julia --startup-file=no

#=
julia localmake.jl [--draft]

Runs the `make.jl` script and serves the docs on `localhost:8000`. 
You have to have `Revise` and `LiveServer` in your global environment for this script to work.

Two modes available:
- Full mode (default): Builds complete docs with examples, prompts to rebuild after each run
- Draft mode (--draft): Uses Documenter's draft mode with automatic rerendering, skips examples

The script will prompt whether to update the docs environment (10 second timeout, defaults to No).
At the end of each full build run, the user is prompted to rerun the make process. Using Revise
this will use updated `*.md` and source files, keeping the Julia session alive for faster builds.
=#

BUILD_DIR = joinpath(@__DIR__, "build")
mkpath(BUILD_DIR) # make sure path exists, otherwise the serve might fail

function readline_timeout(prompt, default, timeout)
    msg = Channel{String}(1)
    task = Task() do
        try eof(stdin); put!(msg, readline(stdin)); catch end
    end
    interrupter = Task() do
        sleep(timeout)
        istaskdone(task) || Base.throwto(task, InterruptException())
    end
    print(prompt)
    schedule(interrupter)
    schedule(task)
    wait(task)
    answer = if isempty(msg)
        println() # close line
        default
    else
        str = take!(msg)
        isempty(str) ? default : str
    end
    close(msg)
    return answer
end

function full_serving()
    @info "Start server..."
    port=8000
    servetask = @async serve(;dir=BUILD_DIR, port)
    errormonitor(servetask)

    run = true
    while run
        revise()
        @info "Start building docs..."
        try
            include("make.jl")
        catch e
            @error "make.jl error" exception=(e, catch_backtrace())
        end

        printstyled("\n\nDocs are served at http://localhost:$port\n\n", color=:blue, bold=true)
        println("Run again? Enter! Exit with 'q'.")
        if readline() == "q"
            run = false
        end
    end
end

function draft_serving()
    ENV["DOCUMENTER_DRAFT"] = "true"
    servedocs(
        foldername=".",
        literate=joinpath(@__DIR__, "examples"),
        skip_dir=joinpath(@__DIR__, "src", "generated")
    )
end

do_update = readline_timeout("Do you want to update docs environment? [y/N] ", "N", 10)

####
#### Set up environment
####
using Pkg
Pkg.activate(@__DIR__)
using Revise
using LiveServer

if VERSION ≤ v"1.11-"
    Pkg.develop(PackageSpec(path=dirname(@__DIR__))) # adds the package this script is called from
end

if do_update[1] == 'y'
    Pkg.update()
end
Pkg.instantiate()

draft_arg = findfirst(a -> a=="--draft", ARGS)
if !isnothing(draft_arg)
    @info "Found --draft, run in draft mode (no examples, automatic rerendering!)"
    draft_serving()
else
    @info "Start full build of docs including all examples"
    full_serving()
end
