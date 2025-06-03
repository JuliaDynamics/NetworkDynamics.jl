#! julia --startup-file=no

#=
julia localmake.jl PORT(optional)

Runs the `make.jl` script and serves the docs on `localhost:8000` (or any other port if specified).
You have to have `Revise` and `LiveServer` in your global environment for this script to work.

At the end of each run the user is prompted to rerun the make process. Using revise this will
use the updated `*.md` and source files. This way the Julia session keeps alive and the
individual builds are much faster.
=#

PORT = isempty(ARGS) ? 8000 : parse(Int, ARGS[1])
@assert 8000 ≤ PORT ≤ 9000 "PORT has to be in range 8000..9000!"

print("Do you want to update docs environment? [y/n] (default: n) ")
answer = readline()
UPDATE_ENV = if !isempty(answer) && answer[1] == 'y'
    true
else
    false
end

print("Do you want to build docs continously on file change? This will enable the `draft=true` flag, and the examples will not run. [y/n] (default: n)")
answer = readline()
SERVEDOCS_DRAFT = if !isempty(answer) && answer[1] == 'y'
    true
else
    false
end


using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

if UPDATE_ENV
    Pkg.update()
end

using Revise
using LiveServer
using NetworkDynamics
cd(pkgdir(NetworkDynamics))


if SERVEDOCS_DRAFT
    servedocs(
        literate_dir = joinpath("docs", "examples"),
        skip_dir = joinpath("docs", "src", "generated")
    )
else
    @info "Start server..."
    @async serve(;dir=joinpath(@__DIR__, "build"), PORT)

    while true
        revise()
        @info "Start building docs..."
        try
            include("make.jl")
        catch e
            @info "make.jl error" e
        end

        println("\nDocs are served at http://localhost:$port")

        println("Press [Enter] to rerun the make process or [Ctrl+C] to exit.")
        readline()
    end
end
