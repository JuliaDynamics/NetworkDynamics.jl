#! julia --startup-file=no

#=
julia localmake.jl PORT(optional)

Runs the `make.jl` script and serves the docs on `localhost:8000` (or any other port if specified).
You have to have `Revise` and `LiveServer` in your global environment for this script to work.

At the end of each run the user is prompted to rerun the make process. Using revise this will
use the updated `*.md` and source files. This way the Julia session keeps alive and the
individual builds are much faster.
=#

using Pkg
Pkg.activate(@__DIR__)
Pkg.develop(PackageSpec(path=dirname(@__DIR__))) # adds the package this script is called from
Pkg.instantiate()
Pkg.update()

using Revise
using LiveServer
using REPL.TerminalMenus

port = isempty(ARGS) ? 8000 : parse(Int, ARGS[1])
@assert 8000 ≤ port ≤ 9000 "port has to be in range 8000..9000!"

@info "Start server..."
@async serve(;dir=joinpath(@__DIR__, "build"), port)

menu = RadioMenu(["Run again!", "Quit!"])
while true
    revise()
    @info "Start building docs..."
    try
        include("make.jl")
    catch e
        @info "make.jl error" e
    end

    println("\nDocs are served at http://localhost:$port")

    if request("What now?", menu) != 1
        break
    end
end
