ENV["GKSwstype"] = "100" # needed for plotting with GitHub Actions and GR (?)

using Documenter
using NetworkDynamics
using Literate

# generate examples
example_dir = joinpath(@__DIR__, "..", "examples")
outdir = joinpath(@__DIR__, "src", "generated")
isdir(outdir) && rm(outdir, recursive=true)
mkpath(outdir)

for example in filter(contains(r".jl$"), readdir(example_dir, join=true))
    Literate.markdown(example, outdir)
    Literate.script(example, outdir)
end

makedocs(; sitename="NetworkDynamics",
         linkcheck=true, # checks if external links resolve
         pages=["General" => "index.md",
             "BasicConstructors.md",
             "parameters.md",
             "Multithreading.md",
             "Library.md",
             "accessing_edge_variables.md",
             "Tutorials" => [
                 "Getting started" => "getting_started_with_network_dynamics.md",
                 "Directed and weighted graphs" => "directed_and_weighted_graphs.md",
                 "Heterogeneous systems" => "heterogeneous_system.md",
                 "Stochastic differential equations" => "SDEVertex.md",
                 "Delay differential equations" => "kuramoto_delay.md"]])

deploydocs(; repo="github.com/PIK-ICoNe/NetworkDynamics.jl.git",
           devbranch="main")
