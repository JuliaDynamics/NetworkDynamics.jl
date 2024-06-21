ENV["GKSwstype"] = "100" # needed for plotting with GitHub Actions and GR (?)

using Documenter
using NetworkDynamics
using SciMLBase
using Literate

# generate examples
example_dir = joinpath(@__DIR__, "examples")
outdir = joinpath(@__DIR__, "src", "generated")
isdir(outdir) && rm(outdir, recursive=true)
mkpath(outdir)

for example in filter(contains(r".jl$"), readdir(example_dir, join=true))
    Literate.markdown(example, outdir)
    Literate.script(example, outdir; keep_comments=true)
end

# TODO: doc on steady state solve https://docs.sciml.ai/NonlinearSolve/stable/native/steadystatediffeq/#SteadyStateDiffEq.SSRootfind
# TODO: doc on parameter & state handling? -> symbolic indexing

makedocs(; root=joinpath(pkgdir(NetworkDynamics), "docs"),
         sitename="NetworkDynamics",
         modules=[NetworkDynamics],
         linkcheck=true, # checks if external links resolve
         pagesonly=true,
         pages=["General" => "index.md",
             # "BasicConstructors.md",
             # "parameters.md",
             # "Multithreading.md",
             "Library.md",
             # "accessing_edge_variables.md",
             "Tutorials" => [
                 "Getting started" => "generated/getting_started_with_network_dynamics.md",
                 "Directed and weighted graphs" => "generated/directed_and_weighted_graphs.md",
                 "Heterogeneous systems" => "generated/heterogeneous_system.md",
                 "Stochastic differential equations" => "generated/StochasticSystem.md",
                 # "Delay differential equations" => "generated/kuramoto_delay.md",
                 "Cascading failure" => "generated/cascading_failure.md",]
         ],
         draft=false,
         format = Documenter.HTML(ansicolor = true),
         warnonly=true)

deploydocs(; repo="github.com/PIK-ICoNe/NetworkDynamics.jl.git",
           devbranch="main", push_preview=true)
