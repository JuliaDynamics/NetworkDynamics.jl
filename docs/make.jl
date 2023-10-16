ENV["GKSwstype"] = "100" # needed for plotting with GitHub Actions and GR (?)

using Documenter
using NetworkDynamics

makedocs(; sitename  = "NetworkDynamics",
         linkcheck = true, # checks if external links resolve
         pages     = ["General" => "index.md",
         "BasicConstructors.md",
         "parameters.md",
         "Multithreading.md",
         "Library.md",
         "accessing_edge_variables.md",
         "Tutorials" => ["Getting started" => "getting_started_with_network_dynamics.md",
         "Directed and weighted graphs" => "directed_and_weighted_graphs.md",
         "Heterogeneous systems" => "heterogeneous_system.md",
         "Stochastic differential equations" => "SDEVertex.md",
         "Delay differential equations" => "kuramoto_delay.md"]])

deploydocs(; repo="github.com/PIK-ICoNe/NetworkDynamics.jl.git",
           devbranch="main")
