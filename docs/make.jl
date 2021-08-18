ENV["GKSwstype"] = "100" # needed for plotting with GitHub Actions and GR (?)

using Documenter
using NetworkDynamics

makedocs(
    sitename = "NetworkDynamics",
    linkcheck=true, # not sure if we need that
    modules = [NetworkDynamics], # creates a warning if a ND.jl docstring is missing
    pages = [
    "General" => "index.md",
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
		"Delay differential equations" => "getting_started_with_DDEs.md",
    	]
    ])

deploydocs(
    repo = "github.com/PIK-ICoNe/NetworkDynamics.jl.git",
    devbranch = "main"
)
