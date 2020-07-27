#=
julia --color=yes localmake.jl
python3 -m http.server --bind localhost
=#


using Pkg
Pkg.activate(@__DIR__)
push!(LOAD_PATH, Base.Filesystem.abspath("../"))
using Documenter
using NetworkDynamics

makedocs(
    sitename = "NetworkDynamics",
    linkcheck = true,
    modules = [NetworkDynamics], # creates a warning if a ND.jl docstring is missing
    pages = [
    "General" => "index.md",
    "BasicConstructors.md",
	"parameters.md",
	"Multithreading.md",
    "Library.md",
    "Tutorials" => [
		"Getting started" => "getting_started_with_network_dynamics.md",
		"Directed and weighted graphs" => "directed_and_weighted_graphs.md",
		"Delay differential equations" => "DDEVertex.md",
<<<<<<< HEAD
		"Heterogeneous systems" => "heterogeneous_system.md",
=======
		"Hetereogeneous systems" => "heterogeneous_system.md",
		"Stochastic differential equations" => "SDEVertex.md",
>>>>>>> 249d9dc092ff267118995bde6dcdb85e790a6709
		]
   	])
