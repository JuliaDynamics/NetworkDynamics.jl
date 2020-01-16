using Pkg
using Documenter
push!(LOAD_PATH, "../")
using NetworkDynamics

makedocs(
    sitename = "NetworkDynamics",
    linkcheck = true,
    modules = [NetworkDynamics], # creates a warning if a ND.jl docstring is missing
    pages = [
    "General" => "index.md",
    "BasicConstructors.md",
    "Library.md",
    "Examples" => [
		"Examples.md",
        "StaticEdges.md",
		"DynamicEdges.md",
		"getting_started_with_network_dynamics.md"
	]
    ]) #Here we have to agree on the Page structure yet.
