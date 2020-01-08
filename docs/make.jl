using Pkg
using Documenter
push!(LOAD_PATH, "../")
using NetworkDynamics

makedocs(
    sitename = "NetworkDynamics",
    # linkcheck=true, # not sure why we would need that
    modules = [NetworkDynamics], # creates a warning if a ND.jl docstring is missing
    pages = [
    "General" => "index.md",
    "Functions_and_Constructors.md",
    "StaticEdges.md",
    "DynamicEdges.md",
    ]) #Here we have to agree on the Page structure yet.

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
# We are experimenting with the GitHubActions CI,
deploydocs(
    repo = "github.com/FHell/NetworkDynamics.jl.git",
#   probably not necessary
#    deploy_config = Documenter.GitHubActions(), # this should work, but it's strange
)
