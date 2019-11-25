using Pkg
# For GitHubActions support, we need to get the latest version (DEV) of Documenter
Pkg.add(PackageSpec(url="https://github.com/JuliaDocs/Documenter.jl"))
using Documenter
push!(LOAD_PATH, "../")
using NetworkDynamics

makedocs(
    sitename = "NetworkDynamics",
    format = Documenter.HTML(canonical = "https://FHell.github.io/NetworkDynamics.jl/stable"),
    linkcheck=true,
    modules = [NetworkDynamics],
    pages = [
    "General" => "index.md",
    "Functions_and_Constructors.md",
    "StaticEdges.md",
    "DynamicEdges.md",
    ]) #Here we have to agree on the Page structure yet.

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
# We are experimenting with the GitHubActions CI, which is not yet supported 
# by a stable branch of Documenter
deploydocs(
    repo = "github.com/FHell/NetworkDynamics.jl.git",
    deploy_config = Documenter.GitHubActions(), # this should work, but it's strange
    target = "build"
)
