using Pkg
using Documenter
push!(LOAD_PATH, "../")
using NetworkDynamics

makedocs(
    sitename = "NetworkDynamics",
    linkcheck=true, # not sure why we need that
    modules = [NetworkDynamics], # creates a warning if a ND.jl docstring is missing
    pages = [
    "General" => "index.md",
    "BasicConstructors.md",
    "Library.md",
    "Examples" => [
		# "Examples.md",
        # "StaticEdges.md",
		# "DynamicEdges.md",
		"getting_started_with_network_dynamics.md"
    ]) #Here we have to agree on the Page structure yet.

# Documenter automatically deploys documentation to gh-pages.
# However at the moment there is a bug when GITHUB_TOKEN is used as
# authentication method by the GithubActions worker: The html pages build is
# not triggered automatically and while the docs are updated in gh-pages branch
# they are not visible on github.io. A workaround is to manually trigger the
# html build by commiting anything to the gh-pages branch. Alternatively the
# SSH_DEPLOY_KEY can be used as authentication method for the worker.
deploydocs(
    repo = "github.com/FHell/NetworkDynamics.jl.git",
#   no longernecessary
#   deploy_config = Documenter.GitHubActions(), # this should work, but it's strange
)
