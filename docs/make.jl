using Documenter
using NetworkDynamics

makedocs(
    sitename = "NetworkDynamics",
    format = :html,
    modules = [NetworkDynamics]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
