ENV["GKSwstype"] = "100" # needed for plotting with GitHub Actions and GR (?)

using Documenter
using NetworkDynamics
using SciMLBase
using Literate
using ModelingToolkit

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

mtkext = Base.get_extension(NetworkDynamics, :MTKExt)
kwargs = (;
    root=joinpath(pkgdir(NetworkDynamics), "docs"),
    sitename="NetworkDynamics",
    modules=[NetworkDynamics, mtkext],
    linkcheck=true, # checks if external links resolve
    pagesonly=true,
    pages=[
        "General" => "index.md",
        "mathematical_model.md",
        "network_construction.md",
        "Features" => [
            "symbolic_indexing.md",
            "metadata.md",
            "initialization.md",
            "callbacks.md",
            "mtk_integration.md",
            "external_inputs.md",
        ],
        "API.md",
        "Tutorials" => [
            "Getting Started" => "generated/getting_started_with_network_dynamics.md",
            "Directed and Weighted Graphs" => "generated/directed_and_weighted_graphs.md",
            "Heterogeneous Systems" => "generated/heterogeneous_system.md",
            # "Stochastic differential equations" => "generated/StochasticSystem.md",
            # "Delay differential equations" => "generated/kuramoto_delay.md",
            "Cascading Failure" => "generated/cascading_failure.md",
            "Stress on Truss" => "generated/stress_on_truss.md",
            "Gas Network" => "generated/gas_network.md",
        ]
    ],
    draft=false,
    format = Documenter.HTML(ansicolor = true),
    # warnonly=true,
)

success = true
thrown_ex = nothing
try
    # strict build
    makedocs(; kwargs...)
catch e
    @info "Strict doc build failed, try again with warnonly=true"
    global success = false
    global thrown_ex = e
    # kwargs = (; kwargs..., warnonly=[:cross_references, :missing_docs, :docs_block])
    kwargs_warnonly = (; kwargs..., warnonly=true)
    makedocs(; kwargs_warnonly...)
end

deploydocs(; repo="github.com/JuliaDynamics/NetworkDynamics.jl.git",
           devbranch="main", push_preview=true)

if !success
    rethrow(thrown_ex)
end

# warnonly options
# :autodocs_block
# :cross_references
# :docs_block
# :doctest
# :eval_block
# :example_block
# :footnote
# :linkcheck_remotes
# :linkcheck
# :meta_block
# :missing_docs
# :parse_error
# :setup_block.
