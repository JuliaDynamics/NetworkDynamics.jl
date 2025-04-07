# ENV["GKSwstype"] = "100" # needed for plotting with GitHub Actions and GR (?)

using Documenter
using NetworkDynamics
using NetworkDynamicsInspector
using NetworkDynamicsInspector: NetworkDynamicsInspector as NDI
using SciMLBase
using Literate
using ModelingToolkit
using DocumenterInterLinks
using Electron


@info "Create global electron window"
NDI.CURRENT_DISPLAY[] = ElectronDisp(resolution=(1200,800)) # hide
NDI.get_electron_window()

links = InterLinks(
    "DiffEq" => "https://docs.sciml.ai/DiffEqDocs/stable/",
    "MTK" => "https://docs.sciml.ai/ModelingToolkit/stable/",
    "SymbolicIndexingInterface" => "https://docs.sciml.ai/SymbolicIndexingInterface/stable/",
    "DiffEqCallbacks" => "https://docs.sciml.ai/DiffEqCallbacks/stable/",
)

# generate examples
example_dir = joinpath(@__DIR__, "examples")
outdir = joinpath(@__DIR__, "src", "generated")
isdir(outdir) && rm(outdir, recursive=true)
mkpath(outdir)

for example in filter(contains(r".jl$"), readdir(example_dir, join=true))
    Literate.markdown(example, outdir)
    Literate.script(example, outdir; keep_comments=true)
end

mtkext = Base.get_extension(NetworkDynamics, :NetworkDynamicsMTKExt)
dfext = Base.get_extension(NetworkDynamics, :NetworkDynamicsDataFramesExt)
kwargs = (;
    root=joinpath(pkgdir(NetworkDynamics), "docs"),
    sitename="NetworkDynamics",
    modules=[NetworkDynamics, mtkext, dfext, NetworkDynamicsInspector],
    linkcheck=true, # checks if external links resolve
    pagesonly=true,
    plugins=[links],
    pages=[
        "General" => "index.md",
        "mathematical_model.md",
        "network_construction.md",
        "data_structure.md",
        "Features" => [
            "symbolic_indexing.md",
            "metadata.md",
            "initialization.md",
            "callbacks.md",
            "mtk_integration.md",
            "external_inputs.md",
            "inspector.md",
        ],
        "API.md",
        "Tutorials" => [
            "Getting Started" => "generated/getting_started_with_network_dynamics.md",
            "Heterogeneous Systems" => "generated/heterogeneous_system.md",
            "Initialization" => "generated/init_tutorial.md",
            "Cascading Failure" => "generated/cascading_failure.md",
            "Gas Network" => "generated/gas_network.md",
            "Stress on Truss" => "generated/stress_on_truss.md",
            "Directed and Weighted Graphs" => "generated/directed_and_weighted_graphs.md",
        ]
    ],
    draft=false,
    format = Documenter.HTML(ansicolor = true),
    warnonly=[:missing_docs],
)
kwargs_warnonly = (; kwargs..., warnonly=true)

if haskey(ENV,"GITHUB_ACTIONS")
    success = true
    thrown_ex = nothing
    try
        makedocs(; kwargs...)
    catch e
        @info "Strict doc build failed, try again with warnonly=true"
        global success = false
        global thrown_ex = e
        makedocs(; kwargs_warnonly...)
    end

    deploydocs(; repo="github.com/JuliaDynamics/NetworkDynamics.jl.git",
            devbranch="main", push_preview=true)

    success || throw(thrown_ex)
else # local build
    makedocs(; kwargs_warnonly...)
end

@info "Close global electron window"
NDI.close_display()


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
