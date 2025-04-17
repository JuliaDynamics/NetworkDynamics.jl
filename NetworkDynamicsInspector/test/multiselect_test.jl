using Bonito
using Electron
using NetworkDynamicsInspector
using NetworkDynamicsInspector: OptionGroup, TomSelect
using NetworkDynamicsInspector: NetworkDynamicsInspector as NDI

function electron(app)
    try
        close.(Electron.applications())
    catch
    end
    disp = Bonito.use_electron_display(devtools = true)
    display(disp, app)
    disp
end

@testset "Test TomSelect" begin
    ENV["JULIA_DEBUG"] = NetworkDynamicsInspector
    ENV["JULIA_DEBUG"] = ""
    gui = (;
        options = Observable{Vector{Union{Symbol,OptionGroup{Symbol}}}}([
            OptionGroup("Programming Languages", [:Julia, :Rust, :Java]),
            OptionGroup("Languages", [:French, :Spanish, :German]),
            :foo,
            :bar
        ]),
        selection = Observable{Vector{Symbol}}(Symbol[:Rust, :French]),
        selection_single = Observable{Vector{Symbol}}(Symbol[:Julia])
    )
    disp = App(;) do session
        NetworkDynamicsInspector.clear_obs_and_close!(gui)

        ms1 = TomSelect(gui.options, gui.selection; placeholder="multi", T=Symbol)
        ms2 = TomSelect(gui.options, gui.selection_single; placeholder="single", multi=false, T=Symbol)

        return Bonito.Grid(ms1, ms2; columns="100%", width="500px")
    end |> electron

    @test gui.selection[] == [:Rust, :French]
    gui.selection[] = [:Rust, :French, :Julia, :foo]
    gui.selection[] = []

    @test gui.selection_single[] == [:Julia]
    gui.selection_single[] = [:Rust]
    gui.selection_single[] = [:French, :Julia]
    @test gui.selection_single[] == [:French]
    gui.selection_single[] = []

    gui.selection[] = [:German, :Rust, :French]
    gui.selection_single[] = [:French]


    # change options
    gui.options[] = [
        OptionGroup("Languages", [:French, :Spanish, :German]),
    ]
    @test gui.selection[] == [:German, :French]
    @test gui.selection_single[] == [:French]


    # change options but keep number the same
    gui.selection[] = [:French]
    gui.selection_single[] = [:French]
    gui.options[] = [
        OptionGroup("Languages", [:French, :Spanish]),
    ]
    @test gui.selection[] == [:French]
    @test gui.selection_single[] == [:French]
    close(disp)
end
