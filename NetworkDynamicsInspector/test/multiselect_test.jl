using Bonito
using NetworkDynamicsInspector
using NetworkDynamicsInspector: OptionGroup, MultiSelect

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

SERVER = Ref{Any}()
let
    app = App(;) do session
        NetworkDynamicsInspector.clear_obs!(gui)

        ms1 = MultiSelect(gui.options, gui.selection; placeholder="pick language", T=Symbol)
        ms2 = MultiSelect(gui.options, gui.selection_single; placeholder="pick language", multi=false, T=Symbol)

        return Grid(ms1, ms2)
    end;
    try
        println("Close existing")
        close(SERVER[])
    catch
    end
    println("Start server")
    SERVER[] = Bonito.Server(app, "0.0.0.0", 8080)
    # Bonito.update_app!
end

@test gui.selection[] == [:Rust, :French]
gui.selection[] = [:Rust, :French, :Julia, :foo]
gui.selection[] = []

@test gui.selection_single[] == [:Julia]
gui.selection_single[] = [:Rust]
gui.selection_single[] = [:French, :Julia]
@test gui.selection_single[] == [:French]
gui.selection_single[] = []

gui.selection[] = [:Rust, :French]
gui.selection_single[] = [:French]
# change options
gui.options[] = [
        OptionGroup("Programming Languages", [:Julia, :Rust, :Java]),
    ]
@test gui.selection[] == [:Rust]
@test gui.selection_single[] == []

# TODO
# - css file for app



#

gui.selection[] = [:Julia, :Rust]
