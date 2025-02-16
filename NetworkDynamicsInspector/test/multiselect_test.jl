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
    selection = Observable{Vector{Symbol}}(Symbol[:Rust]),
)

SERVER = Ref{Any}()
let
    app = App(;) do session
        NetworkDynamicsInspector.clear_obs!(gui)

        ms = MultiSelect(gui.options, gui.selection; placeholder="pick language", T=Symbol)
        @info "compare" ms.selection gui.selection ms.selection===gui.selection

        return DOM.div(
            ms
        )
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

gui.selection[] = [:Rust, :French, :Julia, :foo]

gui.options[] = [
        OptionGroup("Programming Languages", [:Julia, :Rust, :Java]),
    ]

gui.selection[] = [:Julia, :Rust]
