using NetworkDynamicsInspector
using NetworkDynamics
using Bonito
using WGLMakie
using WGLMakie.Makie.ColorSchemes
using GraphMakie
using Graphs: SimpleGraph
using OrdinaryDiffEqTsit5
using Graphs: Graphs

sol = let
    include(joinpath(pkgdir(NetworkDynamics), "test", "ComponentLibrary.jl"))

    g = SimpleGraph([0 1 1 0 1;
                        1 0 1 1 0;
                        1 1 0 1 0;
                        0 1 1 0 1;
                        1 0 0 1 0])
    vs = [Lib.swing_mtk() for _ in 1:5];
    set_default!(vs[1], :Pmech, -1)
    set_default!(vs[2], :Pmech, 1.5)
    set_default!(vs[3], :Pmech, -1)
    set_default!(vs[4], :Pmech, -1)
    set_default!(vs[5], :Pmech, 1.5)
    ls = [Lib.line_mtk() for _ in 1:7];
    nw = Network(g, vs, ls)
    sinit = NWState(nw)
    s0 = find_fixpoint(nw)
    set_defaults!(nw, s0)

    # set_position!(vs[1], (0.0, 0.0))
    set_marker!(vs[1], :dtriangle)
    set_marker!(vs[2], :utriangle)
    set_marker!(vs[3], :dtriangle)
    set_marker!(vs[4], :dtriangle)
    set_marker!(vs[5], :utriangle)

    cond = ComponentCondition([:P, :₋P, :srcθ], [:limit, :K]) do u, p, t
        abs(u[:P]) - p[:limit]
    end
    affect = ComponentAffect([],[:active]) do u, p, ctx
        @info "Trip line $(ctx.eidx) between $(ctx.src) and $(ctx.dst) at t=$(ctx.t)"
        p[:active] = 0
    end
    cb = ContinousComponentCallback(cond, affect)
    set_callback!.(ls, Ref(cb))

    tripfirst = PresetTimeComponentCallback(1.0, affect) # reuse the same affect
    add_callback!(nw[EIndex(5)], tripfirst)

    nwcb = NetworkDynamics.get_callbacks(nw);
    s0 = NWState(nw)
    prob = ODEProblem(nw, uflat(s0), (0,6), copy(pflat(s0)), callback=nwcb)
    sol = solve(prob, Tsit5());
end

app = (;
    sol = Observable{Any}(sol),
    t = Observable{Float64}(0.0),
    tmin = Observable{Float64}(sol.t[begin]),
    tmax = Observable{Float64}(sol.t[end]),
    sel_nodes = Observable{Vector{Int}}(Int[]),
    sel_edges = Observable{Vector{Int}}(Int[]),
    graphplot = (;
        nstate = Observable{Union{Symbol,Nothing}}(:θ),
        estate = Observable{Union{Symbol,Nothing}}(:P),
        nstate_rel = Observable{Bool}(false),
        estate_rel = Observable{Bool}(false),
        ncolorrange = Observable{Tuple{Float32,Float32}}((-1.0, 1.0)),
        ncolorscheme = Observable{ColorScheme}(ColorSchemes.coolwarm),
        ecolorrange = Observable{Tuple{Float32,Float32}}((-1.0, 1.0)),
        ecolorscheme = Observable{ColorScheme}(ColorSchemes.coolwarm),
    )
);

gpfig = Ref{Any}()
App() do session
    WGLMakie.activate!(resize_to=:parent)
    NetworkDynamicsInspector.clear_obs!(app)
    gpfig[] = fig
    Grid(
        NetworkDynamicsInspector.graphplot_card(app; height="400px"),
        NetworkDynamicsInspector.nodestate_control_card(app),
        NetworkDynamicsInspector.timeslider_card(app),
        columns="100%",
        width="500px"
    )
end

app.t[] = 4.0

app.tmin[] = NaN
app.tmin[] = 1.0

app.tmin[]
app.tmax[]

fig = Figure()
ax = Axis(fig[1,1])
graphplot!(ax, Graphs.smallgraph(:karate))



using Bonito
using NetworkDynamicsInspector
using NetworkDynamicsInspector: OptionGroup

options = [
    OptionGroup("Programming Languages", ["Julia", :Rust, "Java"]),
    OptionGroup("Languages", ["French", "Spanish", "German"]),
    :car,
    1.0
]
jsoptions = NetworkDynamicsInspector.options_to_jsoptions(options)

sel = NetworkDynamicsInspector.jsselection_to_selection(options, 1:8)
NetworkDynamicsInspector.selection_to_jsselection(options, sel) == 1:8




multiselect = (;
    options = Observable([
        # OptionGroup("Programming Languages", ["Julia", "Rust", "Java"]),
        # OptionGroup("Languages", ["French", "Spanish", "German"]),
        (;label="Programming Languages", options=[(; id=1, text="Julia"), (;id=2,text="Rust")]),
        (;label="Languages", options=[(; id=3, text="Spanish"), (;id=4,text="French")]),
    ]),
    selection = Observable{Vector{Int}}([]),
    placeholder="Select state(s) to plot"
)

using Bonito
using NetworkDynamicsInspector
using NetworkDynamicsInspector: OptionGroup, MultiSelect

gui = (;
    options = Observable{Vector{Union{Symbol,OptionGroup{Symbol}}}}([
        OptionGroup("Programming Languages", [:Julia, :Rust, :Java]),
        OptionGroup("Languages", [:French, :Spanish, :German]),
    ]),
    selection = Observable{Vector{Symbol}}(Symbol[:Rust]),
)

SERVER = Ref{Any}()
let
    app = App(;) do session
        NetworkDynamicsInspector.clear_obs!(gui)

        ms = MultiSelect(gui.options, gui.selection; placeholder="pick language", T=Symbol)
        @info "compare" ms.selection gui.selection ms.selection===gui.selection

        on(gui.selection) do sel
            @show sel
        end

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
push!(gui.options[], :Baz)
gui.options[]

gui.selection[] = [:Rust, :Julia]
gui.selection[] = [:Julia, :Rust]

NetworkDynamicsInspector.selection_to_jsselection(gui.options[], gui.selection[])

gui.selection[] = [:rust, :julia]

multiselect.selection[] = [1,2]
multiselect.selection[] = [2,1]

multiselect.selection[] = [1,2,3]



close(server)
server




# Start the server
server = Server(app, "0.0.0.0", 8080)
wait(server)


d = DOM.div("foo")
jsrender()


    # jquery = Asset("https://cdn.jsdelivr.net/npm/jquery@3.7.1/dist/jquery.min.js")
        # DOM.head(jquery),

using Bonito
app = App() do session
    return DOM.html(
        DOM.head(
            DOM.meta(name="viewport", content="width=device-width, initial-scale=1.0"),
            DOM.meta(charset="utf-8"),
            Bonito.TailwindCSS,
        ),
        DOM.body(
            DOM.div("Hello")
        )
    )
end;
server = Bonito.Server(app, "0.0.0.0", 8080)

using Bonito
app = App() do session
    return DOM.html(
        DOM.head(),
        DOM.body(
            DOM.div("Hello")
        )
    )
end;
server = Bonito.Server(app, "0.0.0.0", 8080)
close(server)
