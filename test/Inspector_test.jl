using NetworkDynamics
using NetworkDynamics: SymbolicIndex, SII
using NetworkDynamicsInspector
using Bonito
using WGLMakie
using WGLMakie.Makie.ColorSchemes
using GraphMakie
using Graphs: SimpleGraph
using OrdinaryDiffEqTsit5
using Graphs: Graphs
using OrderedCollections
include(joinpath(pkgdir(NetworkDynamics), "test", "ComponentLibrary.jl"))

sol = let
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
end;

# ENV["JULIA_DEBUG"] = NetworkDynamicsInspector
# ENV["JULIA_DEBUG"] = ""

app = (;
    sol = Observable{Any}(sol),
    t = Observable{Float64}(0.0),
    tmin = Observable{Float64}(sol.t[begin]),
    tmax = Observable{Float64}(sol.t[end]),
    active_tsplot = Observable{String}("a"),
    graphplot = (;
        nstate = Observable{Vector{Symbol}}([:θ]),
        estate = Observable{Vector{Symbol}}([:P]),
        nstate_rel = Observable{Bool}(false),
        estate_rel = Observable{Bool}(false),
        ncolorrange = Observable{Tuple{Float32,Float32}}((-1.0, 1.0)),
        ncolorscheme = Observable{ColorScheme}(ColorSchemes.coolwarm),
        ecolorrange = Observable{Tuple{Float32,Float32}}((-1.0, 1.0)),
        ecolorscheme = Observable{ColorScheme}(ColorSchemes.coolwarm),
        _selcomp = Observable{Vector{SymbolicIndex}}(SymbolicIndex[]),
        _hoverel = Observable{Union{EIndex{Int,Nothing},VIndex{Int,Nothing},Nothing}}(nothing),
        _lastclickel = Observable{Union{EIndex{Int,Nothing},VIndex{Int,Nothing},Nothing}}(nothing),
    ),
    tsplots = Observable{Any}(OrderedDict(
        "a" => (;
            selcomp = Observable{Vector{SymbolicIndex}}(SymbolicIndex[]),
            states = Observable{Vector{Symbol}}(Symbol[]),
            rel = Observable{Bool}(false),
        ),
        "b" => (;
            selcomp = Observable{Vector{SymbolicIndex}}(SymbolicIndex[]),
            states = Observable{Vector{Symbol}}(Symbol[]),
            rel = Observable{Bool}(false),
        ),
    ))
);

let
    isfile(WGLMakie.WGL.bundle_file) && rm(WGLMakie.WGL.bundle_file)
    bonitobundle = joinpath(pkgdir(Bonito), "js_dependencies", "Bonito.bundled.js")
    ispath(bonitobundle) && rm(bonitobundle)
    _app = App() do session
        @info "start new session"
        WGLMakie.activate!(resize_to=:parent)
        NetworkDynamicsInspector.clear_obs!(app)

        resize_gp = js"""
        const graphplotCard = document.querySelector(".graphplot-card");

        // Function to update the width dynamically
        function updateResizeWithGpWidth() {
            const graphplotWidth = getComputedStyle(graphplotCard).width;  // Get the computed width of .graphplot-card
            const resizeWithGpElements = document.querySelectorAll(".resize-with-gp");

            // Dynamically update the styles for .resize-with-gp elements
            let styleTag = document.getElementById("dynamic-resize-style");

            // Create the <style> tag if it doesn't exist
            if (!styleTag) {
                styleTag = document.createElement("style");
                styleTag.id = "dynamic-resize-style";
                document.head.appendChild(styleTag);
            }

            // Generate new CSS rule based on the width of .graphplot-card
            styleTag.innerHTML = `.resize-with-gp { width: ${graphplotWidth} !important; }`;

            // Manually trigger the resize event on the window
            // const resizeEvent = new Event('resize');
            // window.dispatchEvent(resizeEvent);
        };

        // Use ResizeObserver for live resizing feedback
        const updateResizeWithGpWidth_throttled = Bonito.throttle_function(updateResizeWithGpWidth, 10);
        const resizeObserver = new ResizeObserver(updateResizeWithGpWidth_throttled);
        resizeObserver.observe(graphplotCard);

        // Initial update
        updateResizeWithGpWidth();
        """
        Bonito.evaljs(session, resize_gp)

        DOM.div(
            NetworkDynamicsInspector.APP_CSS,
            DOM.div(
                NetworkDynamicsInspector.graphplot_card(app, session),
                NetworkDynamicsInspector.gpstate_control_card(app, :vertex),
                NetworkDynamicsInspector.gpstate_control_card(app, :edge),
                NetworkDynamicsInspector.element_info_card(app, session),
                class="graphplot-col"
            ),
            DOM.div(
                NetworkDynamicsInspector.timeslider_card(app),
                NetworkDynamicsInspector.timeseries_cards(app, session),
                class="timeseries-col"
            ),
            class="maingrid"
        )
    end;
    serve_app(_app)
end
