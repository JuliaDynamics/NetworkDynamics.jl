module NetworkDynamicsInspector

using Bonito: Bonito, @js_str, Asset, CSS, Styles,
              Grid, Card, DOM, Session
using NetworkDynamics: NetworkDynamics, SII, EIndex, VIndex, Network,
                       get_metadata, has_metadata, get_position, has_position,
                       obssym, psym, sym, extract_nw

using Graphs: nv, ne
using WGLMakie: WGLMakie
using WGLMakie.Makie: Makie, @lift, MouseEvent, Point2f, with_theme,
                      Cycle, lines!, vlines!, Theme, Figure, Colorbar, Axis,
                      xlims!, ylims!, autolimits!, hidespines!, hidedecorations!,
                      register_interaction!, MouseEventTypes, Consume, events,
                      mouseposition

using GraphMakie: GraphMakie, EdgeClickHandler, EdgeHoverHandler,
                  NodeClickHandler, NodeHoverHandler, graphplot!
using GraphMakie.NetworkLayout: Stress

using OrderedCollections: OrderedDict
using Observables: Observable, on, onany
using Colors: Colors, @colorant_str, RGB, color
using ColorSchemes: ColorSchemes, ColorScheme
using SciMLBase: SciMLBase

# defined based on the julia version
using NetworkDynamics: AnnotatedIOBuffer, AnnotatedString

include("utils.jl")

const ASSETS = joinpath(dirname(@__DIR__),"assets")
download_assets() # download assets if not present

const JQUERY = Asset(joinpath(ASSETS, "jquery.js"))
const SELECT2_CSS = Asset(joinpath(ASSETS, "select2.css"))
const SELECT2_JS = Asset(joinpath(ASSETS, "select2.js"))
const APP_CSS = Asset(joinpath(ASSETS, "app.css"))

include("widgets.jl")
include("graphplot.jl")
include("timeseries.jl")

const SymbolicCompIndex = Union{VIndex{Int,Nothing}, EIndex{Int,Nothing}}

export inspect

@kwdef struct GraphPlot
    nstate::Observable{Vector{Symbol}} = [:nothing]
    estate::Observable{Vector{Symbol}} = [:nothing]
    nstate_rel::Observable{Bool} = false
    estate_rel::Observable{Bool} = false
    ncolorrange::Observable{Tuple{Float32,Float32}} = (-1.0, 1.0)
    ncolorscheme::Observable{ColorScheme} = ColorSchemes.coolwarm
    ecolorrange::Observable{Tuple{Float32,Float32}} = (-1.0, 1.0)
    ecolorscheme::Observable{ColorScheme} = ColorSchemes.coolwarm
    _selcomp::Observable{Vector{SymbolicCompIndex}} = SymbolicCompIndex[]
    _hoverel::Observable{Union{SymbolicCompIndex,Nothing}} = nothing
    _lastclickel::Observable{Union{SymbolicCompIndex,Nothing}} = nothing
end
function GraphPlot(sol)
    nw = extract_nw(sol)
    estate = [_most_common_output_state(nw.im.edgem)]
    nstate = [_most_common_output_state(nw.im.vertexm)]
    GraphPlot(; nstate, estate)
end
function _most_common_output_state(models)
    states = mapreduce(NetworkDynamics.outsym_flat, vcat, models)
    unique_states = unique(states)
    counts = map(unique_states) do s
        count(isequal(s), states)
    end
    unique_states[argmax(counts)]
end

@kwdef struct TimeseriesPlot
    selcomp::Observable{Vector{SymbolicCompIndex}} = SymbolicCompIndex[]
    states::Observable{Vector{Symbol}} = Symbol[]
    rel::Observable{Bool} = false
end

struct AppState
    sol::Observable{Any}
    t::Observable{Float64}
    tmin::Observable{Float64}
    tmax::Observable{Float64}
    active_tsplot::Observable{String}
    graphplot::GraphPlot
    tsplots::Observable{OrderedDict{String, TimeseriesPlot}}
end
function AppState(sol::SciMLBase.AbstractODESolution)
    t = sol.t[begin]
    tmin = sol.t[begin]
    tmax = sol.t[end]
    graphplot = GraphPlot(sol)
    tsplots = OrderedDict(
        "ts-1" => TimeseriesPlot(),
        "ts-2" => TimeseriesPlot()
    )
    active_tsplot = "ts-1"
    AppState(sol, t, tmin, tmax, active_tsplot, graphplot, tsplots)
end

const APPSTATE = Ref{Union{Nothing,AppState}}(nothing)
const SERVER = Ref{Any}(nothing)
const SESSION = Ref{Union{Nothing,Session}}(nothing)

function reset!(sol=nothing)
    if isnothing(APPSTATE[]) && isnothing(sol)
        error("No appstate to reset")
    else
        clear_obs!(APPSTATE[])
        if isnothing(sol)
            APPSTATE[] = AppState(APPSTATE[].sol[])
        else
            APPSTATE[] = AppState(sol)
        end
    end
end

server_running() = !isnothing(SERVER[]) && Bonito.HTTPServer.isrunning(SERVER[])

function stop_server!()
    if !isnothing(SESSION[])
        if Base.isopen(SESSION[])
            @info "Close running session..."
            close(SESSION[])
        end
        SESSION[] = nothing
    end
    if server_running()
        @info "Stop running server..."
        close(SERVER[])
    end
end

function start_server!(restart=true)
    if server_running() && !restart
        error("Server already running")
    end
    stop_server!()

    if isnothing(APPSTATE[])
        error("No appstate to restart")
    end

    app = APPSTATE[]

    webapp = Bonito.App() do session
        @info "New GUI Session started"
        if !isnothing(SESSION[]) && Base.isopen(SESSION[])
            @info "Close previous session..."
            close(SESSION[])
        end
        SESSION[] = session

        WGLMakie.activate!(resize_to=:parent)
        clear_obs!(app)

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
            APP_CSS,
            DOM.div(
                graphplot_card(app, session),
                gpstate_control_card(app, :vertex),
                gpstate_control_card(app, :edge),
                element_info_card(app, session),
                class="graphplot-col"
            ),
            DOM.div(
                timeslider_card(app),
                timeseries_cards(app, session),
                class="timeseries-col"
            ),
            class="maingrid"
        )
    end;

    SERVER[] = Bonito.Server(webapp, "localhost", 8080)
    url = SERVER[].url
    port = SERVER[].port
    @info "Visit $url:$port to launch App"
end

"""
    inspect(sol; restart=false, reset=false)

Main entry point for gui. Starts the server and serves the app for
soution `sol`.

- `restart`: If `true`, stop the server if it is running and start a new one.
- `reset`: If `true`, reset the appstate with the new solution `sol`.
"""
function inspect(sol; restart=false, reset=false)
    if restart && server_running()
        stop_server!()
    end
    if isnothing(APPSTATE[]) || reset
        reset!(sol)
    else
        APPSTATE[].sol[] = sol
    end

    if !server_running()
        start_server!()
    else
        @info "App still served at $(SERVER[].url):$(SERVER[].port)"
    end
    nothing
end

function apptheme()
    Theme(
        fontsize=10,
        palette = (;
           linestyle = [:solid, :dot, :dash, :dashdot, :dashdotdot],
        ),
        Lines = (;
            cycle = Cycle([:color, :linestyle], covary=true),
            linewidth = 3,
        )
    )
end

function timeslider_card(app)
    on(app.sol) do _sol
        app.tmin[] = _sol.t[begin]
        app.tmax[] = _sol.t[end]
    end

    trange_sol = @lift ($(app.sol).t[begin], $(app.sol).t[end])
    tw_slider = ContinuousSlider(trange_sol, app.tmin, app.tmax)
    twindow = @lift ($(app.tmin), $(app.tmax))
    # if window changes, keep t[] inside
    onany(app.tmin, app.tmax; update=true) do tmin, tmax
        _t = clamp(app.t[], tmin, tmax)
        if _t != app.t[]
            @debug "app.tmin, app.tmax => clamp app.t[]"
            app.t[] = _t
        end
        nothing
    end
    t_slider = ContinuousSlider(twindow, app.t; arrowkeys=true)
    Card(
        Grid(
            DOM.span("Time"), t_slider, RoundedLabel(t_slider.value_r),
            RoundedLabel(tw_slider.value_l; style=Styles("text-align"=>"right")),
            tw_slider,
            RoundedLabel(tw_slider.value_r; style=Styles("text-align"=>"left"));
            columns="70px auto 70px",
            justify_content="begin",
            align_items="center",
        );
        class="timeslider-card"
    );
end

function element_info_card(app, session)
    eltext = Observable{AnnotatedString}("")
    onany_throttled(app.graphplot._lastclickel, app.t; update=true, delay=0.1) do el, t
        if isnothing(el)
            eltext[] = ""
        else
            buf = AnnotatedIOBuffer()
            NetworkDynamics.dump_state(buf, app.sol[], t, el)
            s = read(seekstart(buf), AnnotatedString)
            htmlbuf = IOBuffer()
            show(htmlbuf, MIME"text/html"(), s)
            html = replace(String(take!(htmlbuf)), r"\n" => "<br>")

            t = NetworkDynamics.str_significant(app.t[]; sigdigits=3)
            fake_prompt = "<span class=julia-prompt>julia&gt; </span><span class=julia-command>dump_state(sol, $t, $(repr(el)))</span>"
            eltext[] = "<pre>$fake_prompt<br><br>$html</pre>"
        end
    end

    js = js"""
        const div = document.getElementById("element-info-box");

        function replaceDivContentWithHTML(htmlString) {
            if (!div) return;
            div.innerHTML = htmlString;
        }

        $(eltext).on(replaceDivContentWithHTML);
        replaceDivContentWithHTML($(eltext).value);
    """
    Bonito.evaljs(session, js)

    Card(
        DOM.div(;id="element-info-box");
        class="element-info-card resize-with-gp"
    )
end



end # module NetworkDynamicsInspector
