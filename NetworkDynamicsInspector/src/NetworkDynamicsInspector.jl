module NetworkDynamicsInspector

using Bonito: Bonito, @js_str, Asset, CSS, Styles,
              Grid, Card, DOM, Session, ES6Module
using NetworkDynamics: NetworkDynamics, SII, EIndex, VIndex, Network,
                       get_metadata, has_metadata, get_position, has_position,
                       obssym, psym, sym, extract_nw

using Graphs: nv, ne
using WGLMakie: WGLMakie
using WGLMakie.Makie: Makie, @lift, MouseEvent, Point2f, with_theme,
                      lines!, vlines!, Theme, Figure, Colorbar, Axis,
                      xlims!, ylims!, autolimits!, hidespines!, hidedecorations!,
                      register_interaction!, MouseEventTypes, Consume, events,
                      mouseposition

using GraphMakie: GraphMakie, EdgeClickHandler, EdgeHoverHandler,
                  NodeClickHandler, NodeHoverHandler, graphplot!
using GraphMakie.NetworkLayout: Stress

using OrderedCollections: OrderedDict
using Observables: Observables, Observable, on, onany
using Colors: Colors, @colorant_str, RGB, color
using ColorSchemes: ColorSchemes, ColorScheme
using SciMLBase: SciMLBase

# defined based on the julia version
using NetworkDynamics: AnnotatedIOBuffer, AnnotatedString, @styled_str

include("utils.jl")

const ASSETS = joinpath(dirname(@__DIR__),"assets")
download_assets() # download assets if not present
const APP_CSS = Asset(joinpath(ASSETS, "app.css"))
const TOMSELECT_ESS = ES6Module(joinpath(ASSETS, "tomselect.js"))
const TOMSELECT_CSS = Asset(joinpath(ASSETS, "tomselect.css"))
const ELECTRON_JS = Asset(joinpath(ASSETS, "electron.js"))

include("widgets.jl")
include("graphplot.jl")
include("timeseries.jl")
export BrowserDisp, ServerDisp, ElectronDisp
include("serving.jl")

const SymbolicCompIndex = Union{VIndex{Int,Nothing}, EIndex{Int,Nothing}}

export inspect, dump_app_state
export set_sol!, set_state!, set_graphplot!, set_timeseries!, define_timeseries!
include("appstate.jl")

const APPSTATE = Ref{Union{Nothing,AppState}}(nothing)
const SESSION = Ref{Union{Nothing,Session}}(nothing)
const CURRENT_DISPLAY = Ref{Union{Nothing, NDIDisplay}}(BrowserDisp())
const CURRENT_WEBAPP = Ref{Union{Nothing, Bonito.App}}(nothing)

function get_webapp(app)
    webapp = Bonito.App() do session
        this_session = isnothing(session.parent) ? session : session.parent
        if SESSION[] != this_session
            @info "New GUI Session started"
            close_session(SESSION[])
            SESSION[] = this_session
        else
            @info "GUI Session updated"
        end

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
            window.dispatchEvent(new Event('resize'));
        };

        // Use ResizeObserver for live resizing feedback
        const updateResizeWithGpWidth_throttled = Bonito.throttle_function(updateResizeWithGpWidth, 100);
        const resizeObserver = new ResizeObserver(updateResizeWithGpWidth_throttled);
        resizeObserver.observe(graphplotCard);

        // Initial update
        updateResizeWithGpWidth();
        """
        Bonito.evaljs(session, resize_gp)

        DOM.div(
            APP_CSS,
            ELECTRON_JS,
            DOM.div(
                graphplot_card(app, session),
                gpstate_control_card(app, :vertex),
                gpstate_control_card(app, :edge),
                element_info_card(app, session),
                class="graphplot-col"
            ),
            timeseries_col(app, session),
            class="maingrid"
        )
    end;
end

"""
    inspect(sol; restart=false, reset=false, display=nothing)

Main entry point for gui. Starts the server and serves the app for
soution `sol`.

- `restart`: If `true`, the display will be restartet (i.e. new Electron window, new server or new Browser tab)
- `reset`: If `true`, reset the appstate with the new solution `sol`.
- `display=CURRENT_DISPLAY[]`: Can be `BrowserDisp()`, `ServerDisp()` or `ElectronDisp()`.
   Per default, the current display will be used (defaults to`BrowserDisp()`).
"""
function inspect(sol; restart=false, reset=false, display=CURRENT_DISPLAY[])
    if typeof(display) != typeof(CURRENT_DISPLAY[]) || restart
        close_session(SESSION[])
        close_display(CURRENT_DISPLAY[]; strict=true) # make sure to completely close electron
        CURRENT_WEBAPP[] = nothing
    end

    CURRENT_DISPLAY[] = display

    appstate = if reset || isnothing(APPSTATE[])
        AppState(sol)
    else
        APPSTATE[].sol[] = sol
        APPSTATE[]
    end

    webapp = if restart || isnothing(CURRENT_WEBAPP[]) || appstate != APPSTATE[]
        get_webapp(appstate)
    else
        CURRENT_WEBAPP[]
    end

    serve_app(display, webapp)

    APPSTATE[] = appstate
    CURRENT_WEBAPP[] = webapp

    nothing
end

function apptheme()
    Theme(
        fontsize=10,
        Lines = (;
            linewidth = 2,
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

    help = HoverHelp(html"""
    <ul>
    <li>Adjust the time for the graph coloring.</li>
    <li>Use <strong>Arrow Keys</strong> to adjust the time in small increments.</li>
    <li>Use <strong>Shift + Arrow Keys</strong> to adjust the time in larger increments.</li>
    <li>Use second slider to adjust the time window (i.e. to focus on a specific timeframe in the plots below)</li>
    <ul>
    """)

    Card(
        [Grid(
            DOM.span("Time"), t_slider, RoundedLabel(t_slider.value_r),
            RoundedLabel(tw_slider.value_l; style=Styles("text-align"=>"right")),
            tw_slider,
            RoundedLabel(tw_slider.value_r; style=Styles("text-align"=>"left"));
            columns="70px auto 70px",
            justify_content="begin",
            align_items="center",
        ), help];
        class="bonito-card timeslider-card"
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

    help = HoverHelp(html"""
    Show details on last clicked element.
    <strong>Shift + Click</strong> element in graphplot to update details pane
    without adding/removing timeseries.
    """)

    Card(
        [DOM.div(;id="element-info-box"), help];
        class="bonito-card element-info-card resize-with-gp"
    )
end



end # module NetworkDynamicsInspector
