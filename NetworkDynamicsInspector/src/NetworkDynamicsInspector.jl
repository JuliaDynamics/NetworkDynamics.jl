module NetworkDynamicsInspector

using Bonito
using NetworkDynamics
using Bonito.Observables
using Bonito: Grid, @js_str, onload, jsrender
using WGLMakie.Makie
using WGLMakie.Makie.Colors
using WGLMakie.Makie.ColorSchemes
using NetworkDynamics: extract_nw, SymbolicIndex
using NetworkDynamics: SII
using Graphs: nv, ne
using GraphMakie
using GraphMakie.NetworkLayout
using Bonito.Hyperscript
using Bonito.Tables.OrderedCollections

const JQUERY = Asset("https://cdn.jsdelivr.net/npm/jquery@3.7.1/dist/jquery.min.js")
const SELECT2_CSS = Asset("https://cdn.jsdelivr.net/npm/select2@4.1.0-rc.0/dist/css/select2.min.css")
const SELECT2_JS = Asset("https://cdn.jsdelivr.net/npm/select2@4.1.0-rc.0/dist/js/select2.min.js")
const APP_CSS = Asset(joinpath(pkgdir(NetworkDynamicsInspector), "assets", "app.css"))

export ContinuousSlider, RoundedLabel
include("widgets.jl")

export wrap_assets
include("utils.jl")

function apptheme()
    Theme(
        fontsize=10,
        palette = (;
           linestyle = [:solid, :dot, :dash, :dashdot, :dashdotdot],
        ),
        Lines = (;
            cycle = Cycle([:color, :linestyle], covary=true)
        )
    )
end

function graphplot_card(app, session)
    nw = map!(extract_nw, Observable{Network}(), app.sol)
    NV = nv(nw[])
    NE = ne(nw[])

    node_marker = map(nw) do _nw
        @debug "GP: app.ng => node_marker"
        if any(m->has_metadata(m, :marker), _nw.im.vertexm)
            [has_metadata(m, :marker) ? get_metadata(m, :marker) : :circle for m in _nw.im.vertexm]
        else
            markerset = [:circle, :rect, :utriangle, :cross, :diamond, :dtriangle, :pentagon, :xcross]
            getmarker(i) = markerset[(i-1) % length(markerset) + 1]
            map(1:NV) do i
                group = findfirst(batch -> i ∈ batch.indices, _nw.vertexbatches)
                getmarker(group)
            end
        end
    end

    node_state = Observable(Vector{Float32}(undef, NV))
    onany(app.sol, app.t, app.graphplot.nstate, app.graphplot.nstate_rel) do _sol, _t, _state, _rel
        @debug "GP: app.sol, app.t, app.gp.nstate, app.gp.nstate_rel => node_state"
        if length(_state) == 0
            fill!(node_state[], NaN)
        elseif length(_state) == 1
            idxs = VIndex.(1:NV, only(_state))
            _gracefully_extract_states!(node_state[], _sol, _t, idxs, _rel)
        else
            error("Received more than one node state to plot...")
        end
        notify(node_state)
        nothing
    end;

    edge_state = Observable(Vector{Float32}(undef, NE))
    onany(app.sol, app.t, app.graphplot.estate, app.graphplot.estate_rel) do _sol, _t, _state, _rel
        @debug "GP: app.sol, app.t, app.gp.estate, app.gp.estate_rel => edge_state"
        if length(_state) == 0
            fill!(edge_state[], NaN)
        elseif length(_state) == 1
            idxs = EIndex.(1:NE, only(_state))
            _gracefully_extract_states!(edge_state[], _sol, _t, idxs, _rel)
        else
            error("Received more than one edge state to plot...")
        end
        notify(edge_state)
        nothing
    end;

    node_color = Observable(Vector{RGB{Float64}}(undef, NV))
    onany(node_state, app.graphplot.ncolorrange, app.graphplot.ncolorscheme) do statevec, range, scheme
        @debug "GP: node_state, app.gp.ncolorrange, app.gp.ncolorrange => node_color"
        for i in 1:NV
            node_color[][i] = isnan(statevec[i]) ? RGB(0,0,0) : get(scheme, statevec[i], range)
        end
        notify(node_color)
        nothing
    end

    edge_color = Observable(Vector{RGB{Float64}}(undef, NE))
    onany(edge_state, app.graphplot.ecolorrange, app.graphplot.ecolorscheme) do statevec, range, scheme
        @debug "GP: edge_state, app.gp.ecolorrange, app.gp.ecolorrange => edge_color"
        for i in 1:NE
            edge_color[][i] = isnan(statevec[i]) ? RGB(0,0,0) : get(scheme, statevec[i], range)
        end
        notify(edge_color)
        nothing
    end

    notify(app.t) # trigger updates

    SMALL = 30
    BIG = 50
    node_size = Observable(fill(SMALL, NV))
    THIN = 5
    THICK = 8
    edge_width = Observable(fill(THIN, NE))

    onany(app.graphplot._selcomp; update=true) do selcomp
        @debug "GP: Sel comp => node_size, edge_width"
        fill!(node_size[], SMALL)
        fill!(edge_width[], THIN)
        for s in selcomp
            if s isa VIndex
                node_size[][s.compidx] = BIG
            else
                edge_width[][s.compidx] = THICK
            end
        end
        notify(node_size)
        notify(edge_width)
        nothing
    end

    g = @lift $(nw).im.g

    layout = map(nw) do _nw
        @debug "GP: app.nw => layout"
        pin = Dict{Int,Point2f}()
        for i in 1:NV
            vm = _nw[VIndex(i)]
            if has_position(vm)
                pin[i] = get_position(vm)
            end
        end
        Stress(; pin)
    end

    fig, ax = with_theme(apptheme()) do
        fig = Figure(; figure_padding=0)
        ax = Axis(fig[1,1])
        graphplot!(ax, g; layout, node_marker, node_color, node_size, edge_color, edge_width,
            node_attr=(;colorrange=app.graphplot.ncolorrange, colormap=app.graphplot.ncolorscheme),
            edge_attr=(;colorrange=app.graphplot.ncolorrange, colormap=app.graphplot.ncolorscheme))

        hidespines!(ax)
        hidedecorations!(ax)
        fig, ax
    end
    xratio = Ref{Float64}(1.0)
    yratio = Ref{Float64}(1.0)
    on(ax.scene.viewport) do lims
        @debug "GP: viewport => adapt xy scaling"
        adapt_xy_scaling!(xratio, yratio, ax)
        nothing
    end

    ####
    #### Interactions
    ####
    tooltip_text = Observable{String}("")

    js = js"""
    const gpcard = document.querySelector(".graphplot-card");

    // Create tooltip element
    const tooltip = document.createElement("div");
    tooltip.classList.add("tooltip"); // Apply styles from CSS
    document.body.appendChild(tooltip);

    // let tooltipTimeout; // Store timeout reference
    let latestX = 0, latestY = 0; // Store latest cursor position
    let isUpdating = false;

    gpcard.addEventListener("mousemove", (event) => {
        latestX = event.clientX;
        latestY = event.clientY;

        if (!isUpdating && tooltip.style.display === "block") {
            isUpdating = true;
            requestAnimationFrame(() => {
                tooltip.style.left = `${latestX}px`;
                tooltip.style.top = `${latestY}px`;
                isUpdating = false;
            });
        }
    });

    // Show tooltip after delay when idx > 0
    $(tooltip_text).on((text) => {
        // clearTimeout(tooltipTimeout); // Prevent multiple pending tooltips

        if (text.length > 0) {
            gpcard.style.cursor = "pointer";

            // tooltipTimeout = setTimeout(() => { // Delay tooltip display
                tooltip.textContent = text;
                tooltip.style.display = "block";
                tooltip.style.left = `${latestX}px`;
                tooltip.style.top = `${latestY}px`;
            // }, 300); // 0.3s delay
        } else {
            gpcard.style.cursor = "default";
            tooltip.style.display = "none"; // Hide tooltip immediately
            // clearTimeout(tooltipTimeout); // Cancel any pending tooltip show
        }
    });
    """
    evaljs(session, js)

    nhh = NodeHoverHandler() do hstate, idx, event, axis
        tooltip_text[] = hstate ? _sidx_to_str(VIndex(idx), app) : ""
        app.graphplot._hoverel[] = hstate ? VIndex(idx) : nothing
    end
    ehh = EdgeHoverHandler() do hstate, idx, event, axis
        tooltip_text[] = hstate ? _sidx_to_str(EIndex(idx), app) : ""
        app.graphplot._hoverel[] = hstate ? EIndex(idx) : nothing
    end

    # interactions
    clickaction = (i, type) -> begin
        idx = type == :vertex ? VIndex(i) : EIndex(i)
        selcomp = app.tsplots[][app.active_tsplot[]].selcomp
        if idx ∉ selcomp[]
            push!(selcomp[], idx)
        else
            filter!(x -> x != idx, selcomp[])
        end
        app.graphplot._lastclickel[] = idx
        notify(selcomp)
    end
    nch = NodeClickHandler((i, _, _) -> clickaction(i, :vertex))
    ech = EdgeClickHandler((i, _, _) -> clickaction(i, :edge))

    register_interaction!(ax, :nodeclick, nch)
    register_interaction!(ax, :nodehover, nhh)
    register_interaction!(ax, :edgeclick, ech)
    register_interaction!(ax, :edgehover, ehh)


    Card(fig; class="graphplot-card")
end
function _gracefully_extract_states!(vec, sol, t, idxs, rel)
    isvalid(s) = SII.is_variable(sol, s) || SII.is_parameter(sol, s) || SII.is_observed(sol, s)
    fill!(vec, NaN)
    mask = isvalid.(idxs)
    # cannot acces single time point, therefore we use [t] hack
    if rel
        vec[mask] .= reduce(-, sol([t, sol.t[begin]]; idxs=idxs[mask]).u)
    else
        vec[mask] .= only(sol([t]; idxs=idxs[mask]).u)
    end
end
function adapt_xy_scaling!(xratio, yratio, ax)
    # in theory, this should work without autolimits, but it doesnt
    autolimits!(ax)

    vp = ax.scene.viewport[]
    fl = ax.finallimits[]
    vp_aspect = reduce(/, Float64.(vp.widths))
    fl_aspect = reduce(/, Float64.(fl.widths))
    ratio = vp_aspect/fl_aspect

    potential_x_factor = xratio[]*ratio
    potential_y_factor = yratio[]/ratio

    if potential_x_factor >= 1 && potential_y_factor >= 1
        # both valid, pick smaller
        resize = potential_x_factor < potential_y_factor ? :x : :y
    elseif potential_x_factor > 1
        resize = :x
    elseif potential_y_factor > 1
        resize = :y
    else
        resize = potential_x_factor > potential_y_factor ? :x : :y
        @warn "No valid scaling factor found, chose $resize"
    end

    if resize == :x
        xratio[] = potential_x_factor
        low = Float64(ax.finallimits[].origin[1])
        hi = Float64(low + ax.finallimits[].widths[1])
        xlims!(ax, ratio*low, ratio*hi)
    else
        yratio[] = potential_y_factor
        low = Float64(ax.finallimits[].origin[2])
        hi = Float64(low + ax.finallimits[].widths[2])
        ylims!(ax, low/ratio, hi/ratio)
    end
end

function timeslider_card(app)
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

function gpstate_control_card(app, type)
    VEIndex = type == :vertex ? VIndex : EIndex
    label = type == :vertex ? "Node state:" : "Edge state:"
    stateobs = type == :vertex ? app.graphplot.nstate : app.graphplot.estate
    stateobs_rel = type ==:vertex ? app.graphplot.nstate_rel : app.graphplot.estate_rel
    colorrange = type == :vertex ? app.graphplot.ncolorrange : app.graphplot.ecolorrange
    colorscheme = type == :vertex ? app.graphplot.ncolorscheme : app.graphplot.ecolorscheme
    class = type == :vertex ? "gpstate-control-card vertex" : "gpstate-control-card edge"
    class *= " resize-with-gp"

    ####
    #### State selection
    ####
    options = Observable{Vector{OptionGroup{Symbol}}}()
    on(app.sol; update=true) do _sol
        _nw = extract_nw(_sol)
        idxs = VEIndex.(1:nv(_nw))
        options[] = gen_state_options(_nw, idxs)
        nothing
    end
    multisel = MultiSelect(options, stateobs; placeholder="Select state for coloring", multi=false, T=Symbol)
    reltoggle = ToggleSwitch(value=stateobs_rel, label="Rel to u0")
    selector = Grid(
        DOM.span(label),
        multisel,
        reltoggle;
        columns = "70px 1fr auto",
        align_items = "center"
    )

    ####
    #### Color Bar & Color Bar Scaling
    ####
    thumb_pos_cache = Dict{UInt64, Tuple{Float32,Float32}}()
    function thumb_pos_key()
        hash((app.sol[], stateobs[], stateobs_rel[]))
    end

    maxrange = Observable{Tuple{Float32,Float32}}()
    thumb_l = Observable{Float32}()
    thumb_r = Observable{Float32}()

    onany(app.sol, stateobs, stateobs_rel; update=true) do _sol, _state, _rel
        @debug "GP: $type state sel: app.sol, app.gp.state, app.gp.state_rel => $type color maxrange"
        if length(_state) == 0
            maxrange[] = (-1.0, 1.0)
        elseif length(_state) == 1
            idxs = VEIndex.(1:nv(extract_nw(_sol)), only(_state))
            r = Float32.(_maxrange(_sol, idxs, _rel))

            if r[1] == r[2]
                r = r[1] < 0 ? (r[1], 0.0f0) : (0.0f0, r[2])
                colorscheme[] = ColorScheme([colorant"gray"])
            elseif r[1] < 0 && r[2] > 0
                r = (-maximum(abs.(r)), maximum(abs.(r)))
                colorscheme[] = ColorSchemes.coolwarm
            elseif r[1] ≥ 0
                r = (0.0f0, r[2])
                colorscheme.val = ColorSchemes.thermal
            end

            maxrange[] = r
            # adjust thumb position
            new_thumbs = get(thumb_pos_cache, thumb_pos_key(), r)
            thumb_l[], thumb_r[] = new_thumbs
        else
            error("More than one state for maxrange calculation...")
        end
        nothing
    end;

    onany(thumb_l, thumb_r; update=true) do _thumb_l, _thumb_r
        @debug "GP: $(type) color slider => app.gp.colorrange"
        # store the thumb position
        thumb_pos_cache[thumb_pos_key()] = (_thumb_l, _thumb_r)
        colorrange[] = (_thumb_l, _thumb_r)
        nothing
    end

    fig = with_theme(apptheme()) do
        fig = Figure(size=(400,40); figure_padding=10)
        Colorbar(fig[1,1];
            colormap=colorscheme,
            colorrange=colorrange,
            vertical=false, flipaxis=true)
        fig
    end

    cslider = ContinuousSlider(maxrange, thumb_l, thumb_r)

    Card(
        DOM.div(
            selector,
            DOM.div(fig; style=Styles("height" => "40px")),
            # fig,
            # RoundedLabel(@lift $maxrange[1]; style=Styles("text-align"=>"right")),
            cslider,
            # RoundedLabel(@lift $maxrange[2]; style=Styles("text-align"=>"left"));
            class="gpstate-control-card-content",
        );
        class
    )
end
function _maxrange(sol, idxs, rel)
    _isvalid(s) = SII.is_variable(sol, s) || SII.is_parameter(sol, s) || SII.is_observed(sol, s)
    mask = _isvalid.(idxs)

    u_for_t = sol(sol.t; idxs=idxs[mask]).u
    if rel
        for u in u_for_t
            u .-= u_for_t[1]
        end
    end
    extrema(Iterators.flatten(u_for_t))
end

function timeseries_cards(app, session)
    cards = OrderedDict{String,Hyperscript.Node{Hyperscript.HTMLSVG}}()
    container = Observable{Hyperscript.Node{Hyperscript.HTMLSVG}}()

    on(app.tsplots; update=true) do _tsplots
        @debug "TS: app.tsplots => update timeseries cards"
        newkeys = keys(_tsplots)
        knownkeys = keys(cards)

        for delkey in setdiff(knownkeys, newkeys)
            delete!(cards, delkey)
        end
        for newkey in setdiff(newkeys, knownkeys)
            cards[newkey] = timeseries_card(app, newkey, session)
            # cards[newkey] = DOM.div(scatter(rand(100)))
        end
        if keys(cards) != keys(_tsplots)
            @warn "The keys do not match: $(keys(cards)) vs $(keys(_tsplots))"
        end

        container[] = DOM.div(values(cards)...; class="timeseries-stack")

        nothing
    end

    on(app.active_tsplot; update=true) do active
        activesel = app.tsplots[][active].selcomp[]
        app.graphplot._selcomp[] = activesel
    end

    return container[]
end

function timeseries_card(app, key, session)
    tsplot = app.tsplots[][key]

    comp_options = Observable{Vector{OptionGroup{SymbolicIndex}}}()
    on(app.sol; update=true) do _sol
        @debug "TS: app.sol => comp_options"
        g = extract_nw(_sol).im.g
        vg = OptionGroup{SymbolicIndex}("Nodes", VIndex.(1:nv(g)))
        eg = OptionGroup{SymbolicIndex}("Edges", EIndex.(1:ne(g)))
        comp_options[] = [vg, eg]
        nothing
    end

    state_options = Observable{Vector{OptionGroup{Symbol}}}()
    onany(app.sol, tsplot.selcomp; update=true) do _sol, _sel
        @debug "TS: app.sol, tsplot.selcomp => state_options"
        _nw = extract_nw(_sol)
        state_options[] = gen_state_options(_nw, _sel)
        nothing
    end

    comp_sel = MultiSelect(comp_options, tsplot.selcomp;
        placeholder="Select components",
        multi=true,
        option_to_string=s -> _sidx_to_str(s, app),
        T=SymbolicIndex,
        id=gendomid("compsel"))
    # comp_sel_dom = Grid(DOM.span("Components"), comp_sel; columns = "70px 1fr", align_items = "center")
    state_sel = MultiSelect(state_options, tsplot.states;
        placeholder="Select states",
        multi=true,
        T=Symbol,
        id=gendomid("statesel"))

    reset_button = Bonito.Button("Reset Color", style=Styles("margin-left"=>"10px"))
    on(reset_button.value) do _
        empty!(color_cache)
        empty!(linestyle_cache)
        notify(tsplot.selcomp)
        notify(tsplot.states)
    end

    rel_toggle = ToggleSwitch(value=tsplot.rel, label="Rel to u0")

    comp_state_sel_dom = Grid(
        DOM.span("Components"), comp_sel, reset_button,
        DOM.span("States"), state_sel, rel_toggle;
        columns = "min-content auto min-content",
        align_items = "center",
        class = "comp-state-sel-grid"
    )

    # hl choice of elements in graphplot
    on(tsplot.selcomp; update=true) do _sel
        @debug "TS: comp selection => graphplot selection"
        app.graphplot._selcomp[] = _sel
        nothing
    end

    ####
    #### actual plot
    ####
    COLORS = Makie.wong_colors()
    LINESTYLES = [:solid, :dot, :dash, :dashdot, :dashdotdot]
    LINESTYLES = ["─", "⋯", "--", "-⋅-", "-⋅⋅"]
    color_cache = Dict{Union{EIndex{Int,Nothing},VIndex{Int,Nothing}}, Int}()
    linestyle_cache = Dict{Symbol,Int}()

    colorpairs = Observable{Vector{@NamedTuple{title::String,color::String}}}()
    lstylepairs = Observable{Vector{@NamedTuple{title::String,linestyle::String}}}()

    on(tsplot.selcomp; update=true) do _sel
        @debug "TS: comp selection => update color_cache"
        for unused in setdiff(keys(color_cache), _sel)
            delete!(color_cache, unused)
        end
        for new in setdiff(_sel, keys(color_cache))
            i = _smallest_free(color_cache)
            color_cache[new] = i
        end
        colorpairs[] = [(; title=_sidx_to_str(k, app),
                         color="#"*Colors.hex(getcycled(COLORS, v)))
                        for (k,v) in color_cache]
        nothing
    end
    on(tsplot.states; update=true) do _states
        @debug "TS: state selection => update linestyle_cache"
        for unused in setdiff(keys(linestyle_cache), _states)
            delete!(linestyle_cache, unused)
        end
        for new in setdiff(_states, keys(linestyle_cache))
            i = _smallest_free(linestyle_cache)
            linestyle_cache[new] = i
        end
        lstylepairs[] = [(; title=repr(s), linestyle=getcycled(LINESTYLES, i))
                         for (s,i) in linestyle_cache]
        nothing
    end

    comp_ms_id = Bonito.JSString(comp_sel.id)
    state_ms_id = Bonito.JSString(state_sel.id)

    js = js"""
    function colorListItems(items) {
        let styleTag = document.getElementById("dynamic-style-$(comp_ms_id)");

        // Create the <style> tag if it doesn't exist
        if (!styleTag) {
            styleTag = document.createElement("style");
            styleTag.id = "dynamic-style-$(comp_ms_id)";
            document.head.appendChild(styleTag);
        }

        // Clear previous styles
        styleTag.innerHTML = "";

        // Generate new styles
        let styleContent = "";
        items.forEach(({ title, color }) => {
            styleContent += `#$(comp_ms_id) +span li[title='${title}']::after {
                content: 'xx';
                display: inline-block;
                padding: 0px 4px;
                background-color: ${color} !important;
                color: ${color} !important;
                border-left: 1px solid #aaa;
                cursor: default;
            }\n`;
        });

        // Insert new styles
        styleTag.innerHTML = styleContent;
    }

    colorListItems($(colorpairs).value);
    $(colorpairs).on((c) => {
        colorListItems(c);
    });

    function linestyleListItems(items) {
        let styleTag = document.getElementById("dynamic-style-$(state_ms_id)");

        // Create the <style> tag if it doesn't exist
        if (!styleTag) {
            styleTag = document.createElement("style");
            styleTag.id = "dynamic-style-$(state_ms_id)";
            document.head.appendChild(styleTag);
        }

        // Clear previous styles
        styleTag.innerHTML = "";

        // Generate new styles
        let styleContent = "";

        if(items.length > 1) {
            items.forEach(({ title, linestyle }) => {
                styleContent += `#$(state_ms_id) +span li[title='${title}']::after {
                    content: '${linestyle}';
                    display: inline-block;
                    padding: 0px 4px;
                    color: inherit;
                    border-left: 1px solid #aaa;
                    font-size: smaller;
                    cursor: default;
                }\n`;
            });

            // Insert new styles
            styleTag.innerHTML = styleContent;
        }
    }
    linestyleListItems($(lstylepairs).value);
    $(lstylepairs).on((c) => {
        linestyleListItems(c);
    });
    """
    evaljs(session, js)

    fig, ax = with_theme(apptheme()) do
        fig = Figure()
        ax = Axis(fig[1, 1])
        fig, ax
    end

    # set axis limits according to time range
    onany(app.tmin, app.tmax, update=true) do _tmin, _tmax
        @debug "TS: tim/tmax => update ax limits"
        xlims!(ax, (_tmin, _tmax))
        nothing
    end

    ts = Observable(range(app.sol[].t[begin], app.sol[].t[end], length=1000))
    # last update of ts range
    lastupdate = Ref(time())
    on(ax.finallimits) do lims
        lastupdate[] = time()
        nothing
    end
    # every 0.5 seconds trigger resampling
    # timer = Timer(.5; interval=.5) do _
    #     if time() > lastupdate[] + 0.5
    #         lims = ax.finallimits[]
    #         tmin = max(app.sol[].t[begin], lims.origin[1])
    #         tmax = min(app.sol[].t[end], tmin + lims.widths[1])
    #         ts[] = range(tmin, tmax, length=1000)
    #         lastupdate[] = Inf
    #     end
    # end
    # on(session.on_close) do _
    #     @info "Session closed, time to clean up"
    #     # close(timer)
    # end

    # collect all the states wie might want to plot
    valid_idxs = Observable(
        Union{VIndex{Int,Symbol},EIndex{Int,Symbol}}[]
    )
    onany(tsplot.selcomp, tsplot.states; update=true) do _selcomp, _states
        @debug "TS: sel comp/states => update valid_idxs"
        isvalid(s) = SII.is_variable(app.sol[], s) || SII.is_parameter(app.sol[], s) || SII.is_observed(app.sol[], s)
        empty!(valid_idxs[])
        for c in _selcomp, s in _states
            idx = c isa VIndex ? VIndex(c.compidx, s) : EIndex(c.compidx, s)
            isvalid(idx) && push!(valid_idxs[], idx)
        end
        notify(valid_idxs)
    end

    # extract the data
    data = Observable{Vector{Vector{Float32}}}(Vector{Float32}[])
    onany(ts, valid_idxs, tsplot.rel, app.sol; update=true) do _ts, _valid_idxs, _rel, _sol
        @debug "TS: t, valid_idx, rel, sol => update data"
        _dat = _sol(_ts, idxs=_valid_idxs)
        if _rel
            u0 = _sol(_sol.t[begin], idxs=_valid_idxs)
            for row in _dat
                row .-= u0
            end
        end
        resize!(data[], length(_valid_idxs))
        for i in 1:length(_valid_idxs)
            data.val[i] = _dat[i,:]
        end
        notify(data)
    end

    # plot the thing
    on(data; update=true) do _dat
        @async begin
            try
                empty!(ax)
                vlines!(ax, app.t; color=:black)
                for (idx, y) in zip(valid_idxs[], data[])
                    color = begin
                        key = idx isa VIndex ? VIndex(idx.compidx) : EIndex(idx.compidx)
                        Cycled(color_cache[key])
                    end
                    linestyle = Cycled(linestyle_cache[idx.subidx])
                    lines!(ax, ts[], y; label=string(idx), color, linestyle)
                end
            catch e
                @error "Plotting failed" e
            end
        end
        nothing
    end

    ####
    #### Click interaction to set time
    ####
    set_time_interaction = (event::MouseEvent, axis) -> begin
        if event.type === MouseEventTypes.leftclick
            @debug "TS: click on axis => update t"
            pos = mouseposition(axis.scene)[1]
            app.t[] = pos
            return Consume(true)
        end
        return Consume(false)
    end
    register_interaction!(set_time_interaction, ax, :set_time)

    cardclass = "timeseries-card"
    if key == app.active_tsplot[]
        cardclass *= " active-tseries"
    end

    card = Card(
        DOM.div(
            comp_state_sel_dom,
            DOM.div(fig; class="timeseries-axis-container");
            class="timeseries-card-container"
        );
        class=cardclass
    )

    # on click set active-tseries class
    click = js"""
    (card) => {
        card.addEventListener("click", function(event) {
            $(app.active_tsplot).notify($(key));

            document.querySelectorAll(".timeseries-card").forEach(element => {
                element.classList.remove("active-tseries");
            });

            // Add "active-tseries" to the given target element
            card.classList.add("active-tseries");
        }, { capture: true });
    }
    """
    Bonito.onload(session, card, click)

    return card
end

function _sidx_to_str(s, app)
    if s isa VIndex
        "v$(s.compidx)"
    else
        edge = extract_nw(app.sol[]).im.edgevec[s.compidx]
        src, dst = edge.src, edge.dst
        "e$(s.compidx): $src→$dst"
    end
end
function _smallest_free(d::Dict)
    vals = values(d)
    i = 1
    while i ∈ vals
        i += 1
    end
    return i
end

end # module NetworkDynamicsInspector
