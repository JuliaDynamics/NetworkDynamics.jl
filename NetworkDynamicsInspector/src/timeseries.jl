function timeseries_col(app, session)
    col = DOM.div(
        timeslider_card(app),
        timeseries_cards(app, session),
        add_timeseries_button(app),
        class="timeseries-col"
    )

    # add js to update the active ts plot
    onload_js = js"""
    (tscol) => {
        tscol.addEventListener("click", function(event) {
            // Find the closest .timeseries-card ancestor of the clicked element
            const clickedCard = event.target.closest('.timeseries-card');

            // If a .timeseries-card was found, handle the click
            if (clickedCard) {
                // console.log("click on timeseries card", clickedCard.id);

                // Notify the app with the id of the clicked card
                $(app.active_tsplot).notify(clickedCard.id);

                // Remove "active-tseries" class from all .timeseries-card elements
                document.querySelectorAll(".timeseries-card").forEach(element => {
                    element.classList.remove("active-tseries");
                });

                // Add "active-tseries" class to the clicked card
                clickedCard.classList.add("active-tseries");
            }
        }, { capture: true });

        function triggerWindowResize() {
            window.dispatchEvent(new Event('resize'));
        }
        const triggerWindowResize_throttled = Bonito.throttle_function(triggerWindowResize, 100);
        const resizeObserver = new ResizeObserver(triggerWindowResize_throttled);
        resizeObserver.observe(tscol);
    }
    """
    Bonito.onload(session, col, onload_js)

    return col
end

function timeseries_cards(app, session)
    cards = OrderedDict{String,Bonito.Hyperscript.Node{Bonito.Hyperscript.HTMLSVG}}()
    observables = OrderedDict{String, Vector{Observables.ObserverFunction}}()
    container = Observable{Bonito.Hyperscript.Node{Bonito.Hyperscript.HTMLSVG}}()

    on(app.tsplots; update=true) do _tsplots
        @debug "TS: app.tsplots => update timeseries cards"
        newkeys = collect(keys(_tsplots)) # collect to preserv order on setdiff
        knownkeys = keys(cards)

        for delkey in setdiff(knownkeys, newkeys)
            Observables.off.(observables[delkey]) # deactivate allobservables from card
            delete!(observables, delkey)
            delete!(cards, delkey)
        end
        for newkey in setdiff(newkeys, knownkeys)
            card, obsf = timeseries_card(app, newkey, session)
            cards[newkey] = card
            observables[newkey] = obsf
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

    return container
end

function add_timeseries_button(app)
    button = Bonito.Button("Add Timeseries", class="add-ts-button")
    on(button.value) do _
        newkey = free_ts_key()
        @debug "TS: Add Timeseries button clicked. Add $newkey"
        app.tsplots[][newkey] = TimeseriesPlot()
        notify(app.tsplots)
    end
    return button
end

function timeseries_card(app, key, session)
    # store all observer functions to be able to remove them later
    obsf = Observables.ObserverFunction[]
    track_obsf = function (f)
        if f isa AbstractVector
            append!(obsf, f)
        else
            push!(obsf, f)
        end
    end

    tsplot = app.tsplots[][key]

    comp_options = Observable{Vector{OptionGroup{SymbolicCompIndex}}}()
    on(app.sol; update=true) do _sol
        @debug "TS: app.sol => comp_options"
        g = extract_nw(_sol).im.g
        vg = OptionGroup{SymbolicCompIndex}("Nodes", VIndex.(1:nv(g)))
        eg = OptionGroup{SymbolicCompIndex}("Edges", EIndex.(1:ne(g)))
        comp_options[] = [vg, eg]
        nothing
    end |> track_obsf

    state_options = Observable{Vector{OptionGroup{Symbol}}}()
    onany(app.sol, tsplot.selcomp; update=true) do _sol, _sel
        @debug "TS: app.sol, tsplot.selcomp => state_options"
        _nw = extract_nw(_sol)
        state_options[] = gen_state_options(_nw, _sel)
        nothing
    end |> track_obsf

    comp_sel = TomSelect(comp_options, tsplot.selcomp;
        placeholder="Select components",
        multi=true,
        option_to_string=s -> sidx_to_str(s, app),
        T=SymbolicCompIndex,
        id=gendomid("compsel"))
    # comp_sel_dom = Grid(DOM.span("Components"), comp_sel; columns = "70px 1fr", align_items = "center")
    state_sel = TomSelect(state_options, tsplot.states;
        placeholder="Select states",
        multi=true,
        T=Symbol,
        id=gendomid("statesel"))

    reset_style_button = Bonito.Button("Reset Styles", style=Styles("margin-left"=>"10px"))
    on(reset_style_button.value) do _
        empty!(color_cache)
        empty!(linestyle_cache)
        notify(tsplot.selcomp)
        notify(tsplot.states)
    end |> track_obsf

    reset_axis_button = Bonito.Button("Reset Axis", style=Styles("margin-left"=>"10px"))

    rel_toggle = ToggleSwitch(value=tsplot.rel, label="Rel to u0")

    comp_state_sel_dom = Grid(
        DOM.span("Components"), comp_sel,
        DOM.div(reset_style_button, reset_axis_button, rel_toggle; style=Styles("grid-row"=>"1/3", "grid-column"=>"3")),
        DOM.span("States"), state_sel;
        columns = "min-content auto min-content",
        align_items = "center",
        class = "comp-state-sel-grid"
    )

    # hl choice of elements in graphplot
    on(tsplot.selcomp; update=true) do _sel
        @debug "TS: comp selection => graphplot selection"
        app.graphplot._selcomp[] = _sel
        nothing
    end |> track_obsf

    ####
    #### actual plot
    ####
    COLORS = Makie.wong_colors()
    LINESTYLES = [:solid, :dot, :dash, :dashdot, :dashdotdot]
    LINESTYLES_STR = ["─", "⋯", "--", "-⋅-", "-⋅⋅"]
    color_cache = Dict{SymbolicCompIndex, Int}()
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
        colorpairs[] = [(; title=sidx_to_str(k, app),
                         color="#"*Colors.hex(getcycled(COLORS, v)))
                        for (k,v) in color_cache]
        nothing
    end |> track_obsf
    on(tsplot.states; update=true) do _states
        @debug "TS: state selection => update linestyle_cache"
        for unused in setdiff(keys(linestyle_cache), _states)
            delete!(linestyle_cache, unused)
        end
        for new in setdiff(_states, keys(linestyle_cache))
            i = _smallest_free(linestyle_cache)
            linestyle_cache[new] = i
        end
        lstylepairs[] = [(; title=repr(s), linestyle=getcycled(LINESTYLES_STR, i))
                         for (s,i) in linestyle_cache]
        nothing
    end |> track_obsf

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
            styleContent += `#$(comp_ms_id) +div div[data-value='${title}']::before {
                content: '';
                background-color: ${color};
                width: 12px;
                height: 12px;
                border-radius: 50%;
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
                styleContent += `#$(state_ms_id) +div div[data-value='${title}']::before {
                    content: '${linestyle}';
                    color: inherit;
                    border-right: 1px solid #d0d0d0;
                    padding-right: 5px;
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
    Bonito.evaljs(session, js)

    fig, ax = with_theme(apptheme()) do
        fig = Figure(size=(100, 400))
        ax = Axis(fig[1, 1])
        fig, ax
    end

    # set axis limits according to time range
    onany(app.tmin, app.tmax, update=true) do _tmin, _tmax
        @debug "TS: tim/tmax => update ax limits"
        xlims!(ax, (_tmin, _tmax))
        nothing
    end |> track_obsf

    ts = Observable(collect(range(app.sol[].t[begin], app.sol[].t[end], length=1000)))
    refined_xlims = Ref((NaN, NaN))
    onany_delayed(ax.finallimits; delay=0.5) do axlims
        sollims = (app.sol[].t[begin], app.sol[].t[end])
        xlims = (axlims.origin[1], axlims.origin[1] + axlims.widths[1])
        if xlims != refined_xlims[]
            refined_xlims[] = xlims
            _refine_time_limits!(ts[], sollims, xlims)
            notify(ts)
        end
        nothing
    end |> track_obsf

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
    end |> track_obsf

    # extract the data
    data = Observable{Vector{Vector{Float32}}}(Vector{Float32}[])
    datahash = Ref{UInt64}()
    onany(ts, valid_idxs, tsplot.rel, app.sol; update=true) do _ts, _valid_idxs, _rel, _sol
        # only replot if the data has actually changed
        newhash = hash((_ts, _valid_idxs, _rel, _sol, linestyle_cache, color_cache))
        newhash == datahash[] && return
        datahash[] = newhash

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
    end |> track_obsf

    replot = Observable{Nothing}(nothing)

    # store the idxs for which the autolmits where last set
    last_autolimits = Ref((eltype(valid_idxs)(), tsplot.rel[]))
    # plot the thing
    onany(data, replot; update=true) do _dat, _
        @async begin
            try
                empty!(ax)
                vlines!(ax.scene, app.t; color=:black)
                for (idx, y) in zip(valid_idxs[], data[])
                    color = begin
                        key = idx isa VIndex ? VIndex(idx.compidx) : EIndex(idx.compidx)
                        getcycled(COLORS, color_cache[key])
                    end
                    linestyle = getcycled(LINESTYLES, linestyle_cache[idx.subidx])
                    lines!(ax.scene, ts[], y; label=string(idx), color, linestyle)
                    # scatterlines!(ax, ts[], y; label=string(idx), color, linestyle)
                end
                # if last_autolimits[][1] != valid_idxs[] || last_autolimits[][2] != tsplot.rel[]
                if last_autolimits[] != (valid_idxs[], tsplot.rel[])
                    autolimits!(ax)
                    xlims!(ax, (app.tmin[], app.tmax[]))
                    last_autolimits[] = (copy(valid_idxs[]), tsplot.rel[])
                end
            catch e
                @error "Plotting failed" e
            end
        end
        nothing
    end |> track_obsf

    on(reset_axis_button.value) do _
        autolimits!(ax)
        xlims!(ax, (app.tmin[], app.tmax[]))
    end |> track_obsf

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

    cardclass = "bonito-card timeseries-card"
    if key == app.active_tsplot[]
        cardclass *= " active-tseries"
    end

    help = HoverHelp(html"""
    Select the components and states to plot.
    <ul>
    <li><strong>Click</strong> on the plot to set the time.</li>
    <li><strong>Click and Drag</strong> to zoom in.</li>
    <li><strong>Ctrl + Click</strong> to reset zoom.</li>
    <li>Use timeslider to change the x-axis limits.</li>
    <li>Toggle "Rel to u0" to plot the difference to the initial value.</li>
    <li>Click "Reset Styles" to pick colors and linestyles in order of apperance.</li>
    </ul>
    """)

    card = Card(
        [DOM.div(
            comp_state_sel_dom,
            DOM.div(fig; class="timeseries-axis-container"),
            closebutton(app, key);
            class="timeseries-card-container"
        ), help];
        class=cardclass,
        id=key
    )

    return card, obsf
end
function closebutton(app, key)
    button = Bonito.Button("×", class="close-button")
    on(button.value) do _
        @debug "TS: Close button clicked. Remove $key"
        delete!(app.tsplots[], key)
        notify(app.tsplots)
    end
    button
end

function _refine_time_limits!(ts, sollims, axlims)
    FOCUS = 900
    UNFOCUS = 100
    if length(ts) !== FOCUS + UNFOCUS
        @warn "Resize ts, lenght $(length(ts)) != $(FOCUS + UNFOCUS)"
        resize!(ts, FOCUS + UNFOCUS)
    end
    # never plot outside of solution
    axlims = (max(axlims[1], sollims[1]), min(axlims[2], sollims[2]))

    axrange = range(axlims[1], axlims[2], length=FOCUS)
    outer_low =  max(0, axlims[1] - sollims[1])
    outer_high = max(0, sollims[2] - axlims[2])

    if outer_low == 0
        high_range = range(sollims[2], axlims[2], length=UNFOCUS+1)[2:end]
        ts[1:FOCUS] .= axrange
        ts[FOCUS+1:end] .= high_range
    elseif outer_high == 0
        low_range = range(axlims[1], sollims[1], length=UNFOCUS+1)[1:end-1]
        ts[1:UNFOCUS] .= low_range
        ts[UNFOCUS+1:end] .= axrange
    else
        N_LOW = round(Int, outer_low/(outer_low + outer_high) *UNFOCUS)
        N_HIGH = UNFOCUS - N_LOW
        low_range = range(sollims[1], axlims[1], length=N_LOW+1)[1:end-1]
        high_range = range(axlims[2], sollims[2], length=N_HIGH+1)[2:end]
        ts[1:N_LOW] .= low_range
        ts[N_LOW+1:N_LOW+FOCUS] .= axrange
        ts[N_LOW+FOCUS+1:end] .= high_range
    end
    # ts .= range(axlims[1], axlims[2], length=1000)
    nothing
end

function _smallest_free(d::Dict)
    vals = values(d)
    i = 1
    while i ∈ vals
        i += 1
    end
    return i
end
