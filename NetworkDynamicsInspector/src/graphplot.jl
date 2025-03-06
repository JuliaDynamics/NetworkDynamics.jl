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
    THICK = 10
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
        _adapt_xy_scaling!(xratio, yratio, ax)
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
    Bonito.evaljs(session, js)

    nhh = NodeHoverHandler() do hstate, idx, event, axis
        tooltip_text[] = hstate ? sidx_to_str(VIndex(idx), app) : ""
        app.graphplot._hoverel[] = hstate ? VIndex(idx) : nothing
    end
    ehh = EdgeHoverHandler() do hstate, idx, event, axis
        tooltip_text[] = hstate ? sidx_to_str(EIndex(idx), app) : ""
        app.graphplot._hoverel[] = hstate ? EIndex(idx) : nothing
    end

    # interactions
    clickaction = (i, ax, type) -> begin
        shift_pressed = Makie.Keyboard.left_shift ∈ events(ax).keyboardstate ||
                        Makie.Keyboard.right_shift ∈ events(ax).keyboardstate

        idx = type == :vertex ? VIndex(i) : EIndex(i)
        if !shift_pressed # only update info window when click + shift
            selcomp = app.tsplots[][app.active_tsplot[]].selcomp
            if idx ∉ selcomp[]
                push!(selcomp[], idx)
            else
                filter!(x -> x != idx, selcomp[])
            end
            notify(selcomp)
        end
        app.graphplot._lastclickel[] = idx
    end
    nch = NodeClickHandler((i, _, ax) -> clickaction(i, ax, :vertex))
    ech = EdgeClickHandler((i, _, ax) -> clickaction(i, ax, :edge))

    register_interaction!(ax, :nodeclick, nch)
    register_interaction!(ax, :nodehover, nhh)
    register_interaction!(ax, :edgeclick, ech)
    register_interaction!(ax, :edgehover, ehh)

    Card(fig; class="bonito-card graphplot-card")
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
function _adapt_xy_scaling!(xratio, yratio, ax)
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

function gpstate_control_card(app, type)
    NEL = type == :vertex ? nv : ne
    VEIndex = type == :vertex ? VIndex : EIndex
    label = type == :vertex ? "Node state:" : "Edge state:"
    stateobs = type == :vertex ? app.graphplot.nstate : app.graphplot.estate
    stateobs_rel = type ==:vertex ? app.graphplot.nstate_rel : app.graphplot.estate_rel
    colorrange = type == :vertex ? app.graphplot.ncolorrange : app.graphplot.ecolorrange
    colorscheme = type == :vertex ? app.graphplot.ncolorscheme : app.graphplot.ecolorscheme
    class = type == :vertex ? "gpstate-control-card vertex" : "gpstate-control-card edge"
    class *= " resize-with-gp bonito-card"

    ####
    #### State selection
    ####
    options = Observable{Vector{OptionGroup{Symbol}}}()
    on(app.sol; update=true) do _sol
        _nw = extract_nw(_sol)
        idxs = VEIndex.(1:NEL(_nw))
        options[] = gen_state_options(_nw, idxs)
        nothing
    end
    multisel = TomSelect(options, stateobs; placeholder="Select state for coloring", multi=false, T=Symbol)
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
            idxs = VEIndex.(1:NEL(extract_nw(_sol)), only(_state))
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


