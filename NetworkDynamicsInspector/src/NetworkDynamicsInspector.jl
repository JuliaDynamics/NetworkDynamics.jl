module NetworkDynamicsInspector

using Bonito
using NetworkDynamics
using Bonito.Observables
using Bonito: Grid, @js_str, onload, jsrender
using WGLMakie.Makie
using WGLMakie.Makie.Colors
using WGLMakie.Makie.ColorSchemes
using NetworkDynamics: extract_nw
using NetworkDynamics: SII
using Graphs: nv, ne
using GraphMakie
using GraphMakie.NetworkLayout

export ContinuousSlider, RoundedLabel
include("widgets.jl")

export wrap_assets
function wrap_assets(appdom)
    jquery = Asset("https://cdn.jsdelivr.net/npm/jquery@3.7.1/dist/jquery.min.js")
    select2_css = Asset("https://cdn.jsdelivr.net/npm/select2@4.1.0-rc.0/dist/css/select2.min.css")
    select2_js = Asset("https://cdn.jsdelivr.net/npm/select2@4.1.0-rc.0/dist/js/select2.min.js")
    css = Asset(joinpath(pkgdir(NetworkDynamicsInspector), "assets", "app.css"))

    DOM.body(
        jquery,
        select2_css,
        select2_js,
        css,
        appdom;
    )
end

function graphplot_card(app; kwargs...)
    nw = map!(extract_nw, Observable{Network}(), app.sol)
    NV = nv(nw[])
    NE = ne(nw[])

    node_marker = map(nw) do _nw
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
        idxs = VIndex(1:NV, _state)
        _gracefully_extract_states!(node_state[], _sol, _t, idxs, _rel)
        notify(node_state)
    end;

    edge_state = Observable(Vector{Float32}(undef, NE))
    onany(app.sol, app.t, app.graphplot.estate, app.graphplot.nstate_rel) do _sol, _t, _state, _rel
        idxs = EIndex(1:NE, _state)
        _gracefully_extract_states!(edge_state[], _sol, _t, idxs, _rel)
        notify(edge_state)
    end;

    node_color = Observable(Vector{RGB{Float64}}(undef, NV))
    onany(node_state, app.graphplot.ncolorrange, app.graphplot.ncolorscheme) do statevec, range, scheme
        for i in 1:NV
            node_color[][i] = isnan(statevec[i]) ? RGB(0,0,0) : get(scheme, statevec[i], range)
        end
        notify(node_color)
    end

    edge_color = Observable(Vector{RGB{Float64}}(undef, NE))
    onany(edge_state, app.graphplot.ecolorrange, app.graphplot.ecolorscheme) do statevec, range, scheme
        for i in 1:NV
            edge_color[][i] = isnan(statevec[i]) ? RGB(0,0,0) : get(scheme, statevec[i], range)
        end
        notify(edge_color)
    end

    notify(app.t) # trigger updates

    SMALL = 30
    BIG = 50
    node_size = Observable(fill(SMALL, NV))
    onany(app.sel_nodes; update=true) do selected
        fill!(node_size[], SMALL)
        for sel in selected
            node_size[][sel] = BIG
        end
        notify(node_size)
    end

    THIN = 3
    THICK = 6
    edge_width = Observable(fill(THIN, NE))
    onany(app.sel_edges; update=true) do selected
        fill!(edge_width[], THIN)
        for sel in selected
            edge_width[][sel] = THICK
        end
        notify(edge_width)
    end

    g = @lift $(nw).im.g

    layout = map(nw) do _nw
        pin = Dict{Int,Point2f}()
        for i in 1:NV
            vm = _nw[VIndex(i)]
            if has_position(vm)
                pin[i] = get_position(vm)
            end
        end
        Stress(; pin)
    end

    fig = Figure(; figure_padding=0)
    ax = Axis(fig[1,1])
    graphplot!(ax, g; layout, node_marker, node_color, node_size, edge_color, edge_width,
        node_attr=(;colorrange=app.graphplot.ncolorrange, colormap=app.graphplot.ncolorscheme),
        edge_attr=(;colorrange=app.graphplot.ncolorrange, colormap=app.graphplot.ncolorscheme))

    hidespines!(ax)
    hidedecorations!(ax)
    on(ax.scene.viewport) do lims
        adapt_xy_scaling!(ax)
    end
    Card(fig; kwargs...)
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
function adapt_xy_scaling!(ax)
    vp = ax.scene.viewport[]
    fl = ax.finallimits[]
    vp_aspect = reduce(/, vp.widths)
    fl_aspect = reduce(/, fl.widths)
    ratio = vp_aspect/fl_aspect

    if ratio > 1 # need to increase x limits
        low = ax.finallimits[].origin[1]
        hi = low + ax.finallimits[].widths[1]
        xlims!(ax, ratio*low, ratio*hi)
    else # need to increast y limits
        low = ax.finallimits[].origin[2]
        hi = low + ax.finallimits[].widths[2]
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
            app.t[] = _t
        end
    end
    t_slider = ContinuousSlider(twindow, app.t)
    Card(
        Grid(
            DOM.span("Time"), t_slider, RoundedLabel(t_slider.value_r),
            RoundedLabel(tw_slider.value_l; style=Styles("text-align"=>"right")),
            tw_slider,
            RoundedLabel(tw_slider.value_r; style=Styles("text-align"=>"left"));
            columns="10% 80% 10%",
            justify_content="begin",
            align_items="center",
        );
    );
end

function nodestate_control_card(app)
    NV = nv(extract_nw(app.sol[]))
    maxrange = Observable{Tuple{Float32,Float32}}()
    onany(app.sol, app.graphplot.nstate, app.graphplot.nstate_rel; update=true) do _sol, _state, _rel
        idxs = VIndex(1:NV, _state)
        r = _maxrange(_sol, idxs, _rel)
        if r[1] < 0 && r[2] > 0
            r = (-maximum(abs.(r)), maximum(abs.(r)))
        end
        maxrange[] = r
    end;
    clow = Observable{Float32}(maxrange[][1])
    chi = Observable{Float32}(maxrange[][2])

    onany(clow, chi; update=true) do _clow, _chi
        app.graphplot.ncolorrange[] = (_clow, _chi)
    end

    fig = Figure(; figure_padding=0)
    Colorbar(fig[1,1];
        colormap=app.graphplot.ncolorscheme,
        colorrange=app.graphplot.ncolorrange,
        vertical=false, flipaxis=true)

    cslider = ContinuousSlider(maxrange, clow, chi)

    Card(
        Grid(
            DOM.div(fig; style=Styles("height" => "40px")),
            # RoundedLabel(@lift $maxrange[1]; style=Styles("text-align"=>"right")),
            cslider,
            # RoundedLabel(@lift $maxrange[2]; style=Styles("text-align"=>"left"));
            columns="100%",
        )
    )
end
function _maxrange(sol, idxs, rel)
    isvalid(s) = SII.is_variable(sol, s) || SII.is_parameter(sol, s) || SII.is_observed(sol, s)
    mask = isvalid.(idxs)

    u_for_t = sol(sol.t; idxs=idxs[mask]).u
    if rel
        for u in u_for_t
            u .-= u_for_t[1]
        end
    end
    extrema(Iterators.flatten(u_for_t))
end

function clear_obs!(nt::NamedTuple)
    for v in values(nt)
        clear_obs!(v)
    end
end
clear_obs!(obs::Observable) = empty!(obs.listeners)
clear_obs!(x) = x

# TODO: move to grpahmakie
GraphMakie._dimensionality(obs::Observable, g) = GraphMakie._dimensionality(obs[], g)


SERVER = Ref{Any}(nothing)

export serve_app
function serve_app(newapp)
    if !isnothing(SERVER[]) && Bonito.HTTPServer.isrunning(SERVER[])
        @info "Stop running server..."
        close(SERVER[])
    end
    SERVER[] = Bonito.Server(newapp, "0.0.0.0", 8080)
end

# APP = Ref{Any}(nothing)
# function serve_app(newapp)
#     if isnothing(SERVER[]) || !Bonito.HTTPServer.isrunning(SERVER[])
#         @info "Start new SErver"
#         SERVER[] = Bonito.Server(newapp, "0.0.0.0", 8080)
#         APP[] = newapp
#     else
#         oldapp = APP[]
#         APP[] = newapp
#         @info "Update App"
#         Bonito.update_app!(oldapp, newapp)
#     end
#     nothing
# end

end # module NetworkDynamicsInspector
