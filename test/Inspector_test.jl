using NetworkDynamicsInspector
using NetworkDynamics
using Bonito
using WGLMakie
using GraphMakie
using Graphs: SimpleGraph
using OrdinaryDiffEqTsit5

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
    trange_full = Observable{Tuple{Float64, Float64}}((0.0, 10.0)),
    tmin = Observable{Float64}(1.0),
    tmax = Observable{Float64}(9.0),
    t = Observable{Float64}(0.0)
)

App() do session
    tw_slider = ContinuousSlider(app.trange_full, app.tmin, app.tmax)
    twindow = @lift ($(app.tmin), $(app.tmax))
    # if window changes, keep t[] inside
    onany(app.tmin, app.tmax) do tmin, tmax
        _t = clamp(app.t[], tmin, tmax)
        if _t != app.t[]
            app.t[] = _t
        end
    end
    t_slider = ContinuousSlider(twindow, app.t)

    return Card(Grid(
        DOM.div(), t_slider, Bonito.Label(t_slider.value_r),
        Bonito.Label(tw_slider.value_l), tw_slider, Bonito.Label(tw_slider.value_r);
        columns="10% 80% 10%",
        justify_content="begin",
        align_items="center",
    ); width="500px",)
end

app.tmin[] = NaN
app.tmin[] = 1.0

app.tmin[]
app.tmax[]
