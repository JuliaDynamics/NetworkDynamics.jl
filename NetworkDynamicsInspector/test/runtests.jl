using Test
using NetworkDynamics
using NetworkDynamicsInspector
using NetworkDynamicsInspector: NetworkDynamicsInspector as NDI
using ExplicitImports
using Aqua
using Electron
using Bonito
using OrdinaryDiffEqTsit5
using NetworkDynamics.Graphs
using SafeTestsets

include(joinpath(pkgdir(NetworkDynamics), "test", "ComponentLibrary.jl"))

function get_sol(;limit=1.0)
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
    s0.p.e[:, :limit] .= limit

    prob = ODEProblem(nw, uflat(s0), (0,6), copy(pflat(s0)), callback=nwcb)
    sol = solve(prob, Tsit5())
end

@testset "NetworkDynamicsInspector.jl Tests" begin
    @testset "Package Quality Tests" begin
        @test check_no_implicit_imports(NetworkDynamicsInspector) === nothing
        @test check_no_stale_explicit_imports(NetworkDynamicsInspector) === nothing

        Aqua.test_all(NetworkDynamicsInspector;
            ambiguities=false,
            stale_deps=true,
            persistent_tasks=false,
            piracies = (;treat_as_own = [NetworkDynamics.extract_nw])
        )

        @test_broken isempty(Docs.undocumented_names(NetworkDynamicsInspector))
    end

    @testset "Test simple app creation and add + remove of timeseries" begin
        sol = get_sol()
        inspect(sol; restart=true, reset=true, display=ElectronDisp())

        delete!(NDI.APPSTATE[].tsplots[], "ts-1")
        notify(NDI.APPSTATE[].tsplots)

        NDI.APPSTATE[].tsplots[]["ts-2"] = NDI.TimeseriesPlot(selcomp=[VIndex(3), VIndex(2), VIndex(1)], states=[:ω])
        notify(NDI.APPSTATE[].tsplots)

        NDI.APPSTATE[].tsplots[]["ts-1"] = NDI.TimeseriesPlot(selcomp=[EIndex(3), EIndex(2), EIndex(1)], states=[:P])
        notify(NDI.APPSTATE[].tsplots)

        NDI.APPSTATE[].tsplots[]["ts-2"] = NDI.TimeseriesPlot(selcomp=[VIndex(2)], states=[:θ])
        notify(NDI.APPSTATE[].tsplots)

        empty!(NDI.APPSTATE[].tsplots[])
        notify(NDI.APPSTATE[].tsplots)
    end

    @testset "track changes in tsplots" begin
        sol = get_sol()
        inspect(sol; display=BrowserDisp())
        NDI.APPSTATE[].tsplots[]["ts-1"] = NDI.TimeseriesPlot(selcomp=[EIndex(3), EIndex(2), EIndex(1)], states=[:P])
        notify(NDI.APPSTATE[].tsplots) # should update
    end

    @testset "Test different display types" begin
        sol1 = get_sol()
        sol2 = get_sol(; limit = 1.2)

        inspect(sol1; display=BrowserDisp(), restart=true, reset=true)
        # change observable
        empty!(NDI.APPSTATE[].tsplots[]); notify(NDI.APPSTATE[].tsplots)
        NDI.APPSTATE[].tsplots[]["ts-1"] = NDI.TimeseriesPlot(selcomp=[EIndex(i) for i in 1:7], states=[:P])
        # notify(NDI.APPSTATE[].tsplots) # does not work wihtout display
        inspect(sol2) # this should update the plot
        inspect(sol2, reset=true)

        # switch to electron
        inspect(sol1; display=ElectronDisp())
        @test NDI.BROWSER_STATE[:handler] == nothing
        empty!(NDI.APPSTATE[].tsplots[]); notify(NDI.APPSTATE[].tsplots)
        NDI.APPSTATE[].tsplots[]["ts-1"] = NDI.TimeseriesPlot(selcomp=[EIndex(i) for i in 1:7], states=[:P])
        notify(NDI.APPSTATE[].tsplots)
        inspect(sol2) # this should update the plot

        # switch to server
        inspect(sol1; display=ServerDisp(), restart=true)
        @test NDI.BROWSER_STATE[:handler] == nothing
        empty!(NDI.APPSTATE[].tsplots[]); notify(NDI.APPSTATE[].tsplots)
        NDI.APPSTATE[].tsplots[]["ts-1"] = NDI.TimeseriesPlot(selcomp=[EIndex(i) for i in 1:7], states=[:P])
        # notify(NDI.APPSTATE[].tsplots) # does not work without session open
        inspect(sol2) # this should update the plot
    end

    @testset "Widget Tests" begin
        @safetestset "Multiselect Tests" begin include("multiselect_test.jl") end
        @safetestset "Slider Tests" begin include("slider_test.jl") end
    end
end
