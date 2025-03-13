# Interactive Solution Inspection

An interactive solution inspection tool based on [WGLMakie](https://makie.org/website/) and [Bonito](https://github.com/SimonDanisch/Bonito.jl) is provided through the helper package `NetworkDynamicsInspector`.

First, we need to define the system we want to inspect.

!!! details "Define some network, simulate it and get a solution object"
    ```@example ndi
    using NetworkDynamics
    using NetworkDynamicsInspector
    using NetworkDynamicsInspector: NetworkDynamicsInspector as NDI #hide
    using Electron # hide
    using OrdinaryDiffEqTsit5
    using Graphs

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

    sol = get_sol()
    nothing #hide
    ```

Now that we have an `ODESolution` `sol`, we can call [`inspect`](@ref) to open the inspector GUI. The docstring provides several options to customize how the app is displayed.

```@example ndi
inspect(sol; reset=true)
sleep(1) # hide
define_timeseries!([ # hide
    (; selcomp=[EIndex(i) for i in 1:7], states=[:P]) # hide
]) # hide
sleep(3) # hide
NDI.save_electron_screenshot("screenshot.png") #hide
```
![screenshot](screenshot.png)


## Programmatic Access and GUI State Manipulation
Internally, the `NetworkDynamicsInspector` maintains a global reference to an `AppState` object. This AppState reflects changes made to the GUI by the user and can also be modified programmatically.

See the [NetworkDynamicsInspector API](@ref) for a complete list of available functions.
A good starting point is the [`dump_app_state`](@ref) function, which helps you recreate a GUI state that was previously configured manually.

Let's say we've adjusted the AppState to include an additional time series plot for the node states.

```@example ndi
set_state!(; t=1.75) #hide
define_timeseries!([ #hide
    (; selcomp=[VIndex(i) for i in 1:5], states=[:θ, :ω]) #hide
    (; selcomp=[EIndex(i) for i in 1:7], states=[:P]) #hide
]) #hide
sleep(3) #hide
nothing #hide
```

We can dump the code which helps us to recreate the app state:
```@example ndi
dump_app_state()
```

Now we can use this code to recreate the app state even though we've reseted it.
```@example ndi
buf = IOBuffer() #hide
dump_app_state(buf) #hide
code = String(take!(buf)) #hide
inspect(sol; reset=true)
sleep(1) #hide
eval(Meta.parse("begin;"*code*"end;")) #hide
sleep(3) #hide
NDI.save_electron_screenshot("screenshot2.png") #hide
"copy-paste and execute code returned by `dump_app_state` here" #hide
```
![screenshot](screenshot2.png)
