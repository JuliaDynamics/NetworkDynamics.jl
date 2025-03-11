# Interactive Solution Inspection

Interactive solution inspection tool based on [WGLMakie](https://makie.org/website/) and [Bonito](https://github.com/SimonDanisch/Bonito.jl) are provided through the helper package `NetworkDynamicsInspector`.


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
inspect(sol; restart=false, reset=true)
define_timeseries!([
    (; selcomp=[EIndex(i) for i in 1:7], states=[:P])
])
NDI.save_electron_screenshot("screenshot.png") #hide
nothing #hide
```
![screenshot](screenshot.png)

## Programmatric Acces and GUI State manipulation
```@example ndi
set_state!(; t=2.0) #hide
define_timeseries!([ #hide
    (; selcomp=[VIndex(i) for i in 1:5], states=[:θ, :ω]) #hide
    (; selcomp=[EIndex(i) for i in 1:7], states=[:P]) #hide
]) #hide
dump_app_state()
```

```@example ndi
buf = IOBuffer() #hide
dump_app_state(buf) #hide
code = String(take!(buf)) #hide
inspect(sol; reset=true) #hide
eval(Meta.parse("begin;"*code*"end;")) #hide
NDI.save_electron_screenshot("screenshot2.png") #hide
nothing #hide
```
![screenshot](screenshot2.png)

