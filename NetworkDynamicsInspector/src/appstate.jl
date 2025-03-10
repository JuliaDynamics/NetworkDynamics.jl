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
    _lastclickel::Observable{Union{SymbolicCompIndex,Nothing}} = VIndex(1)
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
    _plotqueue::Channel{Task}
end
function AppState(sol::SciMLBase.AbstractODESolution)
    t = sol.t[begin]
    tmin = sol.t[begin]
    tmax = sol.t[end]
    graphplot = GraphPlot(sol)
    tsplots = OrderedDict(
        "ts-1" => TimeseriesPlot(),
    )
    active_tsplot = "ts-1"
    _plotqueue = Channel{Task}(Inf)
    AppState(sol, t, tmin, tmax, active_tsplot, graphplot, tsplots, _plotqueue)
end

function free_ts_key()
    tskeys = keys(APPSTATE[].tsplots[])
    used = Int[parse(Int, s[4:end]) for s in tskeys]
    i = 1
    while i ∈ used
        i += 1
    end
    return "ts-$i"
end

####
#### Programatic API
####
function appstate()
    isnothing(APPSTATE[]) && error("Uninitialized appstate. To initialize call `set_sol!(sol)` or `inspect(sol)`!")
    APPSTATE[]
end

# helper constrcut to mark undefined keyword arguments
struct NotSpecified end
function set_maybe!(obs::Observable, val)
    if obs[] != val
        obs[] = val
    end
end
set_maybe!(obs::Observable, ::NotSpecified) = nothing

"""
    set_sol!(sol)

Set the solution of the current appstate to `sol`.
"""
function set_sol!(sol)
    if isnothing(APPSTATE[])
        APPSTATE[] = AppState(sol)
    else
        APPSTATE[].sol[] = sol
    end
    nothing
end

"""
    set_state!(; sol, t, tmin, tmax)

Set the solution, current time and time limits of the current appstate.

To automaticially create commands see [`dump_app_state()`](@ref).
"""
function set_state!(; sol = NotSpecified(),
                      t = NotSpecified(),
                      tmin = NotSpecified(),
                      tmax = NotSpecified())
    sol != NotSpecified() && set_sol!(sol)
    set_maybe!(appstate().t, t)
    set_maybe!(appstate().tmin, tmin)
    set_maybe!(appstate().tmax, tmax)
    nothing
end

"""
    set_graphplot!(; nstate, estate, nstate_rel, estate_rel, ncolorrange, ecolorrange)

Set the properties of the graphplot of the current appstate.

To automaticially create commands see [`dump_app_state()`](@ref).
"""
function set_graphplot!(; nstate = NotSpecified(),
                          estate = NotSpecified(),
                          nstate_rel = NotSpecified(),
                          estate_rel = NotSpecified(),
                          ncolorrange = NotSpecified(),
                          ecolorrange = NotSpecified())
    gp = appstate().graphplot
    set_maybe!(gp.nstate, nstate)
    set_maybe!(gp.estate, estate)
    set_maybe!(gp.nstate_rel, nstate_rel)
    set_maybe!(gp.estate_rel, estate_rel)
    set_maybe!(gp.ncolorrange, ncolorrange)
    set_maybe!(gp.ecolorrange, ecolorrange)
    nothing
end

"""
    set_timeseries!(key; selcomp, states, rel)

Set properties of the timeseries plot with key `key`. See also [`define_timeseries!`](@ref).

To automaticially create commands see [`dump_app_state()`](@ref).
"""
function set_timeseries!(key; selcomp = NotSpecified(),
                              states = NotSpecified(),
                              rel = NotSpecified())
    if !haskey(appstate().tsplots[], key)
        appstate().tsplots[][key] = TimeseriesPlot()
    end
    tsplot = appstate().tsplots[][key]
    set_maybe!(tsplot.selcomp, selcomp)
    set_maybe!(tsplot.states, states)
    set_maybe!(tsplot.rel, rel)
    nothing
end

"""
    define_timeseries!(tsarray)

Defines timeseries, where `tsarray` is an array of timeseries keyword arguments
(see [`set_timeseries!`](@ref)).

To automaticially create commands see [`dump_app_state()`](@ref).
"""
function define_timeseries!(tsarray)
    if length(tsarray) != length(appstate().tsplots[])
        @warn "Due to current limitations, you need to reload the page if the number of timeseries plots changes"
        empty!(appstate().tsplots[])
        tskeys = [gendomid("ts") for _ in tsarray]
    else
        tskeys = keys(appstate().tsplots[])
    end
    for (key, tsargs) in zip(tskeys, tsarray)
        set_timeseries!(key; tsargs...)
    end
    nothing
end

"""
    dump_app_state()

Generate a list of [`set_sol!`](@ref), [`set_state!`](@ref), [`set_graphplot!`](@ref) and [`define_timeseries!`](@ref)
commands to recreate the current appstate.
The intended usecase is to quickly recreate "starting points" for interactive exploration.
"""
function dump_app_state()
    appstate()
    println("To recreate the current state, run the following commands:\n")
    println(styled"set_sol!({red:sol}) # optional if after inspect(sol)")
    println("set_state!(; t=$(appstate().t[]), tmin=$(appstate().tmin[]), tmax=$(appstate().tmax[]))")
    gp = appstate().graphplot
    println("set_graphplot!(; nstate=$(gp.nstate[]), estate=$(gp.estate[]), nstate_rel=$(gp.nstate_rel[]), estate_rel=$(gp.estate_rel[]), ncolorrange=$(gp.ncolorrange[]), ecolorrange=$(gp.ecolorrange[]))")
    println("define_timeseries!([")
    for ts in values(appstate().tsplots[])
        selstr = replace(repr(ts.selcomp[]), r"^.*\["=>"[")
        println("    (; selcomp=$(selstr), states=$(ts.states[]), rel=$(ts.rel[])),")
    end
    println("])")
    nothing
end
