function gen_state_options(nw::Network, sidxs)
    options = OptionGroup{Symbol}[]
    isempty(sidxs) && return options
    groups = [
        ("Outputs & States", cf -> unique!(vcat(NetworkDynamics.outsym_flat(cf), sym(cf)))),
        ("Inputs", cf -> collect(NetworkDynamics.insym_all(cf))),
        ("Observables", obssym),
        ("Parameters", psym),
    ]
    exclusive_syms = Symbol[]
    for (label, getter) in groups
        common_syms = mapreduce(∩, sidxs) do sidx
            cf = nw[sidx]
            getter(cf)
        end
        push!(options, OptionGroup{Symbol}(label, common_syms))

        all_syms = mapreduce(∪, sidxs) do sidx
            cf = nw[sidx]
            getter(cf)
        end
        append!(exclusive_syms, setdiff(all_syms, common_syms))
    end
    if !isempty(exclusive_syms)
        push!(options, OptionGroup{Symbol}("Partially available syms", exclusive_syms))
    end
    options
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


NetworkDynamics.extract_nw(o::Observable) = extract_nw(o.val)
function NetworkDynamics.extract_nw(o::NamedTuple)
    if haskey(o, :sol)
        return extract_nw(o.sol)
    else
        error("No sol in NamedTuple")
    end
end

getcycled(v::AbstractVector, i) = v[mod1(i, length(v))]
gendomid(s::String) = replace(string(gensym("selectbox")), "#"=>"")

mutable struct Throttle{F}
    f::F
    delay::Float64
    @atomic prev::Float64
    @atomic future_task::Union{Timer,Nothing}
end
function Throttle(f, delay)
    Throttle(f, Float64(delay), 0.0, nothing)
end

function (tr::Throttle)(args...)
    now = time()
    if !isnothing(tr.future_task)
        close(tr.future_task)
        @atomic tr.future_task = nothing
    end
    if now - tr.prev > tr.delay
        @atomic tr.prev = now
        return tr.f(args...)
    else
        waitfor = tr.delay - (now - tr.prev) + 1e-3
        @atomic tr.future_task = Timer(waitfor) do _
            tr(args...)
        end
    end
end

function on_throttled(f, obs; update=false, delay)
    tr = Throttle(f, delay)
    on(obs; update) do args...
        tr(args...)
    end
end
function onany_throttled(f, obs...; update=false, delay)
    tr = Throttle(f, delay)
    onany(obs; update) do args...
        tr(args...)
    end
end
