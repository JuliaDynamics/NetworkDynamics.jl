
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
