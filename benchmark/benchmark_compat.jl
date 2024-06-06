#=
Helper functions to create compatible states between NetworkDynamics.jl and NetworkDynamics.jl v0.8
=#

if pkgversion(NetworkDynamics) < v"0.9.0"
    @info "Define compatibility functions"
    struct Network end
    struct SequentialExecution{P} end
    struct KAExecution{P} end
    struct NNlibScatter; f; end
    struct KAAggregator; f; end
    struct SequentialAggregator; f; end
    struct PolyesterAggregator; f; end
    function Network(g, v, e; execution = SequentialExecution{true}(), aggregator=SequentialAggregator(+))
        if execution isa KAExecution{true} && aggregator isa KAAggregator
            network_dynamics(v, e, g; parallel=true)
        elseif execution isa SequentialExecution && aggregator isa SequentialAggregator
            network_dynamics(v, e, g)
        else
            error("execution type not supported")
        end
    end
end

function _syms_old_order(nd::Network)
    syms = []
    for (i,cf) in enumerate(nd.im.vertexf)
        NetworkDynamics.isdynamic(cf) || continue
        append!(syms, collect(VIndex(i, 1:dim(cf))))
    end
    for (i,cf) in enumerate(nd.im.edgef)
        NetworkDynamics.isdynamic(cf) || continue
        append!(syms, collect(EIndex(i, 1:dim(cf))))
    end
    syms
end
function _psyms_old_order(nd::Network)
    syms = []
    for (i,cf) in enumerate(nd.im.vertexf)
        append!(syms, collect(VPIndex(i, 1:pdim(cf))))
    end
    for (i,cf) in enumerate(nd.im.edgef)
        append!(syms, collect(EPIndex(i, 1:pdim(cf))))
    end
    syms
end
function randx0(nd::Network)
    rng = StableRNG(1)
    _x0 = rand(rng, dim(nd))

    s = NWState(nd; utype=typeof(_x0))
    for (i,sym) in enumerate(_syms_old_order(nd::Network))
        s[sym] = _x0[i]
    end
    s[:]
end
randx0(nd::ODEFunction) = rand(StableRNG(1), length(nd.syms))
function legacy_order(nd::Network, _dx)
    s = NWState(nd, _dx)
    dx = similar(_dx)
    for (i,sym) in enumerate(_syms_old_order(nd::Network))
        dx[i] = s[sym]
    end
    dx
end
legacy_order(nd::ODEFunction, dx) = dx

function randp(nd::Network)
    rng = StableRNG(1)
    _p = rand(rng, pdim(nd))
    p = NWParameter(nd; ptype=typeof(_p))
    for (i, sym) in enumerate(_psyms_old_order(nd))
        p[sym] = _p[i]
    end
    p[:]
end
function randp(nd::ODEFunction)
    rng = StableRNG(1)
    g = nd.f.graph
    _p = rand(rng, nv(g) + ne(g))
    p = (_p[1:nv(g)], _p[nv(g)+1:end])
end
