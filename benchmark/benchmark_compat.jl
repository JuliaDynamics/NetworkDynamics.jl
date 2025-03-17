#=
Helper functions to create compatible states between newer and older versions of NetworkDynamics.jl
for cross version benchmarking.
=#

function randx0(nd::Network)
    rng = StableRNG(1)
    rand(rng, dim(nd))
end

function randp(nd::Network)
    rng = StableRNG(1)
    _p = rand(rng, pdim(nd))
    p = NWParameter(nd; ptype=typeof(_p))
    for (i, sym) in enumerate(_psyms_old_order(nd))
        p[sym] = _p[i]
    end
    p[:]
end
