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
    rand(rng, pdim(nd))
end
