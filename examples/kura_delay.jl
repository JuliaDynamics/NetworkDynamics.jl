using DelayDiffEq
using Random
using LightGraphs
using NetworkDynamics
using Distributions
using BenchmarkTools

N = 6
g = SimpleDiGraph(complete_graph(N))


function kuramoto!(dθ, θ, edges, p, t)
    ω = p
    dθ[1] = ω
    for e in edges
        dθ[1] += e[1]
    end
    return nothing
end

function delay_coupling(e, θ_s, θ_d, uh, θ_s_idxs, θ_d_idxs, p, global_p, t)
    hist1 = uh(t - p[1]; idxs=1)
    e[1] = p[2] * sin(hist1 - θ_d[1])
    return nothing
end

vertex = ODEVertex(f! = kuramoto!, dim = 1)
edge = StaticDelayEdge(f! = delay_coupling, dim=1, coupling=:directed)
nd = network_dynamics(vertex,edge,g)

Random.seed!(1)

ω = rand(nv(g))
τ = rand(ne(g))
minimum(τ)
κ = rand(ne(g))
p = (ω, [τ'; κ'])

θ₀ = rand(Uniform(0, 2π), N)

const past = rand(Uniform(0, 2π), N)
h(p, t; idxs=nothing) = past[idxs]

x0 = ones(N)
dx0 = similar(x0)
@allocated nd.f(dx0, x0,h,p,0) # Works as expected
@allocated nd.f(dx0, x0,h,p,0) # Works as expected
@btime $nd.f($dx0, $x0,$h,$p,0) # Works as expected

prob = DDEProblem(nd, θ₀, h, (0.0, 1.0), p, constant_lags=τ)

@btime solve(prob, MethodOfSteps(BS3()), saveat=0.01)

#### Other system

using DelayDiffEq
using Distributions
using Random

function f!(dθ, θ, h::H, p, t) where H
    ω, A = p
    n = length(θ)
    lags = reshape(lag, n,n)
    @inbounds for j in 1:n
        coupling = 0.0
        @inbounds for i in 1:n
            coupling += foo(dθ, θ, h::H, p, t, i, j)
        end
        dθ[j] = ω[j] + coupling
    end
    nothing
end

function foo(dθ, θ, h::H, p, t, i, j) where H
    ω, A = p
    n = length(θ)
    lags = reshape(lag, n,n)
    hist = h(p, t-lags[i,j]; idxs=i)
    return A[i,j]*sin(hist - θ[j])
end

n = 10
Random.seed!(1)
ω = rand(n)
A = rand(n,n)
const lag = rand(n*n)
θ₀ = rand(Uniform(0, 2π), n)
p = (ω, A)
const past = rand(Uniform(0, 2π), n)
h(p, t; idxs=nothing) = typeof(idxs) <: Number ? past[idxs] : past

prob = DDEProblem(f!, θ₀, h, (0.0, 1.0), p, constant_lags=lag)

@btime solve(prob, MethodOfSteps(BS3()), saveat=0.01, reltol=0, abstol=1e-5)
