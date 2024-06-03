using ForwardDiff
using ReverseDiff
using NetworkDynamics
using Graphs
using Test

N = 10 # number of nodes
k = 3  # average degree'
g = barabasi_albert(N, k) # a little more exciting than a bare random graph

@testset "Automatic Differentiation compatibility" begin
    @inline function diffusionedge!(e, v_s, v_d, p, t)
        e[1] = p[1] * (v_s[1] - v_d[1])
        nothing
    end

    @inline function diffusionvertex!(dv, v, esum, p, t)
        dv[1] = p[1] + esum[1]
        nothing
    end
    ### Constructing the network dynamics

    nd_diffusion_vertex = ODEVertex(; f=diffusionvertex!, dim=1, pdim=1)
    nd_diffusion_edge = StaticEdge(; f=diffusionedge!, dim=1, pdim=1, coupling=AntiSymmetric())
    nd = Network(g, nd_diffusion_vertex, nd_diffusion_edge)

    # p = (p_v, p_e)
    s = NWState(nd; utype=Float64, ptype=Float64)
    s.p.v[:,1] = collect(1:nv(g)) ./ nv(g) .- 0.5
    s.p.e[:,1] = 0.5 .* ones(ne(g))
    s.v[:,1] = randn(N)

    ## Jacobians wrt. states
    @test ForwardDiff.jacobian((out, x) -> nd(out, x, pflat(s), 0.0), ones(N), ones(N)) isa Matrix
    @test ReverseDiff.jacobian((out, x) -> nd(out, x, pflat(s), 0.0), ones(N), ones(N)) isa Matrix

    ## Gradient wrt. vertex parameters
    function gradpv(vp)
        # initializing the output array with the correct type is crucial
        out = similar(vp, N)
        _p = NWParameter(s.p, ptype=Vector{Union{eltype(vp), eltype(s.p)}})
        _p.v[:,1] = vp

        nd(out, s[:], _p[:], 0.0)
        sum(out)
    end

    @test ForwardDiff.gradient(gradpv, ones(N)) isa Vector
    @test ReverseDiff.gradient(gradpv, ones(N)) isa Vector
    vps = rand(N)
    @test ForwardDiff.gradient(gradpv, vps) == ReverseDiff.gradient(gradpv, vps)

    ## Gradient wrt. edge parameters
    function gradpe(ep)
        # initializing the output array with the correct type is crucial
        out = similar(ep, N)
        _p = NWParameter(s.p, ptype=Vector{Union{eltype(ep), eltype(s.p)}})
        _p.e[:,1] = ep

        nd(out, s[:], _p[:], 0.0)
        sum(out)
    end

    @test ForwardDiff.gradient(gradpe, ones(ne(g))) isa Vector
    @test ReverseDiff.gradient(gradpe, ones(ne(g))) isa Vector
    vps = rand(ne(g))
    @test ForwardDiff.gradient(gradpe, vps) == ReverseDiff.gradient(gradpe, vps)

    ## parameter array instead of tuple

    @inline function diffusionedge2!(e, v_s, v_d, p, t)
        e[1] = p[2] * (v_s[1] - v_d[1])
        nothing
    end

    @inline function diffusionvertex2!(dv, v, esum, p, t)
        dv[1] = p[1] - esum[1]
        nothing
    end
    nd_diffusion_vertex2 = ODEVertex(; f=diffusionvertex2!, dim=1, pdim=1)
    nd_diffusion_edge2 = StaticEdge(; f=diffusionedge2!, dim=1, pdim=1, coupling=AntiSymmetric())

    p = [1.0, 2.0]

    nd2 = Network(g, nd_diffusion_vertex2, nd_diffusion_edge2, g)
    nd2.f(x0, zeros(N), p, 0.0)

    function gradpnd2(p)
        # initializing the output array with the correct type is crucial
        out = similar(p, N)
        nd2.f(out, zeros(N), p, 0.0)
        sum(out)
    end

    @test ForwardDiff.gradient(gradpnd2, p) isa Vector
    @test ReverseDiff.gradient(gradpnd2, p) isa Vector
end
