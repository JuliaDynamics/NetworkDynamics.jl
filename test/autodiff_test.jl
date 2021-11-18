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
        e[1] = p * (v_s[1] - v_d[1])
        nothing
    end

    @inline function diffusionvertex!(dv, v, edges, p, t)
        dv[1] = p
        for e in edges
            dv[1] += e[1]
        end
        nothing
    end
    ### Constructing the network dynamics

    nd_diffusion_vertex = ODEVertex(; f=diffusionvertex!, dim=1)
    nd_diffusion_edge = StaticEdge(; f=diffusionedge!, dim=1, coupling=:antisymmetric)

    p_v = collect(1:nv(g)) ./ nv(g) .- 0.5
    p_e = 0.5 .* ones(ne(g))
    p = (p_v, p_e)

    nd = network_dynamics(nd_diffusion_vertex, nd_diffusion_edge, g; parallel=false)
    x0 = randn(N)


    ## Jacobians

    @test ForwardDiff.jacobian((out, x) -> nd(out, x, p, 0.0), ones(N), ones(N)) isa Matrix
    @test ReverseDiff.jacobian((out, x) -> nd(out, x, p, 0.0), ones(N), ones(N)) isa Matrix

    ## Gradient wrt. vertex parameters

    function gradpnd(p)
        # initializing the output array with the correct type is crucial
        out = similar(p, N)
        nd.f(out, zeros(N), (p, p_e), 0.0)
        sum(out)
    end

    @test ForwardDiff.gradient(gradpnd, ones(N)) isa Vector
    @test ReverseDiff.gradient(gradpnd, ones(N)) isa Vector

    ## parameter array instead of tuple

    @inline function diffusionedge2!(e, v_s, v_d, p, t)
        e[1] = p[2] * (v_s[1] - v_d[1])
        nothing
    end

    @inline function diffusionvertex2!(dv, v, edges, p, t)
        dv[1] = p[1]
        for e in edges
            dv[1] -= e[1]
        end

        nothing
    end
    nd_diffusion_vertex2 = ODEVertex(; f=diffusionvertex2!, dim=1)
    nd_diffusion_edge2 = StaticEdge(; f=diffusionedge2!, dim=1, coupling=:antisymmetric)

    p = [1.0, 2.0]

    nd2 = network_dynamics(nd_diffusion_vertex2, nd_diffusion_edge2, g; parallel=false)
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
