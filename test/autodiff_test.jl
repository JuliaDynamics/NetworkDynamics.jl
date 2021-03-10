using ForwardDiff
using ReverseDiff
using NetworkDynamics
using LightGraphs
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

    nd_diffusion_vertex = ODEVertex(f! = diffusionvertex!, dim = 1)
    nd_diffusion_edge = StaticEdge(f! = diffusionedge!, dim = 1)

    p_v = collect(1:nv(g))./nv(g) .- 0.5
    p_e = .5 .* ones(ne(g))
    p  = (p_v, p_e)

    nd = network_dynamics(nd_diffusion_vertex, nd_diffusion_edge, g, parallel=false)
    x0 = randn(N)


    ## Jacobians

    @test ForwardDiff.jacobian((out, x) -> nd(out, x, p, 0.), ones(N), ones(N)) isa Matrix
    @test ReverseDiff.jacobian((out, x) -> nd(out, x, p, 0.), ones(N), ones(N)) isa Matrix

    ## Gradient wrt. vertex parameters

    function gradpnd(p)
        x0int = ones(N)
        nd.f(x0int, zeros(N), (p, p_e), 0.)
        sum(x0int)
    end

    # Forward diff throws a MethodError since multiplication with p transforms the x values
    # to dual numbers and those cant be stored in the type-specified x_array
    # Maybe this could be fixed by using more general types?

    @test_broken ForwardDiff.gradient(gradpnd, ones(N)) isa Vector # type problem
    @test_broken ReverseDiff.gradient(gradpnd, ones(N)) isa Vector

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
    nd_diffusion_vertex2 = ODEVertex(f! = diffusionvertex2!, dim = 1)
    nd_diffusion_edge2 = StaticEdge(f! = diffusionedge2!, dim = 1)

    p = [1., 2.]

    nd2 = network_dynamics(nd_diffusion_vertex2, nd_diffusion_edge2, g, parallel=false)
    nd2.f(x0, zeros(N), p, 0.)

    function gradpnd2(p)
        x0 = ones(N)
        nd2.f(x0, zeros(N), p, 0.)
        sum(x0)
    end

    @test_broken ForwardDiff.gradient(gradpnd2, p) isa Vector # type problem
    @test_borken ReverseDiff.gradient(gradpnd2, p) isa Vector
end
