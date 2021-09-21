using Test
using LightGraphs
using NetworkDynamics

@testset "Test nd_ODE_Static Module" begin
    # create a simpe nd_ODE_Static model
    @inline function diffusionedge!(e, v_s, v_d, p, t)
        e[1] = v_s[1] - v_d[1]
        nothing
    end

    @inline function diffusionvertex!(dv, v, edges, p, t)
        dv[1] = 0.
        for e in edges
            dv[1] += e[1]
        end
        nothing
    end
    nd_diffusion_vertex = ODEVertex(f! = diffusionvertex!, dim = 1)
    nd_diffusion_edge = StaticEdge(f! = diffusionedge!, dim = 1, coupling = :antisymmetric)

    # set up networ_dynamics on graph
    N = 10
    g = barabasi_albert(N, 3, 2)
    nd = network_dynamics(nd_diffusion_vertex, nd_diffusion_edge, g)

    x = rand(N)
    dx = similar(x)
    # first run to compile
    nd(dx, x, nothing, 0.)

    # function barrier, otherwise there are still allocations
    function alloc(nd, dx, x)
        @allocated nd(dx, x, nothing, 0.)
    end
    # test allocation free main loop
    @test 0 == alloc(nd, dx, x)

    # test prep_gd
    gd = nd(x, nothing, 0., GetGD)
    gs = nd(GetGS)
    # first test array of same type
    a = rand(length(x))
    gd_ret = nd_ODE_DDE_combined_mod.prep_gd(dx, a, gd, gs)
    @test gd_ret === gd
    @test gd.gdb.v_array == a

    # check with type missmatch
    a = rand(Int, length(x))
    dx = similar(a)
    gd_ret = nd_ODE_DDE_combined_mod.prep_gd(dx, a, gd, gs)
    @test gd_ret !== gd
    @test gd_ret.gdb.v_array == a
end
