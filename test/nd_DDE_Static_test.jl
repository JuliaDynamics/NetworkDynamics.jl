using Test
using LightGraphs
using NetworkDynamics

@testset "Test nd_DDE_Static Module" begin
    # create a simple nd_DDE_Static model
    @inline function diffusionedge!(e, v_s, v_d, p, t)
        e .= .1 * (v_s - v_d)
        nothing
    end

    @inline function diffusionvertex!(dv, v, edges, h_v, p, t)
        dv .= -h_v
        sum_coupling!(dv, edges)
        nothing
    end

    nd_diffusion_vertex = DDEVertex(f! = diffusionvertex!, dim = 1)
    nd_diffusion_edge = StaticEdge(f! = diffusionedge!, dim = 1, coupling = :antisymmetric)

    # set up network_dynamics on graph
    N = 20
    k = 8
    g = watts_strogatz(N, k, 0.)

    nd = network_dynamics(nd_diffusion_vertex, nd_diffusion_edge, g)

    x = rand(N)
    dx = similar(x)
    p = (nothing, nothing, 1.0)
    h(out, p, t) = (out .= x)

    # first run to compile
    nd(dx, x, h, p, 0.0) # signature found in nd_DDE_Static.jl line 50

    # function barrier, otherwise there are still allocations
    function alloc(nd, dx, x, h)
        @allocated nd(dx, x, h, p, 0.0)
    end
    # test allocation free main loop
    #@test 0 == alloc(nd, dx, x, h)
    alloc(nd, dx, x, h)

    # test prep_gd
    gd = nd(x, p, 0., GetGD)

    gs = nd(GetGS)
    # first test array of same type
    a = rand(length(x))
    gd_ret = nd_DDE_Static_mod.prep_gd(dx, a, gd, gs)
    @test gd_ret === gd
    @test gd.gdb.v_array == a

    # check with type missmatch
    a = rand(Int, length(x))
    dx = similar(a)
    gd_ret = nd_DDE_Static_mod.prep_gd(dx, a, gd, gs)
    @test gd_ret !== gd
    @test gd_ret.gdb.v_array == a
end
