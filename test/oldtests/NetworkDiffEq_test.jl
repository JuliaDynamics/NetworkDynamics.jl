using Test
using Graphs
using NetworkDynamics

@testset "Test ODE_Static" begin
    # create a simpe ODE_Static model
    @inline function diffusionedge!(e, v_s, v_d, p, t)
        e[1] = v_s[1] - v_d[1]
        nothing
    end

    @inline function diffusionvertex!(dv, v, edges, p, t)
        dv[1] = 0.0
        for e in edges
            dv[1] += e[1]
        end
        nothing
    end
    nd_diffusion_vertex = ODEVertex(; f=diffusionvertex!, dim=1)
    nd_diffusion_edge = StaticEdge(; f=diffusionedge!, dim=1, coupling=:antisymmetric)

    # set up networ_dynamics on graph
    N = 10
    g = barabasi_albert(N, 3, 2)
    nd = network_dynamics(nd_diffusion_vertex, nd_diffusion_edge, g)

    x = rand(N)
    dx = similar(x)
    # first run to compile
    nd(dx, x, nothing, 0.0)

    # function barrier, otherwise there are still allocations
    function alloc(nd, dx, x)
        @allocated nd(dx, x, nothing, 0.0)
    end
    # test allocation free main loop
    @test 0 == alloc(nd, dx, x)

    # test prep_gd
    gd = nd(x, nothing, 0.0, GetGD)
    gs = nd(GetGS)
    # first test array of same type
    a = rand(length(x))
    gd_ret = NetworkDynamics.prep_gd(dx, a, gd, gs)
    @test gd_ret === gd
    @test gd.gdb.v_array == a

    # check with type missmatch
    a = rand(Int, length(x))
    dx = similar(a)
    gd_ret = NetworkDynamics.prep_gd(dx, a, gd, gs)
    @test gd_ret !== gd
    @test gd_ret.gdb.v_array == a
end


@testset "Test DDE_Static" begin
    # create a simple DDE_Static model
    @inline function diffusionedge!(e, v_s, v_d, p, t)
        e .= 0.1 * (v_s - v_d)
        nothing
    end

    @inline function diffusionvertex!(dv, v, edges, h_v, p, t)
        dv[1] = -h_v(t - 1., idxs=1)
        sum_coupling!(dv, edges)
        nothing
    end

    nd_diffusion_vertex = DDEVertex(; f=diffusionvertex!, dim=1)
    nd_diffusion_edge = StaticEdge(; f=diffusionedge!, dim=1, coupling=:antisymmetric)

    # set up network_dynamics on graph
    N = 20
    k = 8
    g = watts_strogatz(N, k, 0.0)

    nd = network_dynamics(nd_diffusion_vertex, nd_diffusion_edge, g)

    x = rand(N)
    dx = similar(x)
    p = nothing
    h(p, t; idxs=1) = x[idxs]

    # first run to compile
    nd(dx, x, h, p, 0.0)

    # function barrier, otherwise there are still allocations
    function alloc(nd, dx, x, h)
        @allocated nd(dx, x, h, p, 0.0)
    end
    # test allocation free main loop
    #@test 0 == alloc(nd, dx, x, h)
    alloc(nd, dx, x, h)

    # test prep_gd
    gd = nd(x, p, 0.0, GetGD)

    gs = nd(GetGS)
    # first test array of same type
    a = rand(length(x))
    gd_ret = NetworkDynamics.prep_gd(dx, a, gd, gs)
    @test gd_ret === gd
    @test gd.gdb.v_array == a

    # check with type missmatch
    a = rand(Int, length(x))
    dx = similar(a)
    gd_ret = NetworkDynamics.prep_gd(dx, a, gd, gs)
    @test gd_ret !== gd
    @test gd_ret.gdb.v_array == a
end

@testset "Test ODE_ODE " begin
    # create a simpe ODE_ODE model
    @inline function edge!(de, e, v_s, v_d, p, t)
        de[1] = e[1] + v_s[1] - v_d[1]
        de[2] = e[2] + v_d[1] - v_s[1]
        nothing
    end

    @inline function diffusionvertex!(dv, v, edges, p, t)
        dv[1] = 0.0
        for e in edges
            dv[1] += e[1]
        end
        nothing
    end
    nd_diffusion_vertex = ODEVertex(; f=diffusionvertex!, dim=1)
    nd_edge = ODEEdge(; f=edge!, dim=2, coupling=:fiducial)

    # set up networ_dynamics on graph
    N = 10
    g = barabasi_albert(N, 3, 2)
    nd = network_dynamics(nd_diffusion_vertex, nd_edge, g)

    x = rand(nv(g) + 2ne(g))
    dx = similar(x)
    # first run to compile
    nd(dx, x, nothing, 0.0)

    # function barrier, otherwise there are still allocations
    function alloc(nd, dx, x)
        @allocated nd(dx, x, nothing, 0.0)
    end
    # test allocation free main loop
    @test 0 == alloc(nd, dx, x)

    # test prep_gd
    gd = nd(x, nothing, 0.0, GetGD)
    gs = nd(GetGS)
    # first test array of same type
    a = rand(length(x))
    gd_ret = NetworkDynamics.prep_gd(dx, a, gd, gs)
    @test gd_ret === gd
    @test gd.gdb.v_array == a[1:gs.dim_v]
    @test gd.gdb.e_array == a[gs.dim_v+1:end]

    # check with type missmatch
    a = rand(Int, length(x))
    dx = similar(a)
    gd_ret = NetworkDynamics.prep_gd(dx, a, gd, gs)
    @test gd_ret !== gd
    @test gd_ret.gdb.v_array == a[1:gs.dim_v]
    @test gd_ret.gdb.e_array == a[gs.dim_v+1:end]
end
