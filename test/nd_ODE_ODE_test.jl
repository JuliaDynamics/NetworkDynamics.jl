using Test
using LightGraphs
using NetworkDynamics

@testset "Test nd_ODE_ODE Module" begin
    # create a simpe nd_ODE_ODE model
    @inline function edge!(de, e, v_s, v_d, p, t)
        de[1] = e[1] + v_s[1] - v_d[1]
        de[2] = e[2] + v_d[1] - v_s[1]
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
    nd_edge = ODEEdge(f! = edge!, dim = 2, coupling = :fiducial)

    # set up networ_dynamics on graph
    N = 10
    g = barabasi_albert(N, 3, 2)
    nd = network_dynamics(nd_diffusion_vertex, nd_edge, g)

    x = rand(nv(g)+2ne(g))
    dx = similar(x)
    # first run to compile
    nd(dx, x, nothing, 0.)

    # test allocation free main loop
    @test 0 == @allocated nd(dx, x, nothing, 0.)

    # test prep_gd
    gd = nd(x, nothing, 0., GetGD)
    gs = nd(GetGS)
    # first test array of same type
    a = similar(x)
    gd_ret = nd_ODE_ODE_mod.prep_gd(dx, a, gd, gs)
    @test gd_ret === gd
    @test gd.gdb.v_array == a[1:gs.dim_v]
    @test gd.gdb.e_array == a[gs.dim_v+1:end]

    # check with type missmatch
    a = rand(Int, length(x))
    dx = similar(a)
    gd_ret = nd_ODE_ODE_mod.prep_gd(dx, a, gd, gs)
    @test gd_ret !== gd
    @test gd_ret.gdb.v_array == a[1:gs.dim_v]
    @test gd_ret.gdb.e_array == a[gs.dim_v+1:end]
end
