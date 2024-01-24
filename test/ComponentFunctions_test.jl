using Test
using Graphs
using NetworkDynamics
using OrdinaryDiffEq
using DelayDiffEq
using LinearAlgebra

@testset "Unique edges" begin
    N=3
    g=complete_graph(N)

    f = (e, v_s, v_d, p, t) -> begin
        e .= v_s .- v_d
        nothing
    end

    v = (dv,v,edges,p,t) -> begin
        for e in edges
            dv .+= e
        end
    end

    j = (a,b,c,d,e,f,g) -> nothing

    eundir = StaticEdge(; f=f, dim=2, coupling = :undirected)
    eundef = StaticEdge(; f=f, dim=2)
    deundef = StaticDelayEdge(;f=j, dim=2)

    mixed_edges = [eundef, eundef, deundef]

    vertex = ODEVertex(; f=v, dim=1)

    ndundir = network_dynamics(vertex,eundir, g)
    @test length(ndundir.f.unique_edges!) == 1

    ndundef = network_dynamics(vertex,eundef, g)
    @test length(ndundef.f.unique_edges!) == 1

    nddirundef = network_dynamics(vertex,eundef, SimpleDiGraph(g))
    @test length(nddirundef.f.unique_edges!) == 1

    # The wrapper created for :undefined -> :undirected reconstruction interfers with
    # uniqueness of edge functions
    ndmix = network_dynamics(vertex, mixed_edges, g)
    @test_broken length(ndmix.f.unique_edges!) == 2

    nddirmix = network_dynamics(vertex, repeat(mixed_edges,2), SimpleDiGraph(g))
    @test length(nddirmix.f.unique_edges!) == 2
end


@testset "StaticEdge constructor" begin
    f = (e, v_s, v_d, p, t) -> begin
        e .= v_s .- v_d
        nothing
    end

    @test_throws ErrorException StaticEdge(f, 1, :unspecified, [:e])

    @test_throws ErrorException StaticEdge(f, 0, :undefined, [:e])
    @test StaticEdge(f, 1, :undefined, [:e]) isa StaticEdge


    @test_throws ErrorException StaticEdge(f, 3, :fiducial, [:e, :e, :e])
    @test_throws ErrorException StaticEdge(f, 2, :fiducial, [:e])
    @test StaticEdge(f, 2, :fiducial, [:e, :e]) isa StaticEdge

    @test StaticEdge(f, 1, :undirected, [:e]) isa StaticEdge

    fundir = StaticEdge(; f=f, dim=2, coupling=:undirected)
    fanti = StaticEdge(; f=f, dim=2, coupling=:antisymmetric)
    fdir = StaticEdge(; f=f, dim=2, coupling=:directed)
    fsym = StaticEdge(; f=f, dim=2, coupling=:symmetric)
    ffid = StaticEdge(; f=f, dim=2, coupling=:fiducial)

    x = rand(2)
    y = rand(2)

    eundir = zeros(4)
    eanti  = zeros(4)
    esym   = zeros(4)
    edir   = zeros(2)
    efid   = zeros(2)

    fundir.f(eundir, x, y, nothing, nothing)
    fanti.f(eanti, x, y, nothing, nothing)
    fsym.f(esym, x, y, nothing, nothing)
    fdir.f(edir, x, y, nothing, nothing)
    ffid.f(efid, x, y, nothing, nothing)

    @test eundir == eanti
    @test eanti[1:2] == -eanti[3:4]
    @test esym[1:2] == esym[3:4]
    @test eundir[1:2] == edir
    @test edir == efid
end

@testset "StaticDelayEdge constructor" begin
    f = (e, v_s, v_d, h_v_s, h_v_d, p, t) -> begin
        hist_v_s = h_v_s(t - 1., idxs = 1)
        hist_v_d = h_v_d(t - 1., idxs = 1)
        e[1] = hist_v_s - hist_v_d
        nothing
    end

    @test_throws ErrorException StaticDelayEdge(f, 1, :unspecified, [:e])

    @test_throws ErrorException StaticDelayEdge(f, 0, :undefined, [:e])
    @test StaticDelayEdge(f, 1, :undefined, [:e]) isa StaticDelayEdge


    @test_throws ErrorException StaticDelayEdge(f, 3, :fiducial, [:e, :e, :e])
    @test_throws ErrorException StaticDelayEdge(f, 2, :fiducial, [:e])
    @test StaticDelayEdge(f, 2, :fiducial, [:e, :e]) isa StaticDelayEdge

    @test StaticDelayEdge(f, 1, :undirected, [:e]) isa StaticDelayEdge

    fundir = StaticDelayEdge(; f=f, dim=2, coupling=:undirected)
    fdir = StaticDelayEdge(; f=f, dim=2, coupling=:directed)
    ffid = StaticDelayEdge(; f=f, dim=2, coupling=:fiducial)

    x = rand(2)
    y = rand(2)

    eundir = zeros(4)
    edir   = zeros(2)
    efid   = zeros(2)

    fundir.f(eundir, nothing, nothing, (t; idxs) -> 1., (t; idxs) -> 1., nothing, 0.)
    fdir.f(edir, nothing, nothing, (t; idxs) -> 1., (t; idxs) -> 1., nothing, 0.)
    ffid.f(efid, nothing, nothing, (t; idxs) -> 1., (t; idxs) -> 1., nothing, 0.)


    @test eundir[1:2] == edir
    @test edir == efid
end


@testset "ODEEdge constructor" begin
    f = (de, e, v_s, v_d, p, t) -> begin
        de .= v_s .- v_d
        nothing
    end
    MM = 1

    @test_throws ErrorException ODEEdge(f, 1, :unspecified, MM, [:e])
    @test_throws ErrorException ODEEdge(f, 1, :symmetric, MM, [:e])
    @test_throws ErrorException ODEEdge(f, 1, :antisymmetric, MM, [:e])
    @test_throws ErrorException ODEEdge(f, 1, :undefined, MM, [:e])

    @test ODEEdge(f, 1, :undirected, MM, [:e]) isa ODEEdge
    @test_throws ErrorException ODEEdge(f, 0, :undirected, MM, [:e])

    # MM test

    @test ODEEdge(f, 1, :undirected, 1, [:e]).mass_matrix == 1
    @test ODEEdge(f, 1, :undirected, I, [:e]).mass_matrix == I
    @test ODEEdge(f, 1, :undirected, [0.0], [:e]).mass_matrix == [0.0, 0.0]
    @test ODEEdge(f, 1, :undirected, zeros(1, 1), [:e]).mass_matrix == zeros(2, 2)

    @test_throws ErrorException ODEEdge(f, 3, :fiducial, MM, [:e, :e, :e])
    @test_throws ErrorException ODEEdge(f, 2, :fiducial, MM, [:e])
    @test ODEEdge(f, 2, :fiducial, MM, [:e, :e]) isa ODEEdge

    @test ODEEdge(f, 1, :undirected, MM, [:e]) isa ODEEdge

    fundir = ODEEdge(; f=f, dim=2, coupling=:undirected)
    fdir = ODEEdge(; f=f, dim=2, coupling=:directed)
    ffid = ODEEdge(; f=f, dim=2, coupling=:fiducial)

    x = rand(2)
    y = rand(2)

    eundir = ones(4)
    edir   = ones(2)
    efid   = ones(2)

    deundir = zeros(4)
    dedir   = zeros(2)
    defid   = zeros(2)

    fundir.f(deundir, eundir, x, y, nothing, nothing)
    fdir.f(dedir, edir, x, y, nothing, nothing)
    ffid.f(defid, efid, x, y, nothing, nothing)

    @test deundir[1:2] == -deundir[3:4]
    @test eundir[1:2] == edir
    @test edir == efid
end

@testset "Function Typology" begin
    @inline function diffusion_edge!(e, v_s, v_d, p, t)
        e .= v_s .- v_d
        nothing
    end

    @inline function diff_dyn_edge!(de, e, v_s, v_d, p, t)
        de .= e .- (v_s .- v_d)
        nothing
    end

    @inline function diffusion_vertex!(dv, v, edges, p, t)
        dv .= 0.0
        sum_coupling!(dv, edges) # Oriented sum of the incoming and outgoing edges
        nothing
    end

    @inline function diff_stat_vertex!(v, edges, p, t)
        v .= 0.0
        nothing
    end

    @inline function kuramoto_delay_edge!(e, v_s, v_d, h_v_s, h_v_d, p, t)
        e[1] = p * sin(v_s[1] - h_v_d[1])
        nothing
    end

    @inline function kuramoto_delay_vertex!(dv, v, edges, h_v, p, t)
        dv[1] = p
        for e in edges
            dv[1] -= e[1]
        end
        nothing
    end


    odevertex = ODEVertex(; f=diffusion_vertex!, dim=1)
    staticvertex = StaticVertex(; f=diff_stat_vertex!, dim=1)
    ddevertex = DDEVertex(; f=kuramoto_delay_vertex!, dim=1)

    staticedge = StaticEdge(; f=diffusion_edge!, dim=1, coupling=:antisymmetric)
    odeedge = ODEEdge(; f=diff_dyn_edge!, dim=1, coupling=:undirected)
    sdedge = StaticDelayEdge(; f=kuramoto_delay_edge!, dim=2, coupling=:undirected)


    vertex_list_1 = [odevertex for v in 1:10]
    @test eltype(vertex_list_1) == ODEVertex{typeof(diffusion_vertex!), typeof(nothing)}

    vertex_list_2 = [staticvertex for v in 1:10]
    @test eltype(vertex_list_2) <: ODEVertex

    vertex_list_3 = [ddevertex for v in 1:10]
    @test eltype(vertex_list_3) == DDEVertex{typeof(kuramoto_delay_vertex!)}

    vertex_list_4 = [v % 2 == 0 ? odevertex : staticvertex for v in 1:10]

    vertex_list_5 = Array{VertexFunction}(vertex_list_4)
    @test eltype(vertex_list_5) == VertexFunction

    vertex_list_6 = Array{ODEVertex}(vertex_list_4)
    @test eltype(vertex_list_6) == ODEVertex

    @test_throws MethodError Array{DDEVertex}(vertex_list_4) # Methods for conversion of ODEVertx to DDEVertex are deprecated

    @test_throws TypeError Array{StaticVertex}(vertex_list_4) # this should error out since StaticVertex is not a type


    edge_list_1 = [staticedge for v in 1:25]
    #@test eltype(edge_list_1) == StaticEdge{typeof(diffusion_edge!)}
    @test eltype(edge_list_1) <: StaticEdge

    edge_list_2 = [odeedge for v in 1:25]
    @test eltype(edge_list_2) <: ODEEdge

    edge_list_3 = [sdedge for v in 1:25]
    #@test eltype(edge_list_3) == StaticDelayEdge{typeof(kuramoto_delay_edge!)}
    @test eltype(edge_list_3) <: StaticDelayEdge

    edge_list_4 = [v % 2 == 0 ? odeedge : staticedge for v in 1:10]

    edge_list_5 = Array{EdgeFunction}(edge_list_4)
    @test eltype(edge_list_5) == EdgeFunction

    # To fix this, either allow undefined ODEEdges, or convert static edges earlier to directed
    # respectivel undirected
    @test_broken edge_list_6 = Array{ODEEdge}(edge_list_4)
    #@test eltype(edge_list_6) == ODEEdge

    @test_throws MethodError Array{StaticDelayEdge}(edge_list_4)

    @test_throws MethodError Array{StaticEdge}(edge_list_4) # this should error out

    # network dynamics objects
    g = barabasi_albert(10, 5)

    nd_diff_ode_static = network_dynamics(odevertex, staticedge, g)
    nd_static_ode = network_dynamics(staticvertex, odeedge, g)
    nd_dde_sd = network_dynamics(ddevertex, sdedge, g)

    @test nd_diff_ode_static isa ODEFunction
    @test nd_static_ode isa ODEFunction
    @test nd_dde_sd isa DDEFunction

    nd_1 = network_dynamics(vertex_list_1, edge_list_1, g)
    nd_2 = network_dynamics(vertex_list_2, edge_list_2, g)
    nd_3 = network_dynamics(vertex_list_3, edge_list_3, g)

    nd_4 = network_dynamics(vertex_list_1, staticedge, g)
    nd_5 = network_dynamics(odevertex, edge_list_1, g)

    @test nd_1 isa ODEFunction
    @test nd_2 isa ODEFunction
    @test nd_3 isa DDEFunction
    @test nd_4 isa ODEFunction
    @test nd_5 isa ODEFunction
end


@testset "Allocations of promoted Vertices" begin
    f = (v, edges, p, t) -> v .= pi
    osv = StaticVertex(; f=f, dim=1)
    x = [1.0]
    dx = [0.0]
    foo! = osv.f
    foo!(dx, x, nothing, nothing, nothing)
    @test (@allocated foo!(dx, x, nothing, nothing, nothing)) == 0

    g! = (e, v_s, v_d, p, t) -> e .= pi
    se = StaticEdge(; f=g!, dim=1, coupling=:undirected)
    ose = ODEEdge(se)
    e = [1.0, 1.0]
    de = [0.0, 0.0]
    bar! = ose.f
    bar!(de, e, nothing, nothing, nothing, nothing)
    @test (@allocated bar!(de, e, nothing, nothing, nothing, nothing)) == 0

end


@testset begin "DimWarn"
    f = (v, edges, p, t) -> v .= pi
    g = (e, v_s, v_d, p, t) -> e .= pi
    h = (a,b,c,d,e,f) -> nothing
    j = (a, b, c, d, e, f, g) -> nothing

    @test_throws ErrorException StaticVertex(; f=g, dim=3, sym=[:θ, :ω, :u])
    @test StaticVertex(; f=f, dim=3, sym=[:θ, :ω, :u]) isa VertexFunction

    @test_throws ErrorException StaticVertex(; f=f, dim=5, sym=[:θ, :ω, :u])
    @test_throws ErrorException StaticVertex(; f=f, dim=1, sym=[:θ, :ω, :u])
    @test_throws ErrorException StaticVertex(; f=f, dim=-1)

    @test_throws ErrorException ODEVertex(; f=f, dim=3, sym=[:θ, :ω, :u])
    @test ODEVertex(; f=g, dim=3, sym=[:θ, :ω, :u]) isa VertexFunction


    @test_throws ErrorException ODEVertex(; f=g, dim=5, sym=[:θ, :ω, :u])
    @test_throws ErrorException ODEVertex(; f=g, dim=1, sym=[:θ, :ω, :u])
    @test_throws ErrorException ODEVertex(; f=g, dim=-1)

    @test_throws ErrorException DDEVertex(; f=g, dim=3, sym=[:θ, :ω, :u])
    @test DDEVertex(; f=h, dim=3, sym=[:θ, :ω, :u]) isa VertexFunction


    @test_throws ErrorException DDEVertex(; f=h, dim=5, sym=[:θ, :ω, :u])
    @test_throws ErrorException DDEVertex(; f=h, dim=1, sym=[:θ, :ω, :u])
    @test_throws ErrorException DDEVertex(; f=h, dim=-1)


    @test_throws ErrorException StaticEdge(; f=f, dim=3, sym=[:θ, :ω, :u])
    @test StaticEdge(; f=g, dim=3, sym=[:θ, :ω, :v], coupling=:undirected) isa EdgeFunction

    @test_throws ErrorException StaticEdge(; f=g, dim=5, sym=[:θ, :ω, :u])
    @test_throws ErrorException StaticEdge(; f=g, dim=1, sym=[:θ, :ω, :u])
    @test_throws ErrorException StaticEdge(; f=g, dim=-1)


    @test_throws ErrorException ODEEdge(; f=f, dim=3, sym=[:θ, :ω, :u])
    @test ODEEdge(; f=h, dim=3, sym=[:θ, :ω, :v], coupling=:undirected) isa EdgeFunction

    @test_throws ErrorException ODEEdge(; f=h, dim=5, sym=[:θ, :ω, :u])
    @test_throws ErrorException ODEEdge(; f=h, dim=1, sym=[:θ, :ω, :u])
    @test_throws ErrorException ODEEdge(; f=h, dim=-1)

    @test_throws ErrorException StaticDelayEdge(; f=h, dim=3, sym=[:θ, :ω, :u])
    @test StaticDelayEdge(; f=j, dim=3, sym=[:θ, :ω, :u]) isa EdgeFunction

    @test_throws ErrorException ODEEdge(; f=j, dim=5, sym=[:θ, :ω, :u])
    @test_throws ErrorException ODEEdge(; f=j, dim=1, sym=[:θ, :ω, :u])
    @test_throws ErrorException ODEEdge(; f=j, dim=-1)
end
