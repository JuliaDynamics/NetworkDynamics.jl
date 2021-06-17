using Test
using LightGraphs
using NetworkDynamics
using OrdinaryDiffEq
using DelayDiffEq
using LinearAlgebra

@inline function diffusion_edge!(e,v_s,v_d,p,t)
    e .= v_s .- v_d
    nothing
end

@inline function diff_dyn_edge!(de,e,v_s,v_d,p,t)
    de .= e .- (v_s .- v_d)
    nothing
end

@inline function diffusion_vertex!(dv, v, edges, p, t)
    dv .= 0.
    sum_coupling!(dv, edges) # Oriented sum of the incoming and outgoing edges
    nothing
end

@inline function diff_stat_vertex!(v, edges, p, t)
    v .= 0.
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


odevertex = ODEVertex(f! = diffusion_vertex!, dim = 1)
staticvertex = StaticVertex(f! = diff_stat_vertex!, dim = 1)
ddevertex = DDEVertex(f! = kuramoto_delay_vertex!, dim = 1)

staticedge = StaticEdge(f! = diffusion_edge!, dim = 1)
staticedge2 = StaticEdge(f! = diffusion_edge!, dim = 1, coupling = :antisymmetric)

odeedge = ODEEdge(f! = diff_dyn_edge!, dim = 1, coupling = :undirected)
sdedge = StaticDelayEdge(f! = kuramoto_delay_edge!, dim = 2, coupling = :antisymmetric)


vertex_list_1 = [odevertex for v in 1:10]
@test eltype(vertex_list_1) == ODEVertex{typeof(diffusion_vertex!)}

vertex_list_2 = [staticvertex for v in 1:10]
@test eltype(vertex_list_2) == StaticVertex{typeof(diff_stat_vertex!)}

vertex_list_3 = [ddevertex for v in 1:10]
@test eltype(vertex_list_3) == DDEVertex{typeof(kuramoto_delay_vertex!)}

vertex_list_4 = [v % 2 == 0 ? odevertex : staticvertex for v in 1:10]

vertex_list_5 = Array{VertexFunction}(vertex_list_4)
@test eltype(vertex_list_5) == VertexFunction

vertex_list_6 = Array{ODEVertex}(vertex_list_4)
@test eltype(vertex_list_6) == ODEVertex

vertex_list_7 = Array{DDEVertex}(vertex_list_4)
@test eltype(vertex_list_7) == DDEVertex

@test_throws MethodError Array{StaticVertex}(vertex_list_4) # this should error out


edge_list_1 = [staticedge for v in 1:25]
@test eltype(edge_list_1) == StaticEdge{typeof(diffusion_edge!)}

edge_list_1_2 = [staticedge2 for v in 1:25]
@test eltype(edge_list_1_2) == StaticEdge{typeof(diffusion_edge!)}
@test eltype(edge_list_1_2) <: StaticEdge

edge_list_2 = [odeedge for v in 1:25]
@test eltype(edge_list_2) <: ODEEdge

edge_list_3 = [sdedge for v in 1:25]
@test eltype(edge_list_3) <: StaticDelayEdge

edge_list_4 = [v % 2 == 0 ? odeedge : staticedge for v in 1:10]

edge_list_5 = Array{EdgeFunction}(edge_list_4)
@test eltype(edge_list_5) == EdgeFunction

# To fix this, either allow undefined ODEEdges, or convert static edges earlier to directed
# respectivel undirected
@test_broken edge_list_6 = Array{ODEEdge}(edge_list_4)
#@test eltype(edge_list_6) == ODEEdge

edgelist_7 = Array{StaticDelayEdge}(edge_list_1)
@test eltype(edgelist_7) == StaticDelayEdge

@test_throws MethodError Array{StaticDelayEdge}(edge_list_4) == StaticDelayEdge

@test_throws MethodError Array{StaticEdge}(edge_list_4) # this should error out

# network dynamics objects
g = barabasi_albert(10,5)

nd_diff_ode_static=network_dynamics(odevertex, staticedge, g)
nd_static_ode=network_dynamics(staticvertex, odeedge, g)
nd_dde_sd=network_dynamics(ddevertex, sdedge, g)

@test nd_diff_ode_static isa ODEFunction
@test nd_static_ode isa ODEFunction
@test nd_dde_sd isa DDEFunction

nd_1=network_dynamics(vertex_list_1,edge_list_1,g)
nd_2=network_dynamics(vertex_list_2,edge_list_2,g)
nd_3=network_dynamics(vertex_list_3,edge_list_3,g)

nd_4=network_dynamics(vertex_list_1,staticedge,g)
nd_5=network_dynamics(odevertex,edge_list_1,g)

@test nd_1 isa ODEFunction
@test nd_2 isa ODEFunction
@test nd_3 isa DDEFunction
@test nd_4 isa ODEFunction
@test nd_5 isa ODEFunction
