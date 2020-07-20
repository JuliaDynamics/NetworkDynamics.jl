using Test
import NetworkDynamics.ComponentFunctions
using NetworkDynamics.ComponentFunctions

printstyled("--- Function Typology --- \n", bold=true, color=:white)

@inline function diffusion_edge!(e,v_s,v_d,p,t)
    e .= v_s .- v_d
    nothing
end

@inline function diff_dyn_edge!(de,e,v_s,v_d,p,t)
    de .= e .- (v_s .- v_d)
    nothing
end

@inline function diffusion_vertex!(dv, v, e_s, e_d, p, t)
    dv .= 0.
    oriented_edge_sum!(dv, e_s, e_d) # Oriented sum of the incoming and outgoing edges
    nothing
end

@inline function diff_stat_vertex!(v, e_s, e_d, p, t)
    v .= 0.
    nothing
end

@inline function kuramoto_delay_edge!(e, v_s, v_d, h_v_s, h_v_d, p, t)
    # The coupling is no longer symmetric, so we need to store BOTH values (see tutorials for details)
    e[1] = p * sin(v_s[1] - h_v_d[1])
    e[2] = p * sin(h_v_s[1] - v_d[1])
    nothing
end

@inline function kuramoto_delay_vertex!(dv, v, e_s, e_d, h_v, p, t)
    dv[1] = p
    for e in e_s
        dv[1] -= e[1]
    end
    for e in e_d
        dv[1] -= e[2]
    end
    nothing
end


odevertex = ODEVertex(f! = diffusion_vertex!, dim = 1)
staticvertex = StaticVertex(f! = diff_stat_vertex!, dim = 1)
staticedge = StaticEdge(f! = diffusion_edge!, dim = 1)
odeedge = ODEEdge(f! = diff_dyn_edge!, dim = 1)
# dim of the edge is 2 since the coupling is not symmetric
sdedge = StaticDelayEdge(f! = kuramoto_delay_edge!, dim = 2)
ddevertex = DDEVertex(f! = kuramoto_delay_vertex!, dim = 1)



vertex_list_1 = [odevertex for v in 1:10]
@test eltype(vertex_list_1) == ODEVertex{typeof(diffusion_vertex!)}

vertex_list_2 = [staticvertex for v in 1:10]
@test eltype(vertex_list_2) == StaticVertex{typeof(diff_stat_vertex!)}



vertex_list_3 = [v % 2 == 0 ? odevertex : staticvertex for v in 1:10]
vertex_list_4 = Array{VertexFunction}(vertex_list_3)
@test eltype(vertex_list_4) == VertexFunction

vertex_list_5 = Array{ODEVertex}(vertex_list_3)
@test eltype(vertex_list_5) == ODEVertex


@test_throws MethodError Array{StaticVertex}(vertex_list_3) # this should error out

edge_list_1 = [odeedge for v in 1:10]
@test eltype(edge_list_1) == ODEEdge{typeof(diff_dyn_edge!)}

edge_list_2 = [staticedge for v in 1:10]
edge_list_3 = [v % 2 == 0 ? odeedge : staticedge for v in 1:10]
edge_list_4 = Array{EdgeFunction}(edge_list_3)
edge_list_5 = Array{ODEEdge}(edge_list_3)
@test_throws MethodError Array{StaticEdge}(edge_list_3) # this should error out
