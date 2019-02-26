module StaticLines

using Parameters
using LightGraphs
using LinearAlgebra

export static_line_network_dynamics

"""
Documentation!!
"""
@with_kw struct static_line_network_dynamicst
    nodes!
    lines!
    s_e
    t_e
    l_e
    l_s
    l_t
    len_l
    len_n
end


function (d::static_line_network_dynamicst)(dx, x, p, t)
    l_temp=zeros(d.len_l)
    l=0
    dx_temp=1
    for i in 1:d.len_l
        l_temp[i] = d.lines![i](l,x[d.s_e[i]],x[d.t_e[i]],p,t)
    end

    l_s_temp = [
         [l_temp[j] for j in 1:d.len_l if d.s_e[j] == i]
        for i in 1:d.len_n ]

    l_t_temp = [
        [l_temp[j] for j in 1:d.len_l if d.t_e[j] == i]
        for i in 1:d.len_n ]

    for i in 1:d.len_n
        dx[i]=d.nodes![i](dx_temp,x[i],l_s_temp[i],l_t_temp[i],p,t)
    end
end


function static_line_network_dynamicst(nodes!, lines!, s_e, t_e)
    len_l = length(lines!)
    len_n = length(nodes!)
    l_e = zeros(len_l)

    l_s = [
         [l_e[j] for j in 1:len_l if s_e[j] == i]
        for i in 1:len_n ]
    l_t = [
        [l_e[j] for j in 1:len_l if t_e[j] == i]
        for i in 1:len_n ]

    static_line_network_dynamicst(nodes!, lines!, s_e, t_e, l_e, l_s, l_t, len_l, len_n)
end

"""
When called with a graph, we construct the source and target vectors.
"""
function static_line_network_dynamicst(nodes!, lines!, g::AbstractGraph)
    s_e = [src(e) for e in edges(g)]
    t_e = [dst(e) for e in edges(g)]
    static_line_network_dynamicst(nodes!, lines!, s_e, t_e)
end


g = barabasi_albert(10,3)
line! = (l,x_s,x_t,p,t) -> l = x_s - x_t
lines! = [line! for e in edges(g)]

# l_t and l_s is an array of arrays, we have to sum twice to sum over all elements.

function ssum(a)
    if a == []
        0
    else
        sum(a)
    end
end
node! = (dx,x,l_s,l_t,p,t) -> dx = -x-ssum(l_s) + ssum(l_t)
nodes! = [node! for n in vertices(g)]

sl_nd = static_line_network_dynamicst(nodes!, lines!, g)

x0 = rand(10)
dx0 = ones(10)

sl_nd(dx0,x0,nothing,0.)
sl_nd.l_s

using DifferentialEquations

sl_nd_prob = ODEProblem(sl_nd, x0, (0., 2.))

sl_nd_sol=solve(sl_nd_prob)

using Plots
plot(sl_nd_sol, legend=false)
end
