module StaticLines

using Parameters
using LightGraphs
using LinearAlgebra

export scalar_static_line

"""
nodes!: vector of scalar functions node!(dx,x,l_s,l_t,p,t)
generates the dynamics on the nodes
lines!: vector of scalar functions line!(l,x_s,x_t,p,t)
s_e: vector, i'th entry is the source node of the i'th edge
t_e: vector, i'th entry is the target node of the i'th edge
len_l: number of edges
len_n: number of nodes
"""
@with_kw struct scalar_static_line
    nodes!
    lines!
    s_e
    t_e
    len_l
    len_n
end

'''Calling a struct of type scalar_static_line implements the ODE:
   \frac{dx}{dt} = nodes(x,l_s,l_t)
   l_s/l_t are vectors of arrays where the entries of the ith array
   are the variables on the edges of which the ith node is the source/target'''

function (d::scalar_static_line)(dx, x, p, t)
    l_e=zeros(d.len_l)
    l=0
    dx_temp=1

    for i in 1:d.len_l
        l_e[i] = d.lines![i](l,x[d.s_e[i]],x[d.t_e[i]],p,t)
    end

    l_s = [
         [l_e[j] for j in 1:d.len_l if d.s_e[j] == i]
        for i in 1:d.len_n ]

    l_t = [
        [l_e[j] for j in 1:d.len_l if d.t_e[j] == i]
        for i in 1:d.len_n ]

    for i in 1:d.len_n
        dx[i]=d.nodes![i](dx_temp,x[i],l_s[i],l_t[i],p,t)
    end
end


"""
When called with a graph, we construct the source and target vectors.
"""
function scalar_static_line(nodes!, lines!, g::AbstractGraph)
    s_e = [src(e) for e in edges(g)]
    t_e = [dst(e) for e in edges(g)]
    len_n=length(nodes!)
    len_l=length(lines!)
    scalar_static_line(nodes!, lines!, s_e, t_e, len_l, len_n)
end

end
