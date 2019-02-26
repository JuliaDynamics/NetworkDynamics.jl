module DynamicLines

using Parameters
using LightGraphs
using LinearAlgebra

export scalar_dynamic_line

"""
nodes!: vector of scalar functions node!(dx,x,l_s,l_t,p,t)
generates the dynamics on the nodes
lines!: vector of scalar functions line!(dl,l,x_s,x_t,p,t)
generates the dynamics of the edges
s_e: vector, i'th entry is the source node of the i'th edge
t_e: vector, i'th entry is the target node of the i'th edge
len_l: number of edges
len_n: number of nodes
"""
@with_kw struct scalar_dynamic_line
    nodes!
    lines!
    s_e
    t_e
    len_l
    len_n
end

'''Calling a struct of type scalar_dynamic_lines implements the ODE:
\frac{dx}{dt}=nodes(x,l_s,l_t,p,t)
\frac{dl}{dt}=lines(l,x_s,x_t,p,t)
l_s/l_t are vectors of arrays where the entries of the ith array
are the variables on the edges of which the ith node is the source/target.
dv=(dx_1,dx_2,...,dx_n,dl_1,...,dl_m)'''

function (d::scalar_dynamic_line)(dv, v, p, t)
    dl=0
    dx=1

    for i in 1:d.len_l
        dv[i+d.len_n] = d.lines![i](dl,v[i+d.len_n],v[d.s_e[i]],v[d.t_e[i]],p,t)
    end

    l_s = [
         [v[j+d.len_n] for j in 1:d.len_l if d.s_e[j] == i]
        for i in 1:d.len_n ]

    l_t = [
        [v[j+d.len_n] for j in 1:d.len_l if d.t_e[j] == i]
        for i in 1:d.len_n ]

    for i in 1:d.len_n
        #those are necessary as we sum over the entries of l_s/l_t and the sum
        #of an empty set is ill-defined in Julia unfortunately.
        if l_s[i] == []
            l_s[i] = [0]
        end
        if l_t[i] == []
            l_t[i] = [0]
        end

        dv[i]=d.nodes![i](dx,v[i],l_s[i],l_t[i],p,t)
    end
end


"""
When called with a graph, we construct the source and target vectors.
"""
function scalar_dynamic_line(nodes!, lines!, g::AbstractGraph)
    s_e = [src(e) for e in edges(g)]
    t_e = [dst(e) for e in edges(g)]
    len_n=length(nodes!)
    len_l=length(lines!)
    scalar_dynamic_line(nodes!, lines!, s_e, t_e, len_l, len_n)
end

end
