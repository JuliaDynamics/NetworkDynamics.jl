module DynamicLines

using Parameters
using LightGraphs
using LinearAlgebra

export scalar_dynamic_line

"""
Documentation!!
"""
@with_kw struct scalar_dynamic_line
    nodes!
    lines!
    s_e
    t_e
    len_l
    len_n
end


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
