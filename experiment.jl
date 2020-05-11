using LightGraphs
using OrdinaryDiffEq
using Plots

function node!(dv, v, p, t, inputs)
    dv = -inputs["Layer1"]
end

function interact(v_in, v_out)
    v_in - v_out
end

function aggregate(node_interactions)
  sum(node_interactions)
end

N = 100
g = watts_strogatz(N, 4, 0)


function network!(dx, x, p, t)
    for i in vertices(g)
        ### This may be repeated for several layers
        node_interactions = [] # tbd without allocations
        for j in neighbors(g,i)
            append!(node_interactions, interact(x[i], x[j]))
        end
        ### more inputs are possible
        node!(dx[i], x[i], p, t, Dict("Layer1" => aggregate(node_interactions)))
    end
end


x0    = randn(N)
p     = nothing
tspan = (0., 5.)

prob  = ODEProblem(network!, x0, tspan)
sol   = solve(prob, Tsit5())

plot(sol)


for e in edges(g)
    print(src(e))
end



edges(g)

?edges(g)

using BenchmarkTools

Earr = collect(edges(g))
P = 1=>2

@btime for e in $Parr e[1] end

@btime for e in $Earr src(e) end


Parr = Pair.(edges(g))

smat = sparse(first.(Parr), last.(Parr), randn(ne(g)))

dmat = Array(smat)

@btime sum( $smat[1,:])

@btime sum($dmat[1,:])

@btime findfirst(x -> x == (1=>2), $Parr)

@btime $smat[1,1] = 1
@btime $dmat[1,1] = 1



Parr = first.(Pair.(edges(g)))

create_idxs(last.(Parr), nv(g))
