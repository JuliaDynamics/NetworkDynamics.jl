using NetworkDynamics
using LightGraphs
using LinearAlgebra
using SparseArrays
using DifferentialEquations
using Plots

struct SolAnalytic
    L
end
function (sa::SolAnalytic)(x0, p, t)
    exp(- t * sa.L) * x0
end

g = barabasi_albert(10,5)

L = laplacian_matrix(g)

sol_analytic = SolAnalytic(Symmetric(Array(L)))

diff_network_L = ODEFunction((dx, x, p, t) -> dx .= - L * x, analytic=sol_analytic)

x0 = rand(10)

prob_L = ODEProblem(diff_network_L,x0,(0.,5.))
sol_L = solve(prob_L)

println(sol_L.errors)

plot(sol_L)

sol_ana = zeros(length(sol_L.t), length(x0))
for i in 1:length(sol_L.t)
    sol_ana[i, :] .= sol_analytic(x0, nothing, sol_L.t[i])
end
plot(sol_L.t, sol_ana)

#Now for NetworkDynamics

diffusion_edge! = (e,v_s,v_d,p,t) -> @. e = v_s - v_d

function diffusion_vertex!(dv, v, e_s, e_d, p, t)
    dv .= edge_sum(e_s, e_d) # Oriented sum of the incoming and outgoing edges
end

odevertex = ODEVertex(f! = diffusion_vertex!,dim = 1)
staticedge = StaticEdge(f! = diffusion_edge!, dim = 1)
odeedge = ODEEdge(staticedge) # We promote the static edge to an ODEEdge artifically

vertex_list = [odevertex for v in vertices(g)]
edge_list = [staticedge for e in edges(g)]

diff_network_st = network_dynamics(vertex_list,edge_list,g)

x0 = rand(10)

prob_st = ODEProblem(diff_network_st,x0,(0.,5.))

sol_st = solve(prob_st)

plot(sol_st)


struct Bla{T}
    x::T
end

f = () -> 2.
f2 = () -> 1.
b = Bla(f)

arr_1 = [Bla(f), Bla(f2)]
arr_2 = [Bla(f), Bla(f)]
arr_3 = [f, f]
arr_5 = Bla([f, f2])
arr_4 = Bla(arr_3)
tup_1 = (f, f2)
tup_2 = (f, f)

typeof(arr_1)
typeof(arr_2)

t1 = Tuple(arr_1)

function blaing(b)
    x = b.x
    count = 0.
    for i in 1:10^7
        count += x[1]()
    end
    count
end
function blaing2(b)
    count = 0.
    for i in 1:10^7
        count += b.x[1]()
    end
    count
end
function blaing3(arr)
    count = 0.
    for i in 1:10^7
        count += arr[1].x()
    end
    count
end

function blaing4(arr)
    f_t = arr[1].x
    count = 0.
    for i in 1:10^7
        count += f_t()
    end
    count
end

function blaing5_i(f_t)
    count = 0.
    for i in 1:10^7
        count += f_t()
    end
    count
end

function blaing5(arr)
    f_t = arr[1].x
    blaing5_i(f_t)
end

function tupping(tup)
    count = 0.
    for i in 1:10^7
        count += tup[1]()
    end
    count
end

function tupbla(blatup)
    count = 0.
    for i in 1:10^7
        count += blatup.x[1]()
    end
    count
end

println("Start")
@time for i in 1:10^7 f() end
@time for i in 1:10^7 arr_1[1].x() end
@time for i in 1:10^7 arr_2[1].x() end
@time for i in 1:10^7 arr_3[1]() end
@time for i in 1:10^7 tup_1[1]() end
@time for i in 1:10^7 tup_2[1]() end
@time for i in 1:10^7 t1[1].x() end
@time for i in 1:10^7 arr_4.x[1]() end
@time for i in 1:10^7 b.x() end
println("Start functions")
@time blaing5_i(f)
@time blaing(arr_4)
@time blaing(arr_5)
@time blaing2(arr_4)
@time blaing3(arr_1)
@time blaing3(arr_2)
@time blaing4(arr_1)
@time blaing4(arr_2)
@time blaing5(arr_1)
@time blaing5(arr_2)

@time tupping(tup_1)
@time tupping(tup_2) # This is the only inhomogeneous one that is fast!!

blatup_1 = Bla(tup_1)
blatup_2 = Bla(tup_2)

@time tupbla(blatup_1)
@time tupbla(blatup_2)

arr_1[1].x


import Base.getindex

struct GS_x_View{G}
    gs::G
    idx_offset::Int
end

function getindex(gs_x_view::GS_x_View, idx)
    gs_x_view.gs.x[idx + gs_x_view.idx_offset]
end

mutable struct GS{T}
    x::T
    gs_x_view::GS_x_View{GS{T}}
    function GS(arr::T) where T
        gs = new{T}(arr)
        gs.gs_x_view = GS_x_View{GS{T}}(gs, 1)
        gs
    end
end

arr1 = [1.,2.,3.]
arr2 = [1.,2.,3.] .+ 0.5

gs = GS(arr1)

gs.x[1] # 1.
gs.gs_x_view[1] # 2.

gs.x = arr2

gs.x[1] # 1.5
gs.gs_x_view[1] # 2.5
