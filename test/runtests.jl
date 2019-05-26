using Test
using NetworkDynamics
using LightGraphs
using LinearAlgebra
using DifferentialEquations
using Plots

g = barabasi_albert(10,5)
g = SimpleGraph(10)
for i in 1:9
    add_edge!(g,i,i+1)
    add_edge!(g,i, mod(i+5,10))
end

L = laplacian_matrix(g)
function diff_net!(dx, x, p, t)
    dx .= - L * x
end

function vertex!(dv, v, e_s, e_d, p, t)
    dv .= 0
    for e in e_s
        dv .-= e
    end
    for e in e_d
        dv .+= e
    end
    nothing
end

function ddevertex!(dv,v,h_v,e_s,e_d,h_s,h_d,p,t)
    dv .= 0
    for h in h_s
        dv .-= h
    end
    for h in h_d
        dv .+= h
    end
end

odeedge! = (dl,l,v_s,v_d,p,t) -> dl .= 1000*(v_s - v_d - l)
staticedge! = (l,v_s,v_d,p,t) -> l .= v_s - v_d
ddeedge! = (de,e,h_e,v_s,v_d,h_s,h_d,p,t) -> de .= 1000*(v_s - v_d - e)

odevertices = [ODEVertex(f! = vertex!, dim = 1, mass_matrix = I) for v in vertices(g)]
odeedges = [ODEEdge(f! = odeedge!, dim = 1) for e in edges(g)]
staticedges = [StaticEdge(f! = staticedge!, dim = 1) for e in edges(g)]
# ddevertices = [DDEVertex(f! = ddevertex!, dim = 1, mass_matrix = I) for v in vertices(g)]
# ddeedges = [DDEEdge(f! = ddeedge!, dim = 1) for e in edges(g)]


nd_static = network_dynamics(odevertices, staticedges, g)
nd_ode = network_dynamics(odevertices, odeedges, g)
# nd_dde = network_dynamics(ddevertices,ddeedges,g)

@testset "Network dynamics function" begin
    for i in 1:10
        x = randn(10)
        dx1 = similar(x)
        dx2 = similar(x)
        nd_static(dx2,x,nothing,0.)
        diff_net!(dx1,x,nothing,0.)
        @test isapprox(dx1,dx2)
    end
end


x0 = rand(nv(g)+ne(g))

problem_static_1 = ODEProblem(nd_static, x0[1:nv(g)], (0.,2.))
problem_static_2 = ODEProblem(diff_net!,x0[1:nv(g)],(0.,2.))
solution_static_1 = solve(problem_static_1)
solution_static_2 = solve(problem_static_2)
@testset "Static solution" begin
    @test isapprox(solution_static_1.u, solution_static_2.u)
end

h0(p,t; idxs = 1:nv(g)+ne(g)) = x0

# plot(solution_static_1, legend = false, vars = 1:10)

problem_ode = ODEProblem(nd_ode, x0, (0.,2.))
solution_ode = solve(problem_ode)
# plot(solution_ode,legend = false, vars= 1:10)

# problem_dde = DDEProblem(nd_dde, x0, h0, (0.,2.))
# solution_dde = solve(problem_dde)
# plot(solution_dde,legend = false, vars = 1:10)

@testset "Check if static solution == ode solution" begin
    @test isapprox(solution_static_1[end][1:10], solution_ode[end][1:10], atol=0.001)
end

# for i in 1:10
#     @test isapprox(solution_static_1[end][i], solution_dde[end][i], atol=0.001)
# end
# end
