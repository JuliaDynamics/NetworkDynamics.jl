using NetworkDynamics
using OrdinaryDiffEq
using LightGraphs
using Plots

N = 5
k = 2
g = barabasi_albert(N, k)
t_max = 10
x0 = [rand() for i in 1:N]
tmax = 10
# parameters
omegas = rand()*2π
coupling_constant = 3.0/(N-1)

struct jac_vertex
    Jac_intern
end

Jac_intern(v) = 0.0

function (jv::jac_vertex)(dv, v, edges, p, t) # dv= output; p= z
    dv .= jv.Jac_intern(v) * p
    for e in edges
        dv .+= e
    end
    nothing
end

struct jac_edge_struct
    Jac_src
    Jac_dst
end

Jac_src(v_s, v_d) = coupling_constant * cos.(v_s-v_d)
Jac_dst(v_s, v_d) = coupling_constant * -cos.(v_s-v_d)

function (je::jac_edge_struct)(e, v_s, v_d, p, t) # p = [z_s z_d]
    e .= je.Jac_src(v_s,v_d) * p[1] + je.Jac_dst(v_s,v_d) * p[2]
    nothing
end



vertexfunc! = jac_vertex(Jac_intern)
edgefunc! = jac_edge_struct(Jac_src, Jac_dst)

nd_jac_edge = StaticEdge(f! = edgefunc!, dim = 1)
nd_jac_vertex= ODEVertex(f! = vertexfunc!, dim = 1)

### directional derivative at location of vertex 1
z = zeros(Float64,N)
z[1] = 1

### get the edge indices for vector z in the edge function
edge_list = collect(edges(g))
edge_list = edge_list[:,:]

edge_z = zeros(Float64, ne(g), 2)

for i in 1:ne(g)
    edge_z[i,:] = [z[src(edge_list[i])] z[dst(edge_list[i])]]
end

### put the z vectors into the parameters
p = (z, edge_z')

nd_jac = network_dynamics(nd_jac_vertex, nd_jac_edge, g)


struct JacVecOperator
    punkt::GraphDataBuffer
    vertexlist = [v.Jacintern for v in verties]
    edgelist = ...
end

function update_coefficients!(Jac, u)
    for vertices
        update Jac_intern
    end
    for edges
        update Jac_src
        update Jac_dst
    end
    Jac.punkt = u
end

function Jac_u(out, z)
    for edges
        e_s = Jac_src * z_src
        e_d = Jac_dst * z_dst
    for vertices
        out = Jac_intern_u * z
    end
end


nd_jac(u) = JacVecProductOperator am Punkt u

nd_jac(u)(z) = JacVecProduct Jac(u) * z


prob_jac = ODEProblem(nd_jac, x0, (0.,tmax), p)

sol_jac = solve(prob_jac, Tsit5())

plot(sol_jac)

println(edge_z'[1,:])
