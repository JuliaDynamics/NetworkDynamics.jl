using NetworkDynamics
using ForwardDiff
using LightGraphs


include()

N = 5
k = 2
g = barabasi_albert(N, k)
x0 = [rand() for i in 1:N]
tspan = (0., 10.)
# parameters
omegas = rand()*2π
coupling_constant = 3.0/(N-1)
p = (nothing, nothing)

t_max=10.
edges1=0.

function diffusionedge!(e, v_s, v_d, p, t)
    # usually e, v_s, v_d are arrays, hence we use the broadcasting operator .
    e .= v_s .- v_d
    nothing
end

function diffusionvertex!(dv, v, edges, p, t)
    # usually dv, v, edges are arrays, hence we use the broadcasting operator .
    dv = 0.
    for e in edges
        dv .+= e
    end
    nothing
end


function VertexJacobian!(J, v, p, t)
   J = internal_jacobian(v, p, t)
   # get_vertex_jacobian(gd, i) = something(v,p,t)
   # gd.v_jac_array[i] = something(v,p,t)
end


# signature of StaticVertex: (f!, dim, mass_matrix, sym, z, jac)
nd_jac_vertex = StaticVertex(f! = diffusionvertex!, dim = 1)#, jac=true
# signature of StaticEdge: (f!, dim, coupling, sym, jac)
nd_jac_edge = StaticEdge(f! = diffusionedge!, dim = 1)#, jac=true)

nd_jac = network_dynamics(nd_jac_vertex, nd_jac_edge, g, jac = true)



nd_jac_operator = NDJacVecOperator(x0,p,tspan,g) #gs, GD,bool?

for t in 1:tmax
    nd_jac_operator(x0,p,t) ## Jac*z
end

nd_jac_operator(dx,x0,p,tspan) ## J mul!(dx,Jac,z)

# jac = jacvecprod in source ordner
