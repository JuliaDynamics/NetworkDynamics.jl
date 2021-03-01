# Sketch
struct jac_vertex
    Jac_intern
end
function (jv::jac_vertex)(out, x, edges, z, t)
    out .= jv.Jac_intern * z
    for e in edges
        out .+= e
    end
end

struct jac_edge_struct
    Jac_src
    Jac_dst
end
function (je::jac_edge_struct)(e, v_s, v_d, p, t) # p = [z_s z_d]
    e =  je.Jac_src(v_s-v_d) * p[1] + je.Jac_dst * p[2]

end

Jac_src(v_s, v_d) = cos(v_s, v_d)
Jac_dst(v_s, v_d) = -cos(v_s, v_d)

nd_jac_edge = StaticEdge(jac_edge_struct(Jac_src, Jac_dst))

# diverting the parameter arguments
for i=1:ne(g)
    edgep[i] = [z[src(edge[i])] z[dst(edge[i])]]
end

p = (z, edgep)


# general idea
# 2 node kuramoto dynamics
dx1 = w1 + sin(x1 - x2)
dx2 = w2 + sin(x2 - x1)

# global JacVecProduct
J: (cos  -cos    z1   = cos z1 - cos z2
    -cos cos)    z2   = cos z2 - cos z1

# edge Jacobian
J_e: sin(x_src - xdst) ->     cos
                       ->    -cos





### Implementation strategy in sketched code


# 1. VertexFunction and EdgeFunction will have new fields for Jacobians.

function VertexJacobian!(J, v, p, t)
   J = something(v,p,t)
end

function EdgeJacobian!(J_s, J_d, v_s, v_d, p, t)
   J_s = something(v_s, v_d, p, t)
   J_d = something(v_s, v_d, p, t)
end

# 2.a There will be an internal buffer of Matrices that holds
# the current Jacobian at each point. (Implement as a new field of GraphData)

v_jac_array = [Array{Float64,2}(undef, dim, dim) for dim in v_dims]
# For edges there is another array layer in the buffer since each edge has two Jacobians
e_jac_array = [[zeros(dim,srcdim), zeros(dim,dstdim)] for (dim, srcdim, dstdim) in zip(e_dims, v_src_dims, v_dst_dims)]

# 2.b For convenient access we may define new Array types and/or new accessor functions to be used in the core loop, e.g.

@inline function get_src_edge_jacobian(gd::GraphData, i::Int) = gd.e_jac_array[i][1]



# 3. to update coeffcients we need a loop like the following
struct NDJacVec <: JacVecOperator
   ... # needs to have the right fields, similar to ODEStatic
end

function update_coefficients!(J::NDJacVec,u,p,t)
   # gs and gd are probably fields of J::NDJacVec
   for i in 1:gs.num_e
       edges![i].EdgeJacobian!(
           get_src_edge_jacobian(gd, i),
           get_dst_edge_jacobian(gd, i),
           get_src_vertex(gd, i),
           get_dst_vertex(gd, i),
           p_e_idx(p, i),
           t)
   end

   for i in 1:gs.num_v
       vertices![i].VertexJacobian!(
           get_vertex_jacobian(gd, i),
           get_vertex(gd, i),
           p_v_idx(p, i),
           t)
   end
end

# 4. The JacVecProduct might look like the following
# To begin with Jacobians will only be defined for networks with sum coupling and StaticEdges


function (J::NDJacVec)(output,u,p,t)
   update_coefficients!(J)
   # those loops follow the logic of our prototype
   for i in 1:gs.num_e
       # Im not sure if we want to reuse the internal edge or create a new data structure
       # get_edge accesses the internal edge array
       get_edge(gd,i) .= get_src_edge_jacobian(gd, i) * view(u,get_src_indices(i))
                         get_dst_edge_jacobian(gd, i) * view(u,get_dst_indices(i))
       # Careful to index into the correct subarray of u,
       # e.g. @inline get_src_indices(i) = gs.s_e_idx[i]
   end

   for i in 1:gs.num_v
       view(output, get_vertex_indices(i)) = get_vertex_jacobian(gd,i) * view(u, get_vertex_indices(i))
       # if we use the internal arrays get_dst_edges is probably what we need
       output .+= sum(get_dst_edges(i)) # here we hardcode the assumption of sum coupling
   end
end


# 5. The operator overloading of * and mul! follows a similar logic and should reuse the
# functions defined above
