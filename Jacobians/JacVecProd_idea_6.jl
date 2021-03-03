#nd_jac = network_dynamics(nd_jac_vertex, nd_jac_edge, g)

using ForwardDiff
#using DiffEqOperators
using NetworkDynamics


## TO DO: in Utilities.jl
function internal_jacobian(v, p, t) # jeder Vertex hat den gleiche Formel für den internen Jacobian, wichtig ist nur welches v dann reingegen wird
    # a = v (Werte der vertices, für welche der Jacobian berechnet wird) - Array muss 1 x Anzahl der Vertex Variablen (x_i)
    # p = p, t = t
    # dx = random computer Werte - Array muss 1 x Anzahl der dx_i sein oder umgekehrt

    J = zeros(m, n) # Zeilenanzahl - m: dim der DGL (dx_i), Spaltenanzahl - n: dim der vertices (x_i)
    dx = Array{Float64,2}(undef, m, 1) # dx = Array{Float64,2}(undef, 1, m)
    cfg = ForwardDiff.JacobianConfig((y, x) -> f(y, x, edges, p, t), dx, v)
    ForwardDiff.jacobian!(J, (y, x) -> f(y, x, edges, p, t), dx, v, cfg, Val{false}())
    return J
end



# TO DO
## Ableitung nach source node
function source_jacobian(Jac::NDJacVecOperator, v_s, v_d, p, t)
    # a = v (Werte der vertices, für welche der Jacobian berechnet wird) - Array muss 1 x Anzahl der Vertex Variablen (e_i)
    # p = p, t = t
    # dx = random computer Werte - Array muss 1 x Anzahl der de_i sein oder umgekehrt

    gs = Jac.graph_structure
    checkbounds_p(p, gs.num_v, gs.num_e)
    gd = prep_gd(x, x, Jac.graph_data, Jac.graph_structure)

    m = gs.v_dims
    n = gs.v_variables # Anzahl der Vertex-Variablen

    J = zeros(m, n) # Zeilenanzahl - m: dim der DGL der edges (de_i), Spaltenanzahl - n: dim der edges (Variablenanzahl) (e_i)
    de = Array{Float64,2}(undef, m, 1) # de = Array{Float64,2}(undef, 1, m)
    cfg = ForwardDiff.JacobianConfig((y, x) -> f(y, x, v_d, p, t), de, v_s) # y, x = e, v_s
    ForwardDiff.jacobian!(J, (y, x) -> f(y, x, p, t), de, v_s, cfg, Val{false}())
    return J
end

# TO DO
# Ableitung nach destination node
function dest_jacobian(v_s, v_d, p, t)
    # a = v (Werte der vertices, für welche der Jacobian berechnet wird) - Array muss 1 x Anzahl der Vertex Variablen (e_i)
    # p = p, t = t
    # dx = random computer Werte - Array muss 1 x Anzahl der de_i sein oder umgekehrt
    J = zeros(m, n) # Zeilenanzahl - m: dim der DGL der edges (de_i), Spaltenanzahl - n: dim der edges (Variablenanzahl) (e_i)
    de = Array{Float64,2}(undef, m, 1) # de = Array{Float64,2}(undef, 1, m)
    cfg = ForwardDiff.JacobianConfig((y, x) -> f(y, v_s, x, p, t), de, v_d)
    ForwardDiff.jacobian!(J, (y, x) -> f(y, v_s, x, p, t), de, v_d, cfg, Val{false}())
    return J
end




### TO DO: wo werden diese Fkt integriert?
# VertexFunction and EdgeFunction will have new fields for Jacobians.
function VertexJacobian!(J, v, p, t)
   J = internal_jacobian(v, p, t)
   # get_vertex_jacobian(gd, i) = something(v,p,t)
   # gd.v_jac_array[i] = something(v,p,t)
end

function EdgeJacobian!(J_s, J_d, v_s, v_d, p, t)
   J_s = source_jacobian(v_s, v_d, p, t)
   J_d = dest_jacobian(v_s, v_d, p, t)
end


struct NDJacVecOperator{uType, P, tType,G,GraphStruct,GD,Bool}# <: DiffEqBase.AbstractDiffEqLinearOperator{T}
    x::uType # Punkt u, für welchen das jvp berechnet wird
    p::P
    t::tType
    graph::G
    graph_structure::GraphStruct
    graph_data::GD
    parallel::Bool
end

### TO DO: in NetworkStructures.jl
# field names von graph data
## ?? jacobians müssen nicht quadratisch sein
v_jac_array = [Array{Float64,2}(undef, dim, dim) for dim in v_dims]
# For edges there is another array layer in the buffer since each edge has two Jacobians
e_jac_array = [[zeros(dim,srcdim), zeros(dim,dstdim)] for (dim, srcdim, dstdim) in zip(e_dims, v_src_dims, v_dst_dims)]

e_jac_product = [zeros(e_dims) , zeros(e_dims)]

z = [zeros(num_v,v_dims)]
#e_jac_array = [zeros(3,3), zeros(3,3)]

#e_jac_array = [[1 3 5], [2 4 6]]
#i=1
#for i in 1:3
#    println(e_jac_array[1][i])
    #e_jac_array[2][i]
#end

### TO DO: in NetworkStructures.jl

@inline function get_src_edge_jacobian(gd::GraphData, i::Int) = gd.e_jac_array[i][1] ## ??? Syntax [1][i]

@inline function get_dst_edge_jacobian(gd::GraphData, i::Int) = gd.e_jac_array[i][2]

@inline function get_vertex_jacobian(gd::GraphData, i::Int) = gd.v_jac_array[i]

#@inline function get_z_direction(gd::GraphData, i::Int) = gd.z[i]

# z wird field von VertexFunction, aber auch von GraphData, um es aufzurufen
@inline function get_z_direction(gd::GraphData) = gd.z

#nd_jac_vertex= ODEVertex(f! = vertexfunc!, dim = 1, jacobian=true, z)

#@inline function get_z_direction(vf::VertexFunction, i::Int) = vf.z[i]

#@inline function get_z_direction(vf::VertexFunction) = vf.z

## default case: z = x

function update_coefficients!(Jac::NDJacVecOperator, x, p, t)

    gs = Jac.graph_structure
    checkbounds_p(p, gs.num_v, gs.num_e)
    gd = prep_gd(x, x, Jac.graph_data, Jac.graph_structure)

    for i in 1:gs.num_v
        vertices![i].VertexJacobian!(
          get_vertex_jacobian(gd, i),
          get_vertex(gd, i),
          p_v_idx(p, i),
          t)
    end

    for i in 1:gs.num_e
          edges![i].EdgeJacobian!(
              get_src_edge_jacobian(gd, i),
              get_dst_edge_jacobian(gd, i),
              get_src_vertex(gd, i),
              get_dst_vertex(gd, i),
              p_e_idx(p, i),
              t)
      end

    Jac.x = x
    Jac.p = p
    Jac.t = t
end


function (Jac::NDJacVecOperator)(x, p, t)
    update_coefficients!(Jac, x, p, t)
    Jac*x
end

function (Jac::NDJacVecOperator)(dx,x,p,t::Number)
    update_coefficients!(Jac,x,p,t)
    mul!(dx,Jac,x)
end

Base.:*(L::NDJacVecOperator,z::AbstractVector) = jac_vec_prod(Jac,z)

function LinearAlgebra.mul!(dx::AbstractVector,Jac::NDJacVecOperator,z::AbstractVector)
    let p=Jac.p,t=Jac.t
        jac_vec_prod!(dx,Jac)
    end
    jac_vec_prod!(dx,Jac,z)
end


# loops als Extrafunktion schreiben -> Operator überladen
# TO DO
# bei Chris passiert in der nicht mutating fct das gleiche wie in der anderen, nur dass dx nicht verändert wird
# muss nur x,z  bzw. dx, x, zübergeben werden?
function jac_vec_prod(Jac::NDJacVecOperator,z)

    gs = Jac.graph_structure
    checkbounds_p(p, gs.num_v, gs.num_e)
    gd = prep_gd(x, x, Jac.graph_data, Jac.graph_structure)

    #z = get_z_direction(gd)

    for i in 1:gs.num_e
        e_jac_product[i] .= get_src_edge_jacobian(gd, i) * view(z,get_src_indices(i)) + get_dst_edge_jacobian(gd, i) * view(z,get_dst_indices(i))
    end
    ## neues array erstellen und return
    for i in 1:gs.num_v
        view(dx, get_vertex_indices(i)) = get_vertex_jacobian(gd,i) * view(z, get_vertex_indices(i))
        #dx .+= sum(e_jac_product[i]) ????
    end

end

function jac_vec_prod!(dx,Jac::NDJacVecOperator,z)

    gs = Jac.graph_structure
    checkbounds_p(p, gs.num_v, gs.num_e)
    gd = prep_gd(x, x, Jac.graph_data, Jac.graph_structure)

    z = get_z_direction(gd)

    for i in 1:gs.num_e
        e_jac_product[i] .= get_src_edge_jacobian(gd, i) * view(z,get_src_indices(i)) + get_dst_edge_jacobian(gd, i) * view(z,get_dst_indices(i))
    end

    for i in 1:gs.num_v
        view(dx, get_vertex_indices(i)) = get_vertex_jacobian(gd,i) * view(z, get_vertex_indices(i))
        dx .+= sum(e_jac_product[i])
        #view(dx, get_vertex_indices(i)) .+= sum(e_jac_product[i]) ?
    end
end


#nd_jac(u) = JacVecProductOperator am Punkt u

#nd_jac(u)(z) = JacVecProduct Jac(u) * z
