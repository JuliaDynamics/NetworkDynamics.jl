using DifferentialEquations
using Reexport
using LinearAlgebra
using NetworkDynamics
using LightGraphs
using DiffEqBase
import DiffEqBase.update_coefficients!


include("/Users/lalalulu/Desktop/PIK/github/NetworkDynamics.jl/Jacobians/jac_structs.jl")
@reexport using .jac_structs

#include("/Users/lalalulu/Desktop/PIK/github/NetworkDynamics.jl/src/nd_ODE_Static.jl")
#@reexport using .nd_ODE_Static_mod

N = 4
k = 2
g = barabasi_albert(N, k)

function diffusionedge!(e, v_s, v_d, p, t)
    # usually e, v_s, v_d are arrays, hence we use the broadcasting operator .
    #e .= v_s .- v_d
    e[1] = v_s[1] - v_d[1]
    nothing
end

function diffusionvertex!(dv, v, edges, p, t)
    # usually dv, v, edges are arrays, hence we use the broadcasting operator .
    dv[1] = 0.
    #dv .= 0.
    dv[2] = 1.
    for e in edges
        dv[1] += e[1]
    end
    nothing
end

nd_diffusion_vertex = ODEVertex(f! = diffusionvertex!, dim = 2)

nd_diffusion_edge = StaticEdge(f! = diffusionedge!, dim = 1)

nd = network_dynamics(nd_diffusion_vertex, nd_diffusion_edge, g)

x0 = randn(2N) # random initial conditions
ode_prob = ODEProblem(nd, x0, (0., 4.))
sol = solve(ode_prob, Tsit5());

# collecting the graph infos

graph_structure_ = nd.f.graph_structure
graph_data_ = nd.f.graph_data
#graph_ = nd.f.graph

v_dims = nd.f.graph_structure.v_dims
e_dims = nd.f.graph_structure.e_dims
num_e = nd.f.graph_structure.num_e

# build the jacobian graph data

struct JacGraphData
    v_jac_array::Array{Array{Float64, 2}, 1}
    e_jac_array::Array{Array{Array{Float64, 2}, 1}, 1}
    e_jac_product::Array{Float64, 2}
end

v_jac_array = [Array{Float64,2}(undef, dim, dim) for dim in v_dims]
e_jac_array = [[zeros(dim, srcdim), zeros(dim, dstdim)] for (dim, srcdim, dstdim) in zip(e_dims, v_dims, v_dims)] # homogene Netzwerke: v_src_dim = v_dst_dim = v_dim
e_jac_product =  zeros(e_dims[1], num_e) # Annahme: homogene edges

jac_graph_data_object = JacGraphData(v_jac_array, e_jac_array, e_jac_product)

x = similar(zeros(1), sum(v_dims))
p = nothing
t = 0.0
parallel = false

#### NDJacVecOperator

mutable struct NDJacVecOperator{T, uType, tType, G, GD, JGD} <: DiffEqBase.AbstractDiffEqLinearOperator{T} # mutable da x, p, t geupdated werden
    x::uType
    p
    t::tType
    graph::G
    graph_structure::GraphStruct
    graph_data::GD
    jac_graph_data::JGD
    parallel::Bool

    function NDJacVecOperator{T}(x, p, t, graph, graph_structure, graph_data, jac_graph_data, parallel) where T
        new{T,typeof(x),typeof(t),typeof(graph),typeof(graph_data),typeof(jac_graph_data)}(x, p, t, graph, graph_structure, graph_data, jac_graph_data, parallel)
    end

    function NDJacVecOperator(x, p, t, graph, graph_structure, graph_data, jac_graph_data, parallel)
        NDJacVecOperator{eltype(x)}(x, p, t, graph, graph_structure, graph_data, jac_graph_data, parallel)
    end
end

NDJacVecOperator_object = NDJacVecOperator(x, p, t, g, graph_structure_, graph_data_, jac_graph_data_object, parallel)


### get functions for update_coefficients

@inline get_src_edge_jacobian(jgd::JacGraphData, i::Int) = jgd.e_jac_array[i][1]
@inline get_dst_edge_jacobian(jgd::JacGraphData, i::Int) = jgd.e_jac_array[i][2]
@inline get_vertex_jacobian(jgd::JacGraphData, i::Int) = jgd.v_jac_array[i]

### Vertex, Edge functions for update_coefficients

function vertex_jacobian!(J::AbstractMatrix, v, p, t)
    #J = internal_jacobian(J, v, p, t)
    J[1, 1] = 1.0
    J[1, 2] = 0.0
    J[2, 1] = 1.0
    J[2, 2] = 0.0
end

function edge_jacobian!(J_s::AbstractMatrix, J_d::AbstractMatrix, v_s, v_d, p, t)
   #J_s = source_jacobian(v_s, v_d, p, t)
   #J_d = dest_jacobian(v_s, v_d, p, t)
   J_s[1, 1] = 1.0
   J_s[1, 2] = 0.0
   J_d[1, 1] = 1.0
   J_d[1, 2] = 0.0
end

function update_coefficients!(Jac::NDJacVecOperator, x, p, t)

    gs = Jac.graph_structure
    checkbounds_p(p, gs.num_v, gs.num_e)
    gd = prep_gd(x, x, Jac.graph_data, Jac.graph_structure)
    jgd = Jac.jac_graph_data

    for i in 1:gs.num_v
        #vertices![i].VertexJacobian!(
        #nd_diffusion_vertex[i].VertexJacobian!(
        nd_diffusion_vertex.vertex_jacobian!(
          get_vertex_jacobian(jgd, i),
          get_vertex(gd, i),
          p_v_idx(p, i),
          t)
    end

    for i in 1:gs.num_e
          #edges![i].EdgeJacobian!(
          #nd_diffusion[i].EdgeJacobian!(
          nd_diffusion.edge_jacobian!(
              get_src_edge_jacobian(jgd, i),
              get_dst_edge_jacobian(jgd, i),
              get_src_vertex(gd, i),
              get_dst_vertex(gd, i),
              p_e_idx(p, i),
              t)
      end

    Jac.x = x
    Jac.p = p
    Jac.t = t
end

# functions needed

# If the type of dx and x match, we swap out v_array for x
@inline function prep_gd(dx::AbstractArray{T}, x::AbstractArray{T}, gd::GraphData{GDB, T, T}, gs) where {GDB, T}
    # We construct views into an Array of size dim_v, so when we swap the
    # underlying Array it has to have to same size
    if size(x) == (gs.dim_v,)
         swap_v_array!(gd, x)
         return gd
    else
         error("Size of x does not match the dimension of the system.")
    end
end

# If the type of dx and x do not match, we swap initialize a new GraphData object
# that is based on the type of dx for the edge buffer.
# Some solvers take the derivative with respect to time, thus x will not be dual
# but dx will be, leading to errors otherwise
@inline function prep_gd(dx, x, gd, gs)
    # Type mismatch
    if size(x) == (gs.dim_v,)
        e_array = similar(dx, gs.dim_e)
        return GraphData(x, e_array, gs)
    else
        error("Size of x does not match the dimension of the system.")
   end
end

# function call update_coefficients!

x_test = randn(2*N)
p_test = nothing
t_test = 1.0

update_coefficients!(NDJacVecOperator_object, x_test, p_test, t_test)

# callable structs

function (Jac::NDJacVecOperator)(x, p, t) # auch Number bei t?
    update_coefficients!(Jac, x, p, t)
    Jac*x
end

function (Jac::NDJacVecOperator)(dx, x, p, t::Number)
    update_coefficients!(Jac, x, p, t)
    mul!(dx, Jac, x)
end

Base.:*(Jac::NDJacVecOperator, z::AbstractVector) = jac_vec_prod(Jac, z)

function LinearAlgebra.mul!(dx::AbstractVector, Jac::NDJacVecOperator, z::AbstractVector)
    jac_vec_prod!(dx, Jac, z)
end

# final functions for NDJacVecOperator

function jac_vec_prod(Jac::NDJacVecOperator, z)

    gs = Jac.graph_structure
    checkbounds_p(p, gs.num_v, gs.num_e)
    gd = prep_gd(x, x, Jac.graph_data, Jac.graph_structure)
    jgd = Jac.jac_graph_data

    for i in 1:gs.num_e
        e_jac_product[i] .= get_src_edge_jacobian(jgd, i) * view(z, get_src_indices(i)) + get_dst_edge_jacobian(jgd, i) * view(z, get_dst_indices(i))
    end
    ## neues array wird erstellt und returned

    dx = zeros(gs.v_dims, gs.num_v)

    for i in 1:gs.num_v
        view(dx, gs.v_idx[i]) .= get_vertex_jacobian(gd, i) * view(z, gs.v_idx[i])
        dx .+= sum(e_jac_product[i])
    end
    return dx
end

function jac_vec_prod!(dx, Jac::NDJacVecOperator, z)

    gs = Jac.graph_structure
    checkbounds_p(p, gs.num_v, gs.num_e)
    gd = prep_gd(x, x, Jac.graph_data, Jac.graph_structure)

    for i in 1:gs.num_e
        e_jac_product[i] .= get_src_edge_jacobian(gd, i) * view(z, get_src_indices(i)) + get_dst_edge_jacobian(gd, i) * view(z, get_dst_indices(i))
    end

    for i in 1:gs.num_v
        view(dx, get_vertex_indices(i)) .= get_vertex_jacobian(gd, i) * view(z, get_vertex_indices(i))
        dx .+= sum(e_jac_product[i])
    end
end
