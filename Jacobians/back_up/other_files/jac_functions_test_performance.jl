using DifferentialEquations
using Reexport
using LinearAlgebra
using NetworkDynamics
using LightGraphs
using DiffEqBase
import DiffEqBase.update_coefficients!

N = 4
k = 2
g = barabasi_albert(N, k)

function diffusionedge!(e, v_s, v_d, p, t)
    e[1] = v_s[1] - v_d[1]
    nothing
end

function diffusionvertex!(dv, v, edges, p, t)
    dv[1] = 0.
    dv[2] = 1.
    for e in edges
        dv[1] += e[1]
    end
    nothing
end

function jac_vertex!(J::AbstractMatrix, v, p, t)
    #J = internal_jacobian(J, v, p, t)
    J[1, 1] = 1.0
    J[1, 2] = 0.0
    J[2, 1] = 1.0
    J[2, 2] = 0.0
end

function jac_edge!(J_s::AbstractMatrix, J_d::AbstractMatrix, v_s, v_d, p, t)
   #J_s = source_jacobian(v_s, v_d, p, t)
   #J_d = dest_jacobian(v_s, v_d, p, t)
   J_s[1, 1] = 1.0
   J_s[1, 2] = 0.0
   J_d[1, 1] = 1.0
   J_d[1, 2] = 0.0
end



nd_diffusion_vertex = ODEVertex(f! = diffusionvertex!, dim = 2, vertex_jacobian! = jac_vertex!)

nd_diffusion_edge = StaticEdge(f! = diffusionedge!, dim = 1, edge_jacobian! = jac_edge!)

nd = network_dynamics(nd_diffusion_vertex, nd_diffusion_edge, g)

x0 = randn(2N) # random initial conditions
ode_prob = ODEProblem(nd, x0, (0., 4.))
sol = solve(ode_prob, Tsit5());

# collecting the graph infos

graph_structure_ = nd.f.graph_structure
graph_data_ = nd.f.graph_data

v_dims = nd.f.graph_structure.v_dims
e_dims = nd.f.graph_structure.e_dims
num_e = nd.f.graph_structure.num_e


vertices! = nd.f.vertices!
edges! = nd.f.edges!

# for vertices! and edges! we need to get the index
@inline Base.@propagate_inbounds function maybe_idx(p::T, i) where T <: AbstractArray
    p[i]
end

@inline function maybe_idx(p, i)
    p
end

# build the jacobian graph data

struct JacGraphData{JGDB}
    jgdb::JGDB
    v_jac_array::Array{Array{Float64, 2}, 1}
    e_jac_array::Array{Array{Array{Float64, 2}, 1}, 1}
    #e_jac_product::Array{Float64, 2}
end

mutable struct JacGraphDataBuffer{Tvj, Tej}
    v_Jac_array::Tvj
    e_Jac_array::Tej
    #e_Jac_product_array::Tep
end

function JacGraphData(v_Jac_array::Tvj, e_Jac_array::Tej, gs::GraphStruct) where {Tvj, Tej}
    jgdb = JacGraphDataBuffer{Tvj, Tej}(v_Jac_array, e_Jac_array)
    JGDB = typeof(jgdb)
    #v_jac = [Array{Float64,2}(undef, dim, dim) for dim in gs.v_dims]
    #e_jac = [[zeros(dim, srcdim), zeros(dim, dstdim)] for (dim, srcdim, dstdim) in zip(gs.e_dims, gs.v_dims, gs.v_dims)] # homogene Netzwerke: v_src_dim = v_dst_dim = v_dim
    #e_jac_product =  zeros(gs.e_dims[1], gs.num_e) # Annahme: homogene edges
    v_jac = [Array{Float64,2}(undef, dim, dim) for dim in v_dims]
    e_jac = [[zeros(dim, srcdim), zeros(dim, dstdim)] for (dim, srcdim, dstdim) in zip(e_dims, v_dims, v_dims)] # homogene Netzwerke: v_src_dim = v_dst_dim = v_dim
    #e_jac_product =  zeros(num_e, e_dims[1])

    JacGraphData{JGDB}(jgdb, v_jac, e_jac)
end

v_jac_array = [Array{Float64,2}(undef, dim, dim) for dim in v_dims]
e_jac_array = [[zeros(dim, srcdim), zeros(dim, dstdim)] for (dim, srcdim, dstdim) in zip(e_dims, v_dims, v_dims)] # homogene Netzwerke: v_src_dim = v_dst_dim = v_dim
#e_jac_product =[Matrix{Float64}(zeros(num_e, e_dims[1]))] # Annahme: homogene edges
#e_jac_product =  zeros(num_e, e_dims[1])
#e_jac_product = [Array{Float64,2}(dim, e_dim) for (dim, e_dim) in zip(e_dims[1], num_e)]
#e_jac_product = [[zeros(dim, e_dim)] for (dim, e_dim) in zip(e_dims[1], num_e)]

#jac_graph_data_object = JacGraphData(v_jac_array, e_jac_array, e_jac_product)
jac_graph_data_object = JacGraphData(v_jac_array, e_jac_array,graph_structure_)

x = similar(zeros(1), sum(v_dims))
p = nothing
t = 0.0
parallel = false

#### NDJacVecOperator

mutable struct NDJacVecOperator{T, T1, T2, uType, tType, G, GD, JGDB} <: DiffEqBase.AbstractDiffEqLinearOperator{T} # mutable da x, p, t geupdated werden
    vertices!::T1
    edges!::T2
    x::uType
    p
    t::tType
    graph::G
    graph_structure::GraphStruct
    graph_data::GD
    jac_graph_data::JGDB
    parallel::Bool

    function NDJacVecOperator{T}(vertices!, edges!, x, p, t, graph, graph_structure, graph_data, jac_graph_data, parallel) where T
        new{T, typeof(vertices!), typeof(edges!), typeof(x),typeof(t),typeof(graph),typeof(graph_data),typeof(jac_graph_data)}(vertices!, edges!, x, p, t, graph, graph_structure, graph_data, jac_graph_data, parallel)
    end

    function NDJacVecOperator(vertices!, edges!, x, p, t, graph, graph_structure, graph_data, jac_graph_data, parallel)
        NDJacVecOperator{eltype(x)}(vertices!, edges!, x, p, t, graph, graph_structure, graph_data, jac_graph_data, parallel)
    end
end

NDJacVecOperator_object = NDJacVecOperator(vertices!,edges!,x, p, t, g, graph_structure_, graph_data_, jac_graph_data_object, parallel)


### get functions for update_coefficients

@inline get_src_edge_jacobian(jgd::JacGraphData, i::Int) = jgd.e_jac_array[i][1]
@inline get_dst_edge_jacobian(jgd::JacGraphData, i::Int) = jgd.e_jac_array[i][2]
@inline get_vertex_jacobian(jgd::JacGraphData, i::Int) = jgd.v_jac_array[i]

### Vertex, Edge functions for update_coefficients

function update_coefficients!(Jac::NDJacVecOperator, x, p, t)

    gs = Jac.graph_structure
    checkbounds_p(p, gs.num_v, gs.num_e)
    gd = prep_gd(x, x, Jac.graph_data, Jac.graph_structure)
    jgd = Jac.jac_graph_data

    for i in 1:gs.num_v
        #vertices![i].vertex_jacobian!(
        maybe_idx(vertices!, i).vertex_jacobian!(
          get_vertex_jacobian(jgd, i),
          get_vertex(gd, i),
          p_v_idx(p, i),
          t)
    end

    for i in 1:gs.num_e
          #edges![i].edge_jacobian!(
          maybe_idx(edges!, i).edge_jacobian!(
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


# final functions for NDJacVecOperator

function e_jac_prod!(z, ejp, gs, jgd)
    for i in 1:gs.num_e
        ejp[i, :] .= get_src_edge_jacobian(jgd, i) * view(z, gs.s_e_idx[i]) + get_dst_edge_jacobian(jgd, i) * view(z, gs.d_e_idx[i])
    end
end;

function v_jac_prod!(dx, z, ejp, gs, jgd)
    for i in 1:gs.num_v
        view(dx, gs.v_idx[i]) .= get_vertex_jacobian(jgd, i) * view(z, gs.v_idx[i])
        view(dx, gs.v_idx[i]) .+= vec(sum(view(ejp, gs.d_v[i], :), dims = 1))
    end
end;


function jac_vec_prod(Jac::NDJacVecOperator, z)

    gs = Jac.graph_structure
    checkbounds_p(p, gs.num_v, gs.num_e)
    gd = prep_gd(x, x, Jac.graph_data, Jac.graph_structure)
    jgd = Jac.jac_graph_data

    ejp = zeros(gs.num_e, gs.e_dims[1])
    e_jac_prod!(z, ejp, gs, jgd)

    dx = zeros(gs.v_dims[1]*gs.num_v)
    v_jac_prod!(dx, z, ejp, gs, jgd)

    return dx

end

function jac_vec_prod!(dx, Jac::NDJacVecOperator, z)

    gs = Jac.graph_structure
    checkbounds_p(p, gs.num_v, gs.num_e)
    gd = prep_gd(x, x, Jac.graph_data, Jac.graph_structure)
    jgd = Jac.jac_graph_data

    e_jac_prod!(z, e_jac_product, gs, jgd)
    v_jac_prod!(dx, z, e_jac_product, gs, jgd)

end

function jac_vec_prod(Jac::NDJacVecOperator, z)

    gs = Jac.graph_structure
    checkbounds_p(p, gs.num_v, gs.num_e)
    gd = prep_gd(x, x, Jac.graph_data, Jac.graph_structure)
    jgd = Jac.jac_graph_data
    #e_jac_product = Jac.
    #println(typeof(e_jac_product))

    for i in 1:gs.num_e
        #e_jac_product[i] .= get_src_edge_jacobian(gd, i) * view(z, get_src_indices(i)) + get_dst_edge_jacobian(gd, i) * view(z, get_dst_indices(i))
        e_jac_product[i, :] .= get_src_edge_jacobian(jgd, i) * view(z, gs.s_e_idx[i]) + get_dst_edge_jacobian(jgd, i) * view(z, gs.d_e_idx[i])
    end

    ## neues array wird erstellt und returned
    dx = zeros(gs.v_dims[1]*gs.num_v)

    for i in 1:gs.num_v
        view(dx, gs.v_idx[i]) .= get_vertex_jacobian(jgd, i) * view(z, gs.v_idx[i])
        view(dx, gs.v_idx[i]) .+= vec(sum(view(e_jac_p, gs.d_v[i], :), dims = 1))
    end
    return dx
end

function jac_vec_prod!(dx, Jac::NDJacVecOperator, z)

    gs = Jac.graph_structure
    checkbounds_p(p, gs.num_v, gs.num_e)
    gd = prep_gd(x, x, Jac.graph_data, Jac.graph_structure)
    jgd = Jac.jac_graph_data
    display(@allocated e_jac_p = zeros(gs.num_e, gs.e_dims[1]))

    for i in 1:gs.num_e
        #e_jac_product[i, :] .= get_src_edge_jacobian(jgd, i) * view(z, gs.s_e_idx[i]) + get_dst_edge_jacobian(jgd, i) * view(z, gs.d_e_idx[i])
        e_jac_p[i, :] .= get_src_edge_jacobian(jgd, i) * view(z, gs.s_e_idx[i]) + get_dst_edge_jacobian(jgd, i) * view(z, gs.d_e_idx[i])
    end

    for i in 1:gs.num_v
        view(dx, gs.v_idx[i]) .= get_vertex_jacobian(jgd, i) * view(z, gs.v_idx[i])
        view(dx, gs.v_idx[i]) .+= vec(sum(view(e_jac_p, gs.d_v[i], :), dims = 1)) # display(@allocated)
        #println(sum(view(e_jac_product, gs.d_v[i], :), dims = 1))
    end
end


# functions for callable structs

Base.:*(Jac::NDJacVecOperator, z::AbstractVector) = jac_vec_prod(Jac, z)

function LinearAlgebra.mul!(dx::AbstractVector, Jac::NDJacVecOperator, z::AbstractVector)
    jac_vec_prod!(dx, Jac, z)
end

# callable structs

function (Jac::NDJacVecOperator)(x, p, t) # auch Number bei t?
    update_coefficients!(Jac, x, p, t)
    Jac*x
end

function (Jac::NDJacVecOperator)(dx, x, p, t::Number)
    update_coefficients!(Jac, x, p, t)
    mul!(dx, Jac, x)
end



#x_test = randn(N)
x_test = ones(2N)
p_test = nothing
t_test = 1.0
dx_test = ones(2N)

@time update_coefficients!(NDJacVecOperator_object, x_test, p_test, t_test)
call_update_coefficients! = update_coefficients!(NDJacVecOperator_object, x_test, p_test, t_test)


call_callable_struct_1 = NDJacVecOperator_object(x_test, p_test, t_test)

call_callable_struct_2 = NDJacVecOperator_object(dx_test, x_test, p_test, t_test)

@time NDJacVecOperator_object(x_test, p_test, t_test)
@time NDJacVecOperator_object(dx_test, x_test, p_test, t_test)

@allocated NDJacVecOperator_object(x_test, p_test, t_test)
@allocated NDJacVecOperator_object(dx_test, x_test, p_test, t_test)
