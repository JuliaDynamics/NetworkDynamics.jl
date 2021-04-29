using DifferentialEquations
using Reexport
using LinearAlgebra
using NetworkDynamics
using LightGraphs
using DiffEqBase
import DiffEqBase.update_coefficients!
using Test

N = 4
k = 2
g = barabasi_albert(N, k)

function diffusionedge!(e, v_s, v_d, p, t)
    e[1] = v_s[1] - v_d[1]
    nothing
end

function diffusionvertex!(dv, v, edges, p, t)
    dv[1] = 0.
    for e in edges
        dv[1] += e[1]
    end
    nothing
end

function jac_vertex!(J::AbstractMatrix, v, p, t)
    #J = internal_jacobian(J, v, p, t)
    J[1, 1] = 0.0
end

function jac_edge!(J_s::AbstractMatrix, J_d::AbstractMatrix, v_s, v_d, p, t)
   J_s[1, 1] = 1.0
   #J_s[1, 2] = 0.0
   J_d[1, 1] = -1.0
   #J_d[1, 2] = 0.0
end

nd_jac_vertex = ODEVertex(f! = diffusionvertex!, dim = 1, vertex_jacobian! = jac_vertex!)

nd_jac_edge = StaticEdge(f! = diffusionedge!, dim = 1, edge_jacobian! = jac_edge!)

nd_jac = network_dynamics(nd_jac_vertex, nd_jac_edge, g, jac = true)

@test nd_jac_vertex isa ODEVertex
@test nd_jac_edge isa StaticEdge
@test nd_jac isa ODEFunction

@allocated nd_jac_vertex
@allocated nd_jac_edge
@allocated nd_jac
#@test nd_jac_vertex.vertex_jacobian! isa function jac_vertex!(J::AbstractMatrix, v, p, t)
#    J[1, 1] = 0.0
#end

nd_jac_vertex.vertex_jacobian!
x0 = ones(N)
#ode_prob = ODEProblem(nd, x0, (0., 4.))
#sol = solve(ode_prob, Tsit5());

# collecting the graph infos

graph_structure = nd_jac.f.graph_structure
graph_data = nd_jac.f.graph_data

v_dims = nd_jac.f.graph_structure.v_dims
e_dims = nd_jac.f.graph_structure.e_dims
num_e = nd_jac.f.graph_structure.num_e


vertices! = nd_jac.f.vertices!
edges! = nd_jac.f.edges!



# for vertices! and edges! we need to get the index
@inline Base.@propagate_inbounds function maybe_idx(p::T, i) where T <: AbstractArray
    p[i]
end

@inline function maybe_idx(p, i)
    p
end

struct JacGraphData1
    v_jac_array::Array{Array{Float64, 2}, 1}
    e_jac_array::Array{Array{Array{Float64, 2}, 1}, 1}
    e_jac_product::Array{Array{Float64, 1}, 1}
end

function JacGraphData1(v_jac_array, e_jac_array, e_jac_product_array, gs::GraphStruct)
    v_jac_array = [Array{Float64,2}(undef, dim, dim) for dim in gs.v_dims]
    e_jac_array = [[zeros(dim, srcdim), zeros(dim, dstdim)] for (dim, srcdim, dstdim) in zip(gs.e_dims, gs.v_dims, gs.v_dims)] # homogene Netzwerke: v_src_dim = v_dst_dim = v_dim
    e_jac_product = [zeros(gs.e_dims[1]) for i in 1:gs.num_e]
    JacGraphData1(v_jac_array, e_jac_array, e_jac_product)
end

v_jac_array = [Array{Float64,2}(undef, dim, dim) for dim in v_dims]
e_jac_array = [[zeros(dim, srcdim), zeros(dim, dstdim)] for (dim, srcdim, dstdim) in zip(e_dims, v_dims, v_dims)] # homogene Netzwerke: v_src_dim = v_dst_dim = v_dim
e_jac_product =  [zeros(e_dims[1]) for i in 1:num_e]

@allocated jac_graph_data_object = JacGraphData1(v_jac_array, e_jac_array, e_jac_product, graph_structure)

x = ones(N)
#x = similar(zeros(1), sum(v_dims))
p = nothing
t = 0.0
parallel = false

#### NDJacVecOperator

mutable struct NDJacVecOperator1{T, uType, tType, T1, T2, G, GD, JGD} <: DiffEqBase.AbstractDiffEqLinearOperator{T} # mutable da x, p, t geupdated werden
    x::uType
    p
    t::tType
    vertices!::T1
    edges!::T2
    graph::G
    graph_structure::GraphStruct
    graph_data::GD
    jac_graph_data::JGD
    parallel::Bool

    function NDJacVecOperator1{T}(x, p, t, vertices!, edges!, graph, graph_structure, graph_data, jac_graph_data, parallel) where T
        new{T, typeof(x), typeof(t), typeof(vertices!), typeof(edges!),typeof(graph), typeof(graph_data), typeof(jac_graph_data)}(x, p, t, vertices!, edges!, graph, graph_structure, graph_data, jac_graph_data, parallel)
    end

    function NDJacVecOperator1(x, p, t, vertices!, edges!, graph, graph_structure, graph_data, jac_graph_data, parallel)
        NDJacVecOperator1{eltype(x)}(x, p, t, vertices!, edges!, graph, graph_structure, graph_data, jac_graph_data, parallel)
    end
end

@time NDJacVecOp = NDJacVecOperator1(x, p, t, vertices!, edges!, g, graph_structure, graph_data, jac_graph_data_object, parallel) # 51 allocations

@test NDJacVecOp isa NDJacVecOperator1

### get functions for update_coefficients

@inline get_src_edge_jacobian(jgd::JacGraphData1, i::Int) = jgd.e_jac_array[i][1]
@inline get_dst_edge_jacobian(jgd::JacGraphData1, i::Int) = jgd.e_jac_array[i][2]
@inline get_vertex_jacobian(jgd::JacGraphData1, i::Int) = jgd.v_jac_array[i]

### Vertex, Edge functions for update_coefficients
## Test Block update coefficients

jgd1 = NDJacVecOp.jac_graph_data
i = 1
maybe_idx(NDJacVecOp.vertices!, i).vertex_jacobian!(
  get_vertex_jacobian(jgd1, i),
  get_vertex(graph_data, i),
  p_v_idx(p, i),
  t)

@time get_vertex_jacobian(jgd1, i)
@test get_vertex_jacobian(jgd1, i) == jgd1.v_jac_array[i]
# == 0 ?
get_vertex_jacobian(jgd1, i)
maybe_idx(NDJacVecOp.edges!, i).edge_jacobian!(
    get_src_edge_jacobian(jgd1, i),
    get_dst_edge_jacobian(jgd1, i),
    get_src_vertex(graph_data, i),
    get_dst_vertex(graph_data, i),
    p_e_idx(p, i),
    t)

@time get_src_edge_jacobian(jgd1, i)
@time get_dst_edge_jacobian(jgd1, i)

@test get_src_edge_jacobian(jgd1, i) == [1.0; 0.0]
@test get_dst_edge_jacobian(jgd1, i) == [-1.0; 0.0]

test = Matrix{Float64}(undef, 2, 1)
test[1, 1] = 1.0
test[2, 1] = 0.0

@test get_src_edge_jacobian(jgd1, i) == test

@test get_src_edge_jacobian(jgd1, i) == jgd1.e_jac_array[i][1]
@test get_dst_edge_jacobian(jgd1, i) == jgd1.e_jac_array[i][2]

function update_coefficients!(Jac::NDJacVecOperator, x, p, t)

    gs = Jac.graph_structure
    checkbounds_p(p, gs.num_v, gs.num_e)
    gd = prep_gd(x, x, Jac.graph_data, Jac.graph_structure)
    jgd = Jac.jac_graph_data

    #@test get_src_edge_jacobian(jgd, 1) ==

    for i in 1:gs.num_v
        maybe_idx(Jac.vertices!, i).vertex_jacobian!(
          get_vertex_jacobian(jgd, i),
          get_vertex(gd, i),
          p_v_idx(p, i),
          t)
    end

    for i in 1:gs.num_e
          maybe_idx(Jac.edges!, i).edge_jacobian!(
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

@inline function prep_gd(dx::AbstractArray{T}, x::AbstractArray{T}, gd::GraphData{GDB, T, T}, gs) where {GDB, T}
    if size(x) == (gs.dim_v,)
         swap_v_array!(gd, x)
         return gd
    else
         error("Size of x does not match the dimension of the system.")
    end
end

@inline function prep_gd(dx, x, gd, gs)
    # Type mismatch
    if size(x) == (gs.dim_v,)
        e_array = similar(dx, gs.dim_e)
        return GraphData(x, e_array, gs)
    else
        error("Size of x does not match the dimension of the system.")
   end
end

# final functions for NDJacVecOperator

function jac_vec_prod(Jac::NDJacVecOperator1, z)

    gs = Jac.graph_structure
    p = Jac.p
    x = Jac.x
    checkbounds_p(p, gs.num_v, gs.num_e)
    gd = prep_gd(x, x, Jac.graph_data, Jac.graph_structure)
    jgd = Jac.jac_graph_data

    src_edge_jac= get_src_edge_jacobian(jgd, 2)
    print(src_edge_jac)
    dst_edge_jac= get_dst_edge_jacobian(jgd, 2)
    print(src_edge_jac)
    vertex_jac= get_vertex_jacobian(jgd, 2)
    print(src_edge_jac)

    for i in 1:gs.num_e
        jgd.e_jac_product[i] .= get_src_edge_jacobian(jgd, i) * view(z, gs.s_e_idx[i]) + get_dst_edge_jacobian(jgd, i) * view(z, gs.d_e_idx[i])
    end

    ## neues array wird erstellt und returned
    dx = zeros(gs.v_dims[1], gs.num_v)

    for i in 1:gs.num_v
        #view(dx, get_vertex_indices(i)) .= get_vertex_jacobian(gd, i) * view(z, get_vertex_indices(i))
        view(dx, gs.v_idx[i]) .= get_vertex_jacobian(jgd, i) * view(z, gs.v_idx[i])
        dx .+= sum(jgd.e_jac_product[i])
    end
    return dx
end

### Test Block jac_vec_prod
gs = NDJacVecOp.graph_structure
p = NDJacVecOp.p
x = NDJacVecOp.x
# x = ones(N)
checkbounds_p(p, gs.num_v, gs.num_e)
gd = prep_gd(x, x, NDJacVecOp.graph_data, NDJacVecOp.graph_structure)
jgd1 = NDJacVecOp.jac_graph_data
z = x
dx = zeros(N)
#dx = [0.0,1.0,2.0,3.0]
## for i = 1
src_edge_jac = get_src_edge_jacobian(jgd1, 1)
dst_edge_jac = get_dst_edge_jacobian(jgd1, 1)
view(z, gs.v_idx[1])

@test src_edge_jac == [1.0 ; 0.0]
@test dst_edge_jac == [-1.0 ; 0.0]
@test view(z, gs.v_idx[1]) == [1.0]

@test get_src_edge_jacobian(jgd1, i) * view(z, gs.s_e_idx[i]) == [1.0 ; 0.0]
@test get_dst_edge_jacobian(jgd1, i) * view(z, gs.s_e_idx[i]) == [-1.0 ; 0.0]
# [1.0 ; 0.0] + [-1.0 ; 0.0] = [0.0 ; 0.0]
@test get_src_edge_jacobian(jgd1, i) * view(z, gs.s_e_idx[i]) + get_dst_edge_jacobian(jgd1, i) * view(z, gs.s_e_idx[i]) == [0.0 ; 0.0]

jgd1.e_jac_product isa Vector
jgd1.e_jac_product[i] .= get_src_edge_jacobian(jgd1, i) * view(z, gs.s_e_idx[i]) + get_dst_edge_jacobian(jgd1, i) * view(z, gs.d_e_idx[i])
@test jgd1.e_jac_product[i] == [0.0 ; 0.0]

view(dx, gs.v_idx[i])
@test view(dx, gs.v_idx[i]) == [0.0]

vertex_jac = get_vertex_jacobian(jgd1, i)
@test vertex_jac == [0.0]
@test view(z, gs.v_idx[i]) == [1.0]
@test get_vertex_jacobian(jgd1, i) * view(z, gs.v_idx[i]) == [0.0]
view(dx, gs.v_idx[i]) .= get_vertex_jacobian(jgd1, i) * view(z, gs.v_idx[i])
@test view(dx, gs.v_idx[i]) == [0.0]
dx .+= sum(jgd1.e_jac_product[i])
@test dx .+= sum(jgd1.e_jac_product[i]) == [0.0, 0.0, 0.0, 0.0]
#@test dx .+= sum(jgd1.e_jac_product[i]) == [0.0, 1.0, 2.0, 3.0]
function jac_vec_prod!(dx, Jac::NDJacVecOperator1, z)

    gs = Jac.graph_structure
    p = Jac.p
    x = Jac.x
    checkbounds_p(p, gs.num_v, gs.num_e)
    gd = prep_gd(x, x, Jac.graph_data, Jac.graph_structure)
    jgd = Jac.jac_graph_data

    for i in 1:gs.num_e
        jgd.e_jac_product[i] .= get_src_edge_jacobian(jgd, i) * view(z, gs.s_e_idx[i]) + get_dst_edge_jacobian(jgd, i) * view(z, gs.d_e_idx[i])
    end

    for i in 1:gs.num_v
        view(dx, gs.v_idx[i]) .= get_vertex_jacobian(jgd, i) * view(z, gs.v_idx[i])
        dx .+= sum(jgd.e_jac_product[i])
    end
end

# functions for callable structs

Base.:*(Jac::NDJacVecOperator1, z::AbstractVector) = jac_vec_prod(Jac, z)

function LinearAlgebra.mul!(dx::AbstractVector, Jac::NDJacVecOperator1, z::AbstractVector)
    jac_vec_prod!(dx, Jac, z)
end

# callable structs

function (Jac::NDJacVecOperator1)(x, p, t) # auch Number bei t?
    update_coefficients!(Jac, x, p, t)
    Jac*x
end

function (Jac::NDJacVecOperator1)(dx, x, p, t) # t::Number
    update_coefficients!(Jac, x, p, t)
    mul!(dx, Jac, x)
end

#x_test = randn(N)
x_test = ones(N)
p_test = nothing
t_test = 1.0
dx_test = ones(N)

@time update_coefficients!(NDJacVecOp, x_test, p_test, t_test) # 5 allocations
call_update_coefficients! = update_coefficients!(NDJacVecOp, x_test, p_test, t_test)

call_callable_struct_1 = NDJacVecOp(x_test, p_test, t_test)

call_callable_struct_2 = NDJacVecOp(dx_test, x_test, p_test, t_test)

@time NDJacVecOp(x_test, p_test, t_test)
@time NDJacVecOp(dx_test, x_test, p_test, t_test)

@allocated NDJacVecOp(x_test, p_test, t_test)
@allocated NDJacVecOp(dx_test, x_test, p_test, t_test)
