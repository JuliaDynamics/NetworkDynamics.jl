using NetworkDynamics
using OrdinaryDiffEq
using LightGraphs

N = 5
k = 2
g = barabasi_albert(N,k)
t_max = 10
x0 = [rand() for i in 1:N]

p_vertex = 1.0


@inline function kuramoto_vertex!(dv, v, edges, p, t)
	dv[1] = p
	@inbounds for e in edges
			dv[1] += e[1]
	end
	nothing
end

@inline function kuramoto_edge!(e,v_s,v_d,p,t)
    e[1] = 1.0 * sin(v_s[1] - v_d[1]) # 1.0 statt p
    nothing
end

odevertex = ODEVertex(f! = kuramoto_vertex!, dim = 1)
staticedge = StaticEdge(f! = kuramoto_edge!, dim = 1)

parameters = (p_vertex, nothing)


kuramoto_network_ND! = network_dynamics(odevertex, staticedge, g)

kur_prob_ND = ODEProblem(kuramoto_network_ND!, x0, (0.,t_max), parameters)

kur_sol_ND=solve(kur_prob_ND, DP5())



### now with Jacobians

function vertex_jacobian(J_intern, z, v, node_index)
	J_intern = 0 # internal jacobian of vertex i
	return v = J_intern * z[node_index]
	#nothing
end

internal_jacobian = zeros(Int64, 1) # set the internal Jacobian (firstly) to zero
z = [1,0,0,0,0] # directional derivative at location of vertex 1
v = 0.0
v = vertex_jacobian(internal_jacobian, z, v, 1) # jacobian of vertex 1

function vertex_edge(J_edge, z, e, node_index)
	J_edge[1,1] = 0
	J_edge[1,2] = 1



struct test_struct
	field_function1
	field_function2
end


# callable struct
function(structi::test_struct)(x,y)
	sum = structi.field_function1(x,y) + structi.field_function2(x,y)
	return sum
end

field_function1(x,y) = x*x
field_function2(x,y) = y*y

function(f1::test_struct)(x,y)
	f1.field_function1(x,y)

test = test_struct(field_function1, field_function2)

struct final_struct
	parameter1
	parameter2
end

function(final::final_struct)(x,y)
	sum = final.parameter1 + final.parameter2
	return sum
end

final_object = final_struct(1,2)
final_object(test(1,2))

# Fieldname function1 wird ausgeführt
test.field_function2(1,2)

test.field_function2 = test.field_function2(1,2)

# call callable struct
test(1,2)

print(test)

function1(x,y) = cos(x-y)

test1 = function1(1,3)


struct Jac_vertex
    Jac_intern
end

struct Jac_vertex2
	Jac_intern2
end

#Jac_intern(a,b) = a+b

function Jac_intern2(x, p, t)
	a = x+p
    return a # Funtion von DiffEq.jl, die den Jacobian (automatisch) berechnet
end

j2 = Jac_vertex2(Jac_intern2)

Jac_intern2(1,2,3)


## Test return nothing
mutable struct funktion2
    x
end

b=funktion2(3)

function f1!(b::funktion2)
    b.x=f2(b.x)
    println(b.x)
    nothing
end

function f2(y)
    y=2
	return y
end

f1!(b)
f1!(3)




using ForwardDiff

function attempt(dx, x, p, t)
  dx[1] = x[1]^2+x[2]^3-1
  dx[2] = x[1]^4 - x[2]^4 + x[1]*x[2]
  nothing
end


a = [2.0, 1.0]

dx = Array{Float64,2}(undef, 2, 1)
results = zeros(2, 2)

p = 1.0
t = 0.0


cfg = ForwardDiff.JacobianConfig((y, x) -> attempt(y, x, p, t), dx, a)

ForwardDiff.jacobian!(results, (y, x) -> attempt(y, x, p, t), dx, a, cfg, Val{false}())

results


@Base.kwdef struct struct_with_function
	dimension::Int64
	funktionsname = :F
end

function funktion()
	print("Hallo")
end

test_objekt = struct_with_function(dimension = 1)

test_objekt_with_function = struct_with_function(dimension = 1, funktionsname = funktion)

typeof(test_objekt.funktionsname)
typeof(test_objekt_with_function.funktionsname)

test_objekt.funktionsname = funktion

using LinearAlgebra

typeof(I) <: AbstractArray
################ arrays sizes ##################################################

v_dims = 2
e_dims = 2
num_e = 10

v_jac_array = [Array{Float64,2}(undef, dim, dim) for dim in v_dims]

println(v_jac_array)


# For edges there is another array layer in the buffer since each edge has two Jacobians
#e_jac_array = [[zeros(dim, srcdim), zeros(dim, dstdim)] for (dim, srcdim, dstdim) in zip(gs.e_dims, v_src_dims, v_dst_dims)]
e_jac_array = [[zeros(dim, srcdim), zeros(dim, dstdim)] for (dim, srcdim, dstdim) in zip(e_dims, v_dims, v_dims)]
e_jac_product = [zeros(e_dims)]

println(e_jac_array)
typeof(e_jac_array)

example_v_jac_array = [1.0 3.0; 5.0 7.0] # siehe Konsole

example_e_jac_array = [[1.0 3.0; 5.0 7.0], [2.0 4.0; 6.0 8.0]] # siehe Konsole

example_e_jac_product = [[1.0, 2.0]] # siehe Konsole














filled_array[1][1]

(vcat(e_jac_array, filled_array))[1][2]


e_jac_array_num_e = [[zeros(dim, srcdim), zeros(dim, dstdim)] for (dim, srcdim, dstdim) in zip(e_dims, v_dims, v_dims)]


get_src_edge_jacobian(i::Int) = e_jac_array[i, 1]

get_src_edge_jacobian_micha(i::Int) = e_jac_array[i][1] ## Micha: gd.e_jac_array[i][1]

get_dst_edge_jacobian(i::Int) = e_jac_array[i, 2]


get_src_edge_jacobian(1)

get_dst_edge_jacobian(1)


e_jac_array1 = [ [[1 3 5], [2 4 6]], [[7 8 9], [10 11 12]] ]

############

e_jac_array_part1 = [[zeros(dim, srcdim), zeros(dim, dstdim)] for (dim, srcdim, dstdim) in zip(e_dims, v_dims, v_dims)]


e_jac_array_part2 = [[zeros(dim, srcdim), zeros(dim, dstdim)] for (dim, srcdim, dstdim) in zip(e_dims, v_dims, v_dims)]

e_jac_array_part3 = [[zeros(dim, srcdim), zeros(dim, dstdim)] for (dim, srcdim, dstdim) in zip(e_dims, v_dims, v_dims)]

A = vcat(e_jac_array_part1, e_jac_array_part2, e_jac_array_part3)

A[3][2]



get_src_edge_jacobian_of_parts(i::Int) = A[i][1]
get_dst_edge_jacobian_of_parts(i::Int) = A[i][2]

get_src_edge_jacobian_of_parts(3)
get_dst_edge_jacobian_of_parts(3)

z_edge_1 = [1 3; 2 4]
z_edge_2 = [1 2; 3 4]
z_edge_3 = [1 2; 3 4]
z = vcat(z_edge_1, z_edge_2, z_edge_3)
z[1, :]


e_jac_prod_1 = get_src_edge_jacobian_of_parts(1) * z[1, :]
e_jac_prod_2 = get_src_edge_jacobian_of_parts(2) * z[2, :]

e_jac_prod_sum = hcat(e_jac_prod_1, e_jac_prod_2)


e_jac_prod_sum = Array{Array{Float64,2},1}
e_jac_prod_sum = Array{Float64,2}
e_jac_prod_sum = zeros(2, 2)
for i in 1:2
	#print(size(get_src_edge_jacobian_of_parts(i) * z[i, :]))
	e_jac_prod_sum[i, :] = get_src_edge_jacobian_of_parts(i) * z[i, :]
end

e_jac_prod_sum

#### other struct stuff

abstract type VertexFunction end
"""
Abstract supertype for all edge functions.
"""
abstract type EdgeFunction end

@Base.kwdef struct ODEVertex{T} <: VertexFunction
    f!::T # signature (dx, x, edges, p, t) -> nothing
    dim::Int
	mass_matrix = I
    #vertex_jacobian!::F # signature (J::AbstractMatrix, v, p, t)
end

@Base.kwdef struct ODEVertex1{T} <: VertexFunction
    #f!::T # signature (dx, x, edges, p, t) -> nothing
    dim::Int64
	mass_matrix = I
    #vertex_jacobian!::F # signature (J::AbstractMatrix, v, p, t)
end


object = ODEVertex1(dim = 1, mass_matrix = nothing)

@Base.kwdef struct ODEVertex3{T}
    #f!::T # signature (dx, x, edges, p, t) -> nothing
    dim::Int64
	#mass_matrix = I
    #vertex_jacobian!::F # signature (J::AbstractMatrix, v, p, t)
end

object2 = ODEVertex3(2)

@Base.kwdef struct hallo
	zahl::Int64
	#masse = 1.0
end

test = hallo(zahl = 3)

@Base.kwdef struct hallo1
	zahl::Int64
	masse = 1.0
end

test1 = hallo1(zahl = 2, masse = 2.0)
test2 = hallo1(zahl = 100)

@Base.kwdef struct struct_with_function
	dimension::Int64
	funktionsname = :F
end

function funktion()
	print("Hallo")
end

test_objekt = struct_with_function(dimension = 1)

test_objekt_with_function = struct_with_function(dimension = 1, funktionsname = funktion)
