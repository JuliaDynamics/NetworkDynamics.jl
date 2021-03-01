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

function auto_jacvec_test(x)
	x=2
end

# Operator überladen
import Base
Base.:*(L::funktion2, x::Array) = auto_jacvec_test(L.x)
x=[1 2 3]
b*x

using ForwardDiff

function attempt(dx, b, x, p, t)
  dx[1] = x[1]^2+x[2]^3-1
  dx[2] = x[1]^4 - x[2]^4 + x[1]*x[2]
  nothing
end


a = [1.0, 1.0]

dx = Array{Float64,2}(undef, 2, 1)
results = zeros(2, 2)

p = 1.0
t = 0.0
b = [0.0, 0.0]

cfg = ForwardDiff.JacobianConfig((y, x) -> attempt(y, b, x, p, t), dx, a)

ForwardDiff.jacobian!(results, (y, x) -> attempt(y, b, x, p, t), dx, a, cfg, Val{false}())




Base.:*(L::JacVecOperator,x::AbstractVector) = L.autodiff ? auto_jacvec(_u->L.f(_u,L.p,L.t),L.u,x) : num_jacvec(_u->L.f(_u,L.p,L.t),L.u,x)

function auto_jacvec!(du, f, x, v,
                 cache1 = ForwardDiff.Dual{JacVecTag}.(x, x), # this won't alias
                 cache2 = similar(cache1))
    cache1 .= ForwardDiff.Dual{JacVecTag}.(x, reshape(v, size(x)))
    f(cache2,cache1)
    du .= vec(ForwardDiff.partials.(cache2, 1))
end

function auto_jacvec(f, x, v)
    vv = reshape(v, axes(x))
    ForwardDiff.partials.(vec(f(ForwardDiff.Dual{JacVecTag}.(x, vv))), 1)
end
