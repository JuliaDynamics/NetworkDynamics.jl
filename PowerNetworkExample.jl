include("src/NetworkDynamics.jl")
using .NetworkDynamics
using LightGraphs
using LinearAlgebra
using DifferentialEquations
using Plots
pyplot()

g = barabasi_albert(10,5)

begin
    #=
    Conventon is the following:
    v[1] + 1.j*v[2] is the complex voltage at a vertex.
    e[1] + 1.j*e[2] is the complex current at an edge.
    =#
end
struct complex_admittance_edge!
    admittance
end

function (cae::complex_admittance_edge!)(e,v_s,v_d,p,t)
    source_voltage = v_s[1] + v_s[2]*im
    destination_voltage = v_d[1] + v_d[2]*im
    # If current is flowing away from the source, it is negative at the source.
    complex_current = cae.admittance * (destination_voltage - source_voltage)
    e[1] = real(complex_current)
    e[2] = imag(complex_current)
    nothing
end

begin
    # Making sure the function works as intended
    v1 = [1.,2.]
    v2 = [3.,4.]
    e = zeros(10)
    e1 = view(l, 5:6)
    cae = complex_admittance_edge!(3. + 5*im)
    cae(e1,v1,v2,nothing,0.)
    println(e1)
end

function total_current(e_s, e_d)
    # Keeping with the convention of negative sign for outging current
    current = 0.0im
    for e in e_s
        current += e[1] + e[2]*im
    end
    for e in e_d
        current -= e[1] + e[2]*im
    end
    current
end

struct SwingVertex
    P
    D
end
function (sv::SwingVertex)(dv, v, e_s, e_d, p, t)
    current = total_current(e_s, e_d)
    voltage = v[1] + v[2] * im
    dv[3] = sv.P - sv.D * v[3] + real(voltage * conj(current))
    dvolt = 1.0im * v[3] * voltage - (abs(voltage) - 1) * voltage
    dv[1] = real(dvolt)
    dv[2] = imag(dvolt)
    nothing
end


struct PQVertex
    P_complex
end
function (pq::PQVertex)(dv, v, e_s, e_d, p, t)
    current = total_current(e_s, e_d)
    voltage = v[1] + v[2] * im
    residual = pq.P_complex - voltage * conj(current)
    dv[1] = real(residual)
    dv[2] = imag(residual)
    nothing
end

# Example PQ node:
pq_1 = ODEVertex(f! = PQVertex(randn() + randn()*im),
                 dim = 2,
                 massmatrix = 0.)

using GraphPlot
gplot(g)

pq_list = [ODEVertex(f! = PQVertex(randn() + randn()*im),
                     dim = 2,
                     massmatrix = 0.,
                     sym = [:v_r, :v_i])
           for i in 1:5]

vertex_list = [ODEVertex(f! = SwingVertex(randn(), 1.),
                        dim = 3,
                        sym = [:v_r, :v_i, :ω])
              for i in 1:5]

append!(vertex_list, pq_list)

swing_list = [ODEVertex(f! = SwingVertex(randn(), 1.),
                        dim = 3,
                        sym = [:v_r, :v_i, :ω])
              for i in 1:10]


edge_list = [StaticEdge(f! = complex_admittance_edge!(0.0 - 5.0im),
                        dim = 2)
             for e in edges(g)]

power_network_rhs = network_dynamics(vertex_list, edge_list, g)

begin
    x0 = rand(25)
    test_prob = ODEProblem(power_network_rhs,x0,(0.,50.))
end
test_sol = solve(test_prob, Rosenbrock23(autodiff=false), force_dtmin=true)

struct root_rhs
    rhs
    mm
end
function (rr::root_rhs)(x)
    dx = similar(x)
    rr.rhs(dx, x, nothing, 0.)
    rr.mm * dx .- dx
end

rr = root_rhs(power_network_rhs, power_network_rhs.mass_matrix)

using NLsolve

nl_res = nlsolve(rr, x0)
ic = nl_res.zero
test_prob = ODEProblem(power_network_rhs,ic,(0.,50.))
test_sol = solve(test_prob, Rosenbrock23(autodiff=false))

test_sol
plot(test_sol)

plot(test_sol, vars = [s for s in power_network_rhs.syms if occursin("ω", string(s))])

plot(test_sol, vars=[:ω_1])
