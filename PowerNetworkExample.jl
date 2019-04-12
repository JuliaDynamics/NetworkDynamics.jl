include("src/NetworkDynamics.jl")
using .NetworkDynamics
using LightGraphs
using LinearAlgebra
using DifferentialEquations
using Plots

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

ident = diagm(0 => ones(3))

vertex_list = [ODEVertex(SwingVertex(randn(), 1.), 3, ident, [:v_r, :v_i, :ω]) for v in vertices(g)]
edge_list = [StaticEdge(complex_admittance_edge!(0.0 - 5.0im), 2) for e in edges(g)]

power_network_rhs = network_dynamics(vertex_list, edge_list, g)

pyplot()

begin
    x0 = rand(30)
    test_prob = ODEProblem(power_network_rhs,x0,(0.,50.))
    test_sol = solve(test_prob)
    plot(test_sol, vars = :ω_1)
end

plot(test_sol, vars = [s for s in power_network_rhs.syms if occursin("ω", string(s))])

plot(test_sol, vars=[:ω_1])

using GraphPlot
gplot(g)
