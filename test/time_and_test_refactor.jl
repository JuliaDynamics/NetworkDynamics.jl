using Pkg
Pkg.activate(".")
using BenchmarkTools
using NetworkDynamics
using Test
using LightGraphs
using Random


for N = [20,200]
    # Pure Julia "ground truth"

    g = watts_strogatz(N, 5, 0.)
    println("\nTesting ND with $N nodes and $(ne(g)) edges:\n")

    B = incidence_matrix(g, oriented=true)

    struct kuramoto_dyn{T, T2, U}
        B::T
        B_trans::T2
        ω::U
        N::Int
    end
    function (dd::kuramoto_dyn)(dx, x, p, t)
        dx[1:N] .= x[N+1:2N]
        dx[N+1:2N] .= dd.ω .- x[N+1:2N] .- 5. * dd.B * sin.(dd.B_trans * x[1:N])
        nothing
    end

    ω = Array{Float64}(1:N) ./ N
    kn! = kuramoto_dyn(B, transpose(B), ω, N)

    # NetworkDynamics

    function kuramoto_edge!(e, θ_s, θ_d, K, t)
        e[1] = K * sin(θ_s[1] - θ_d[1])
    end

    function kuramoto_inertia!(dv, v, e_s, e_d, ω, t)
        dv[1] = v[2]
        dv[2] = ω - v[2]
        for e in e_d
            dv[2] += e[1]
        end
    end

    inertia!  = ODEVertex(f! = kuramoto_inertia!, dim = 2, sym= [:θ, :ω]);
    edge!     = StaticEdge(f! = kuramoto_edge!, dim = 1)
    edge_anti! = StaticEdge(f! = kuramoto_edge!, dim = 1, coupling = :antisymmetric )



    nd! = network_dynamics(inertia!, edge!, g)
    nd_anti! = network_dynamics(inertia!, edge_anti!, g)
    nd_dir! = network_dynamics(inertia!, edge!, SimpleDiGraph(g))

    # Testing accuracy

    @testset "Accuracy" begin
        for i in 1:1
            continue
            x0_θ = randn(N)
            x0_ω = randn(N)

            x0_kn = [x0_θ; x0_ω ]

            x0_nd = Vector(vec([x0_θ  x0_ω]'))
            x0_nd_anti = Vector(vec([x0_θ  x0_ω]'))
            x0_nd_dir = Vector(vec([x0_θ  x0_ω]'))

            dx_kn = zeros(2N)
            dx_nd = zeros(2N)
            dx_nd_anti = zeros(2N)
            dx_nd_dir = zeros(2N)

            kn!(dx_kn, x0_kn, nothing, 0.)
            nd!(dx_nd, x0_nd, (ω, 5.), 0.)
            nd_anti!(dx_nd_anti, x0_nd_anti, (ω, 5.), 0.)
            nd_dir!(dx_nd_dir, x0_nd_dir, (ω, 5.), 0.)

            @test isapprox(dx_kn, [dx_nd[1:2:2N]; dx_nd[2:2:2N]])
            @test isapprox(dx_kn, [dx_nd_anti[1:2:2N]; dx_nd_anti[2:2:2N]])
            @test isapprox(dx_kn, [dx_nd_dir[1:2:2N]; dx_nd_dir[2:2:2N]])

        end
    end
    println("")

    dx0   = zeros(2N)
    x0 = randn(2N)
    Random.seed!(42)
    display(@benchmark($nd_anti!($dx0, $x0, ($ω, 5.), 0.)))
end


# For N = [20,200] ND#master is at 1.3 μs, 10.5 μs (median time) on my machine


N = 2

# Pure Julia "ground truth"

g = watts_strogatz(N, Int(round(N/2)), 0.)
println("\nTesting ND with $N nodes and $(ne(g)) edges:\n")


g=SimpleGraph(N,1)
B = incidence_matrix(g, oriented=true)

struct kuramoto_dyn{T, T2, U}
    B::T
    B_trans::T2
    ω::U
    N::Int
end
function (dd::kuramoto_dyn)(dx, x, p, t)
    dx[1:N] .= x[N+1:2N]
    dx[N+1:2N] .= dd.ω .- x[N+1:2N] .- 5. * dd.B * sin.(dd.B_trans * x[1:N])
    nothing
end

ω = Array{Float64}(1:N) ./ N
kn! = kuramoto_dyn(B, transpose(B), ω, N)

# NetworkDynamics

function kuramoto_edge!(e, θ_s, θ_d, K, t)
    e[1] = K * (θ_s[1] - θ_d[1])
end

function kuramoto_inertia!(dv, v, e_s, e_d, ω, t)
    dv[1] = v[2]
    dv[2] = ω - v[2]
    #print(e_d)
    for e in e_d
        dv[2] += e[1]
    end
end

inertia!  = ODEVertex(f! = kuramoto_inertia!, dim = 2, sym= [:θ, :ω]);
edge!     = StaticEdge(f! = kuramoto_edge!, dim = 1)
edge_anti! = StaticEdge(f! = kuramoto_edge!, dim = 1, coupling = :antisymmetric )



nd! = network_dynamics(inertia!, edge!, g)
nd_anti! = network_dynamics(inertia!, edge_anti!, g)
nd_dir! = network_dynamics(inertia!, edge!, SimpleDiGraph(g))

x0_θ = [1.0,2.0]
x0_ω = [1.0,1.0]

x0_kn = [x0_θ; x0_ω ]

x0_nd = Vector(vec([x0_θ  x0_ω]'))
x0_nd_anti = Vector(vec([x0_θ  x0_ω]'))
x0_nd_dir = Vector(vec([x0_θ  x0_ω]'))

dx_kn = zeros(2N)
dx_nd = zeros(2N)
dx_nd_anti = zeros(2N)
dx_nd_dir = zeros(2N)

kn!(dx_kn, x0_kn, nothing, 0.)
@btime $nd!($dx_nd, $x0_nd, ($ω, 5.), 0.)
@btime $nd_anti!($dx_nd_anti, $x0_nd_anti, ($ω, 5.), 0.)
@btime $nd_dir!($dx_nd_dir, $x0_nd_dir, ($ω, 5.), 0.)


dx_nd
dx_nd_anti
dx_nd_dir



dx_nd .- dx_nd_dir
