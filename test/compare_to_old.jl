using Metal
using Adapt
using Random
using OrdinaryDiffEq
using DiffEqGPU
import NetworkDynamics as OldND
include("testnetwork_constructor.jl")

N = 100

Ns = Int[]
at1s = Float64[]
at2s = Float64[]
at3s = Float64[]
ato1s = Float64[]
ato2s = Float64[]
t1s = Float64[]
t2s = Float64[]
t3s = Float64[]
to1s = Float64[]
to2s = Float64[]
st1s = Float64[]
st2s = Float64[]
st3s = Float64[]


# empty!(Ns); empty!(at1s); empty!(at2s); empty!(at3s); empty!(ato1s); empty!(ato2s); empty!(t1s); empty!(t2s); empty!(t3s); empty!(to1s); empty!(to2s); empty!(st1s); empty!(st2s); empty!(st3s);
for N in Int[1e1,1e2,1e3,1e4,1e5,1e6]
    p, v, e, g = heterogeneous(N)
    at1 = @elapsed nd1 = Network(g, v, e;
        execution=SequentialExecution{true}(),
        # execution=KAExecution{true}(),
        # accumulator=NNlibScatter(+)
        accumulator=SequentialAggregator(+)
        # accumulator=KAAggregator(+)
        # accumulator=NaiveAggregator(+)
    );
    x0old = rand(dim(nd1))

    pold, vold, eold, gold = heterogeneous_old(N)
    ato1 = @elapsed old_nd1 = OldND.network_dynamics(vold, eold, g)
    old_vidx = old_nd1(OldND.GetGS).v_idx

    x0 = zeros(length(x0old))
    for i in 1:nv(g)
        x0[nd1.im.v_data[i]] .= x0old[old_vidx[i]]
    end
    @test sort(x0)== sort(x0old)

    dx = similar(x0)
    dxold = similar(x0)

    old_nd1(dxold, x0old, pold, NaN)
    nd1(dx, x0, p, NaN)

    @testset "du equal?" begin
        for i in 1:nv(g)
            @test dx[nd1.im.v_data[i]] ≈ dxold[old_vidx[i]]
        end
    end

    at2 = @elapsed nd2 = Network(g, v, e; execution=KAExecution{true}(), accumulator=PolyesterAggregator(+))
    # at2 = @elapsed nd2 = Network(g, v, e; execution=KAExecution{true}(), accumulator=KAAggregator(+))
    ato2 = @elapsed old_nd2 = OldND.network_dynamics(vold, eold, g; parallel=true)

    backend = MetalBackend()
    at3 = @elapsed nd3 = let
        _nd = Network(g, v, e; execution=KAExecution{true}(), accumulator=KAAggregator(+))
        adapt(backend, _nd)
    end
    dx_d = adapt(backend, convert.(Float32, dx))
    x0_d = adapt(backend, convert.(Float32, x0))
    p_d  = adapt(backend, convert.(Float32, p))

    b1 = @b $nd1($dx, $x0, $p, NaN)
    b2 = @b $nd2($dx, $x0, $p, NaN)
    b3 = @b $nd3($dx_d, $x0_d, $p_d, NaN)

    bo1 = @b $old_nd1($dxold, $x0old, $pold, NaN)
    bo2 = @b $old_nd2($dxold, $x0old, $pold, NaN)

    prob = ODEProblem(nd1, x0, (0,1),p)
    sol = solve(prob, Tsit5())
    st1 = @elapsed sol = solve(prob, Tsit5())
    prob = ODEProblem(nd2, x0, (0,1),p)
    sol = solve(prob, Tsit5())
    st2 = @elapsed sol = solve(prob, Tsit5())
    prob = ODEProblem(nd3, x0_d, (Float32(0),Float32(1)),p_d)
    sol = solve(prob, Tsit5())
    st3 = @elapsed sol = solve(prob, Tsit5())

    push!(Ns, N)
    push!(at1s, at1)
    push!(at2s, at2)
    push!(at3s, at3)
    push!(ato1s, ato1)
    push!(ato2s, ato2)
    push!(t1s, b1.time)
    push!(t2s, b2.time)
    push!(t3s, b3.time)
    push!(to1s, bo1.time)
    push!(to2s, bo2.time)

    push!(st1s, st1)
    push!(st2s, st2)
    push!(st3s, st3)
end

# ax = Axis(fig[1,1]; xscale=log10, yscale=log10, ylabel="construction time",xticks=Ns)
# ax = Axis(fig[1,1]; xscale=log10, yscale=log10, ylabel="construction time", title="Benchmark: coreloop for inhomogenious kuramoto network")
# scatterlines!(ax, Ns, at1s, label="new, sequential")
# scatterlines!(ax, Ns, at2s, label="new, multithreaded")
# scatterlines!(ax, Ns, at3s, label="new, GPU")
# scatterlines!(ax, Ns, ato1s, label="old, sequential", linestyle=:dash)
# scatterlines!(ax, Ns, ato2s, label="old, multithreaded", linestyle=:dash)
# axislegend(ax; position=:lt)

using GLMakie
fig = Figure()
ax = Axis(fig[1,1]; xscale=log10, yscale=log10, ylabel="coreloop time")
scatterlines!(ax, Ns, t1s, label="new, sequential")
scatterlines!(ax, Ns, t2s, label="new, multithreaded")
scatterlines!(ax, Ns, t3s, label="new, GPU")
scatterlines!(ax, Ns, to1s, label="old, sequential", linestyle=:dash)
scatterlines!(ax, Ns, to2s, label="old, multithreaded", linestyle=:dash)
axislegend(ax; position=:lt)

ax = Axis(fig[2,1]; xscale=log10, yscale=log10, ylabel="solve time",xlabel="network size")
scatterlines!(ax, Ns, st1s, label="new, sequential")
scatterlines!(ax, Ns, st2s, label="new, multithreaded")
scatterlines!(ax, Ns, st3s, label="new, GPU")
