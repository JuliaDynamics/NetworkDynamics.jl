import NetworkDynamics as OldND
include("testnetwork_constructor.jl")

N = 100

Ns = Int[]
at1s = Float64[]
at2s = Float64[]
ato1s = Float64[]
ato2s = Float64[]
t1s = Float64[]
t2s = Float64[]
to1s = Float64[]
to2s = Float64[]


for N in Int[1e1,1e2,1e3,1e4,1e5,1e6,1e7]
    p, v, e, g = heterogeneous(N)
    at1 = @elapsed nd1 = Network(g, v, e;
        execution=SequentialExecution{false}(),
        # execution=ThreadedExecution{true}(),
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

    at2 = @elapsed nd2 = Network(g, v, e; execution=ThreadedExecution{true}(), accumulator=KAAggregator(+))
    ato2 = @elapsed old_nd2 = OldND.network_dynamics(vold, eold, g; parallel=true)

    b1 = @b $nd1($dx, $x0, $p, NaN)

    b2 = @b $nd2($dx, $x0, $p, NaN)

    bo1 = @b $old_nd1($dxold, $x0old, $pold, NaN)
    bo2 = @b $old_nd2($dxold, $x0old, $pold, NaN)
    push!(Ns, N)
    push!(at1s, at1)
    push!(at2s, at2)
    push!(ato1s, ato1)
    push!(ato2s, ato2)
    push!(t1s, b1.time)
    push!(t2s, b2.time)
    push!(to1s, bo1.time)
    push!(to2s, bo2.time)
end

using GLMakie
fig = Figure()
ax = Axis(fig[1,1]; xscale=log10, yscale=log10, ylabel="construction time",xticks=Ns,
title="Benchmark: coreloop for inhomogenious kuramoto network")
scatterlines!(ax, Ns, at1s, label="new, sequential")
scatterlines!(ax, Ns, at2s, label="new, multithreaded")
scatterlines!(ax, Ns, ato1s, label="old, sequential", linestyle=:dash)
scatterlines!(ax, Ns, ato2s, label="old, multithreaded", linestyle=:dash)
axislegend(ax; position=:lt)

ax = Axis(fig[2,1]; xscale=log10, yscale=log10, ylabel="coreloop time",xticks=Ns,xlabel="network size")
scatterlines!(ax, Ns, t1s, label="new, sequential")
scatterlines!(ax, Ns, t2s, label="new, multithreaded")
scatterlines!(ax, Ns, to1s, label="old, sequential", linestyle=:dash)
scatterlines!(ax, Ns, to2s, label="old, multithreaded", linestyle=:dash)

using GLMakie
fig = Figure()
ax = Axis(fig[1,1]; xscale=log10, yscale=log10, ylabel="construction time",xticks=Ns,
title="Benchmark: coreloop for inhomogenious kuramoto network")
scatterlines!(ax, Ns, at1s, label="new, sequential")
scatterlines!(ax, Ns, at2s, label="new, multithreaded")
scatterlines!(ax, Ns, ato1s, label="old, sequential", linestyle=:dash)
scatterlines!(ax, Ns, ato2s, label="old, multithreaded", linestyle=:dash)
axislegend(ax; position=:lt)

ax = Axis(fig[2,1]; xscale=log10, yscale=log10, ylabel="coreloop time",xticks=Ns,xlabel="network size")
scatterlines!(ax, Ns, t1s, label="new, sequential")
scatterlines!(ax, Ns, t2s, label="new, multithreaded")
scatterlines!(ax, Ns, to1s, label="old, sequential", linestyle=:dash)
scatterlines!(ax, Ns, to2s, label="old, multithreaded", linestyle=:dash)
